//===------ PPCGCodeGeneration.cpp - Polly Accelerator Code Generation. ---===//
//
//                     The LLVM Compiler Infrastructure
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//
//
// Take a scop created by ScopInfo and map it to GPU code using the ppcg
// GPU mapping strategy.
//
//===----------------------------------------------------------------------===//

#include "polly/CodeGen/IslNodeBuilder.h"
#include "polly/CodeGen/Utils.h"
#include "polly/DependenceInfo.h"
#include "polly/LinkAllPasses.h"
#include "polly/Options.h"
#include "polly/ScopInfo.h"
#include "llvm/Analysis/AliasAnalysis.h"
#include "llvm/Analysis/BasicAliasAnalysis.h"
#include "llvm/Analysis/GlobalsModRef.h"
#include "llvm/Analysis/PostDominators.h"
#include "llvm/Analysis/ScalarEvolutionAliasAnalysis.h"

#include "isl/union_map.h"

extern "C" {
#include "ppcg/cuda.h"
#include "ppcg/gpu.h"
#include "ppcg/gpu_print.h"
#include "ppcg/ppcg.h"
#include "ppcg/schedule.h"
}

#include "llvm/Support/Debug.h"

using namespace polly;
using namespace llvm;

#define DEBUG_TYPE "polly-codegen-ppcg"

static cl::opt<bool> DumpSchedule("polly-acc-dump-schedule",
                                  cl::desc("Dump the computed GPU Schedule"),
                                  cl::Hidden, cl::init(false), cl::ZeroOrMore,
                                  cl::cat(PollyCategory));

static cl::opt<bool>
    DumpCode("polly-acc-dump-code",
             cl::desc("Dump C code describing the GPU mapping"), cl::Hidden,
             cl::init(false), cl::ZeroOrMore, cl::cat(PollyCategory));

/// Create the ast expressions for a ScopStmt.
///
/// This function is a callback for to generate the ast expressions for each
/// of the scheduled ScopStmts.
static __isl_give isl_id_to_ast_expr *pollyBuildAstExprForStmt(
    void *Stmt, isl_ast_build *Build,
    isl_multi_pw_aff *(*FunctionIndex)(__isl_take isl_multi_pw_aff *MPA,
                                       isl_id *Id, void *User),
    void *UserIndex,
    isl_ast_expr *(*FunctionExpr)(isl_ast_expr *Expr, isl_id *Id, void *User),
    void *User_expr) {

  // TODO: Implement the AST expression generation. For now we just return a
  // nullptr to ensure that we do not free uninitialized pointers.

  return nullptr;
}

/// Generate code for a GPU specific isl AST.
///
/// The GPUNodeBuilder augments the general existing IslNodeBuilder, which
/// generates code for general-prupose AST nodes, with special functionality
/// for generating GPU specific user nodes.
///
/// @see GPUNodeBuilder::createUser
class GPUNodeBuilder : public IslNodeBuilder {
public:
  GPUNodeBuilder(PollyIRBuilder &Builder, ScopAnnotator &Annotator, Pass *P,
                 const DataLayout &DL, LoopInfo &LI, ScalarEvolution &SE,
                 DominatorTree &DT, Scop &S)
      : IslNodeBuilder(Builder, Annotator, P, DL, LI, SE, DT, S) {}

private:
  /// Create code for user-defined AST nodes.
  ///
  /// These AST nodes can be of type:
  ///
  ///   - ScopStmt:      A computational statement (TODO)
  ///   - Kernel:        A GPU kernel call (TODO)
  ///   - Data-Transfer: A GPU <-> CPU data-transfer (TODO)
  ///
  /// @param UserStmt The ast node to generate code for.
  virtual void createUser(__isl_take isl_ast_node *UserStmt);
};

void GPUNodeBuilder::createUser(__isl_take isl_ast_node *UserStmt) {
  isl_ast_node_free(UserStmt);
  return;
}

namespace {
class PPCGCodeGeneration : public ScopPass {
public:
  static char ID;

  /// The scop that is currently processed.
  Scop *S;

  LoopInfo *LI;
  DominatorTree *DT;
  ScalarEvolution *SE;
  const DataLayout *DL;
  RegionInfo *RI;

  PPCGCodeGeneration() : ScopPass(ID) {}

  /// Construct compilation options for PPCG.
  ///
  /// @returns The compilation options.
  ppcg_options *createPPCGOptions() {
    auto DebugOptions =
        (ppcg_debug_options *)malloc(sizeof(ppcg_debug_options));
    auto Options = (ppcg_options *)malloc(sizeof(ppcg_options));

    DebugOptions->dump_schedule_constraints = false;
    DebugOptions->dump_schedule = false;
    DebugOptions->dump_final_schedule = false;
    DebugOptions->dump_sizes = false;

    Options->debug = DebugOptions;

    Options->reschedule = true;
    Options->scale_tile_loops = false;
    Options->wrap = false;

    Options->non_negative_parameters = false;
    Options->ctx = nullptr;
    Options->sizes = nullptr;

    Options->tile_size = 32;

    Options->use_private_memory = false;
    Options->use_shared_memory = false;
    Options->max_shared_memory = 0;

    Options->target = PPCG_TARGET_CUDA;
    Options->openmp = false;
    Options->linearize_device_arrays = true;
    Options->live_range_reordering = false;

    Options->opencl_compiler_options = nullptr;
    Options->opencl_use_gpu = false;
    Options->opencl_n_include_file = 0;
    Options->opencl_include_files = nullptr;
    Options->opencl_print_kernel_types = false;
    Options->opencl_embed_kernel_code = false;

    Options->save_schedule_file = nullptr;
    Options->load_schedule_file = nullptr;

    return Options;
  }

  /// Get a tagged access relation containing all accesses of type @p AccessTy.
  ///
  /// Instead of a normal access of the form:
  ///
  ///   Stmt[i,j,k] -> Array[f_0(i,j,k), f_1(i,j,k)]
  ///
  /// a tagged access has the form
  ///
  ///   [Stmt[i,j,k] -> id[]] -> Array[f_0(i,j,k), f_1(i,j,k)]
  ///
  /// where 'id' is an additional space that references the memory access that
  /// triggered the access.
  ///
  /// @param AccessTy The type of the memory accesses to collect.
  ///
  /// @return The relation describing all tagged memory accesses.
  isl_union_map *getTaggedAccesses(enum MemoryAccess::AccessType AccessTy) {
    isl_union_map *Accesses = isl_union_map_empty(S->getParamSpace());

    for (auto &Stmt : *S)
      for (auto &Acc : Stmt)
        if (Acc->getType() == AccessTy) {
          isl_map *Relation = Acc->getAccessRelation();
          Relation = isl_map_intersect_domain(Relation, Stmt.getDomain());

          isl_space *Space = isl_map_get_space(Relation);
          Space = isl_space_range(Space);
          Space = isl_space_from_range(Space);
          Space = isl_space_set_tuple_id(Space, isl_dim_in, Acc->getId());
          isl_map *Universe = isl_map_universe(Space);
          Relation = isl_map_domain_product(Relation, Universe);
          Accesses = isl_union_map_add_map(Accesses, Relation);
        }

    return Accesses;
  }

  /// Get the set of all read accesses, tagged with the access id.
  ///
  /// @see getTaggedAccesses
  isl_union_map *getTaggedReads() {
    return getTaggedAccesses(MemoryAccess::READ);
  }

  /// Get the set of all may (and must) accesses, tagged with the access id.
  ///
  /// @see getTaggedAccesses
  isl_union_map *getTaggedMayWrites() {
    return isl_union_map_union(getTaggedAccesses(MemoryAccess::MAY_WRITE),
                               getTaggedAccesses(MemoryAccess::MUST_WRITE));
  }

  /// Get the set of all must accesses, tagged with the access id.
  ///
  /// @see getTaggedAccesses
  isl_union_map *getTaggedMustWrites() {
    return getTaggedAccesses(MemoryAccess::MUST_WRITE);
  }

  /// Collect parameter and array names as isl_ids.
  ///
  /// To reason about the different parameters and arrays used, ppcg requires
  /// a list of all isl_ids in use. As PPCG traditionally performs
  /// source-to-source compilation each of these isl_ids is mapped to the
  /// expression that represents it. As we do not have a corresponding
  /// expression in Polly, we just map each id to a 'zero' expression to match
  /// the data format that ppcg expects.
  ///
  /// @returns Retun a map from collected ids to 'zero' ast expressions.
  __isl_give isl_id_to_ast_expr *getNames() {
    auto *Names = isl_id_to_ast_expr_alloc(
        S->getIslCtx(),
        S->getNumParams() + std::distance(S->array_begin(), S->array_end()));
    auto *Zero = isl_ast_expr_from_val(isl_val_zero(S->getIslCtx()));
    auto *Space = S->getParamSpace();

    for (int I = 0, E = S->getNumParams(); I < E; ++I) {
      isl_id *Id = isl_space_get_dim_id(Space, isl_dim_param, I);
      Names = isl_id_to_ast_expr_set(Names, Id, isl_ast_expr_copy(Zero));
    }

    for (auto &Array : S->arrays()) {
      auto Id = Array.second->getBasePtrId();
      Names = isl_id_to_ast_expr_set(Names, Id, isl_ast_expr_copy(Zero));
    }

    isl_space_free(Space);
    isl_ast_expr_free(Zero);

    return Names;
  }

  /// Create a new PPCG scop from the current scop.
  ///
  /// The PPCG scop is initialized with data from the current polly::Scop. From
  /// this initial data, the data-dependences in the PPCG scop are initialized.
  /// We do not use Polly's dependence analysis for now, to ensure we match
  /// the PPCG default behaviour more closely.
  ///
  /// @returns A new ppcg scop.
  ppcg_scop *createPPCGScop() {
    auto PPCGScop = (ppcg_scop *)malloc(sizeof(ppcg_scop));

    PPCGScop->options = createPPCGOptions();

    PPCGScop->start = 0;
    PPCGScop->end = 0;

    PPCGScop->context = S->getContext();
    PPCGScop->domain = S->getDomains();
    PPCGScop->call = nullptr;
    PPCGScop->tagged_reads = getTaggedReads();
    PPCGScop->reads = S->getReads();
    PPCGScop->live_in = nullptr;
    PPCGScop->tagged_may_writes = getTaggedMayWrites();
    PPCGScop->may_writes = S->getWrites();
    PPCGScop->tagged_must_writes = getTaggedMustWrites();
    PPCGScop->must_writes = S->getMustWrites();
    PPCGScop->live_out = nullptr;
    PPCGScop->tagged_must_kills = isl_union_map_empty(S->getParamSpace());
    PPCGScop->tagger = nullptr;

    PPCGScop->independence = nullptr;
    PPCGScop->dep_flow = nullptr;
    PPCGScop->tagged_dep_flow = nullptr;
    PPCGScop->dep_false = nullptr;
    PPCGScop->dep_forced = nullptr;
    PPCGScop->dep_order = nullptr;
    PPCGScop->tagged_dep_order = nullptr;

    PPCGScop->schedule = S->getScheduleTree();
    PPCGScop->names = getNames();

    PPCGScop->pet = nullptr;

    compute_tagger(PPCGScop);
    compute_dependences(PPCGScop);

    return PPCGScop;
  }

  /// Collect the array acesses in a statement.
  ///
  /// @param Stmt The statement for which to collect the accesses.
  ///
  /// @returns A list of array accesses.
  gpu_stmt_access *getStmtAccesses(ScopStmt &Stmt) {
    gpu_stmt_access *Accesses = nullptr;

    for (MemoryAccess *Acc : Stmt) {
      auto Access = isl_alloc_type(S->getIslCtx(), struct gpu_stmt_access);
      Access->read = Acc->isRead();
      Access->write = Acc->isWrite();
      Access->access = Acc->getAccessRelation();
      isl_space *Space = isl_map_get_space(Access->access);
      Space = isl_space_range(Space);
      Space = isl_space_from_range(Space);
      Space = isl_space_set_tuple_id(Space, isl_dim_in, Acc->getId());
      isl_map *Universe = isl_map_universe(Space);
      Access->tagged_access =
          isl_map_domain_product(Acc->getAccessRelation(), Universe);
      Access->exact_write = Acc->isWrite();
      Access->ref_id = Acc->getId();
      Access->next = Accesses;
      Accesses = Access;
    }

    return Accesses;
  }

  /// Collect the list of GPU statements.
  ///
  /// Each statement has an id, a pointer to the underlying data structure,
  /// as well as a list with all memory accesses.
  ///
  /// TODO: Initialize the list of memory accesses.
  ///
  /// @returns A linked-list of statements.
  gpu_stmt *getStatements() {
    gpu_stmt *Stmts = isl_calloc_array(S->getIslCtx(), struct gpu_stmt,
                                       std::distance(S->begin(), S->end()));

    int i = 0;
    for (auto &Stmt : *S) {
      gpu_stmt *GPUStmt = &Stmts[i];

      GPUStmt->id = Stmt.getDomainId();

      // We use the pet stmt pointer to keep track of the Polly statements.
      GPUStmt->stmt = (pet_stmt *)&Stmt;
      GPUStmt->accesses = getStmtAccesses(Stmt);
      i++;
    }

    return Stmts;
  }

  /// Derive the extent of an array.
  ///
  /// The extent of an array is defined by the set of memory locations for
  /// which a memory access in the iteration domain exists.
  ///
  /// @param Array The array to derive the extent for.
  ///
  /// @returns An isl_set describing the extent of the array.
  __isl_give isl_set *getExtent(ScopArrayInfo *Array) {
    isl_union_map *Accesses = S->getAccesses();
    Accesses = isl_union_map_intersect_domain(Accesses, S->getDomains());
    isl_union_set *AccessUSet = isl_union_map_range(Accesses);
    isl_set *AccessSet =
        isl_union_set_extract_set(AccessUSet, Array->getSpace());
    isl_union_set_free(AccessUSet);

    return AccessSet;
  }

  /// Derive the bounds of an array.
  ///
  /// For the first dimension we derive the bound of the array from the extent
  /// of this dimension. For inner dimensions we obtain their size directly from
  /// ScopArrayInfo.
  ///
  /// @param PPCGArray The array to compute bounds for.
  /// @param Array The polly array from which to take the information.
  void setArrayBounds(gpu_array_info &PPCGArray, ScopArrayInfo *Array) {
    if (PPCGArray.n_index > 0) {
      isl_set *Dom = isl_set_copy(PPCGArray.extent);
      Dom = isl_set_project_out(Dom, isl_dim_set, 1, PPCGArray.n_index - 1);
      isl_pw_aff *Bound = isl_set_dim_max(isl_set_copy(Dom), 0);
      isl_set_free(Dom);
      Dom = isl_pw_aff_domain(isl_pw_aff_copy(Bound));
      isl_local_space *LS = isl_local_space_from_space(isl_set_get_space(Dom));
      isl_aff *One = isl_aff_zero_on_domain(LS);
      One = isl_aff_add_constant_si(One, 1);
      Bound = isl_pw_aff_add(Bound, isl_pw_aff_alloc(Dom, One));
      Bound = isl_pw_aff_gist(Bound, S->getContext());
      PPCGArray.bound[0] = Bound;
    }

    for (unsigned i = 1; i < PPCGArray.n_index; ++i) {
      isl_pw_aff *Bound = Array->getDimensionSizePw(i);
      auto LS = isl_pw_aff_get_domain_space(Bound);
      auto Aff = isl_multi_aff_zero(LS);
      Bound = isl_pw_aff_pullback_multi_aff(Bound, Aff);
      PPCGArray.bound[i] = Bound;
    }
  }

  /// Create the arrays for @p PPCGProg.
  ///
  /// @param PPCGProg The program to compute the arrays for.
  void createArrays(gpu_prog *PPCGProg) {
    int i = 0;
    for (auto &Element : S->arrays()) {
      ScopArrayInfo *Array = Element.second.get();

      std::string TypeName;
      raw_string_ostream OS(TypeName);

      OS << *Array->getElementType();
      TypeName = OS.str();

      gpu_array_info &PPCGArray = PPCGProg->array[i];

      PPCGArray.space = Array->getSpace();
      PPCGArray.type = strdup(TypeName.c_str());
      PPCGArray.size = Array->getElementType()->getPrimitiveSizeInBits() / 8;
      PPCGArray.name = strdup(Array->getName().c_str());
      PPCGArray.extent = nullptr;
      PPCGArray.n_index = Array->getNumberOfDimensions();
      PPCGArray.bound =
          isl_alloc_array(S->getIslCtx(), isl_pw_aff *, PPCGArray.n_index);
      PPCGArray.extent = getExtent(Array);
      PPCGArray.n_ref = 0;
      PPCGArray.refs = nullptr;
      PPCGArray.accessed = true;
      PPCGArray.read_only_scalar = false;
      PPCGArray.has_compound_element = false;
      PPCGArray.local = false;
      PPCGArray.declare_local = false;
      PPCGArray.global = false;
      PPCGArray.linearize = false;
      PPCGArray.dep_order = nullptr;

      setArrayBounds(PPCGArray, Array);
      i++;

      collect_references(PPCGProg, &PPCGArray);
    }
  }

  /// Create an identity map between the arrays in the scop.
  ///
  /// @returns An identity map between the arrays in the scop.
  isl_union_map *getArrayIdentity() {
    isl_union_map *Maps = isl_union_map_empty(S->getParamSpace());

    for (auto &Item : S->arrays()) {
      ScopArrayInfo *Array = Item.second.get();
      isl_space *Space = Array->getSpace();
      Space = isl_space_map_from_set(Space);
      isl_map *Identity = isl_map_identity(Space);
      Maps = isl_union_map_add_map(Maps, Identity);
    }

    return Maps;
  }

  /// Create a default-initialized PPCG GPU program.
  ///
  /// @returns A new gpu grogram description.
  gpu_prog *createPPCGProg(ppcg_scop *PPCGScop) {

    if (!PPCGScop)
      return nullptr;

    auto PPCGProg = isl_calloc_type(S->getIslCtx(), struct gpu_prog);

    PPCGProg->ctx = S->getIslCtx();
    PPCGProg->scop = PPCGScop;
    PPCGProg->context = isl_set_copy(PPCGScop->context);
    PPCGProg->read = isl_union_map_copy(PPCGScop->reads);
    PPCGProg->may_write = isl_union_map_copy(PPCGScop->may_writes);
    PPCGProg->must_write = isl_union_map_copy(PPCGScop->must_writes);
    PPCGProg->tagged_must_kill =
        isl_union_map_copy(PPCGScop->tagged_must_kills);
    PPCGProg->to_inner = getArrayIdentity();
    PPCGProg->to_outer = getArrayIdentity();
    PPCGProg->may_persist = compute_may_persist(PPCGProg);
    PPCGProg->any_to_outer = nullptr;
    PPCGProg->array_order = nullptr;
    PPCGProg->n_stmts = std::distance(S->begin(), S->end());
    PPCGProg->stmts = getStatements();
    PPCGProg->n_array = std::distance(S->array_begin(), S->array_end());
    PPCGProg->array = isl_calloc_array(S->getIslCtx(), struct gpu_array_info,
                                       PPCGProg->n_array);

    createArrays(PPCGProg);

    return PPCGProg;
  }

  struct PrintGPUUserData {
    struct cuda_info *CudaInfo;
    struct gpu_prog *PPCGProg;
    std::vector<ppcg_kernel *> Kernels;
  };

  /// Print a user statement node in the host code.
  ///
  /// We use ppcg's printing facilities to print the actual statement and
  /// additionally build up a list of all kernels that are encountered in the
  /// host ast.
  ///
  /// @param P The printer to print to
  /// @param Options The printing options to use
  /// @param Node The node to print
  /// @param User A user pointer to carry additional data. This pointer is
  ///             expected to be of type PrintGPUUserData.
  ///
  /// @returns A printer to which the output has been printed.
  static __isl_give isl_printer *
  printHostUser(__isl_take isl_printer *P,
                __isl_take isl_ast_print_options *Options,
                __isl_take isl_ast_node *Node, void *User) {
    auto Data = (struct PrintGPUUserData *)User;
    auto Id = isl_ast_node_get_annotation(Node);

    if (Id) {
      bool IsUser = !strcmp(isl_id_get_name(Id), "user");

      // If this is a user statement, format it ourselves as ppcg would
      // otherwise try to call pet functionality that is not available in
      // Polly.
      if (IsUser) {
        P = isl_printer_start_line(P);
        P = isl_printer_print_ast_node(P, Node);
        P = isl_printer_end_line(P);
        isl_id_free(Id);
        isl_ast_print_options_free(Options);
        return P;
      }

      auto Kernel = (struct ppcg_kernel *)isl_id_get_user(Id);
      isl_id_free(Id);
      Data->Kernels.push_back(Kernel);
    }

    return print_host_user(P, Options, Node, User);
  }

  /// Print C code corresponding to the control flow in @p Kernel.
  ///
  /// @param Kernel The kernel to print
  void printKernel(ppcg_kernel *Kernel) {
    auto *P = isl_printer_to_str(S->getIslCtx());
    P = isl_printer_set_output_format(P, ISL_FORMAT_C);
    auto *Options = isl_ast_print_options_alloc(S->getIslCtx());
    P = isl_ast_node_print(Kernel->tree, P, Options);
    char *String = isl_printer_get_str(P);
    printf("%s\n", String);
    free(String);
    isl_printer_free(P);
  }

  /// Print C code corresponding to the GPU code described by @p Tree.
  ///
  /// @param Tree An AST describing GPU code
  /// @param PPCGProg The PPCG program from which @Tree has been constructed.
  void printGPUTree(isl_ast_node *Tree, gpu_prog *PPCGProg) {
    auto *P = isl_printer_to_str(S->getIslCtx());
    P = isl_printer_set_output_format(P, ISL_FORMAT_C);

    PrintGPUUserData Data;
    Data.PPCGProg = PPCGProg;

    auto *Options = isl_ast_print_options_alloc(S->getIslCtx());
    Options =
        isl_ast_print_options_set_print_user(Options, printHostUser, &Data);
    P = isl_ast_node_print(Tree, P, Options);
    char *String = isl_printer_get_str(P);
    printf("# host\n");
    printf("%s\n", String);
    free(String);
    isl_printer_free(P);

    for (auto Kernel : Data.Kernels) {
      printf("# kernel%d\n", Kernel->id);
      printKernel(Kernel);
    }
  }

  // Generate a GPU program using PPCG.
  //
  // GPU mapping consists of multiple steps:
  //
  //  1) Compute new schedule for the program.
  //  2) Map schedule to GPU (TODO)
  //  3) Generate code for new schedule (TODO)
  //
  // We do not use here the Polly ScheduleOptimizer, as the schedule optimizer
  // is mostly CPU specific. Instead, we use PPCG's GPU code generation
  // strategy directly from this pass.
  gpu_gen *generateGPU(ppcg_scop *PPCGScop, gpu_prog *PPCGProg) {

    auto PPCGGen = isl_calloc_type(S->getIslCtx(), struct gpu_gen);

    PPCGGen->ctx = S->getIslCtx();
    PPCGGen->options = PPCGScop->options;
    PPCGGen->print = nullptr;
    PPCGGen->print_user = nullptr;
    PPCGGen->build_ast_expr = &pollyBuildAstExprForStmt;
    PPCGGen->prog = PPCGProg;
    PPCGGen->tree = nullptr;
    PPCGGen->types.n = 0;
    PPCGGen->types.name = nullptr;
    PPCGGen->sizes = nullptr;
    PPCGGen->used_sizes = nullptr;
    PPCGGen->kernel_id = 0;

    // Set scheduling strategy to same strategy PPCG is using.
    isl_options_set_schedule_outer_coincidence(PPCGGen->ctx, true);
    isl_options_set_schedule_maximize_band_depth(PPCGGen->ctx, true);
    isl_options_set_schedule_whole_component(PPCGGen->ctx, false);

    isl_schedule *Schedule = get_schedule(PPCGGen);

    int has_permutable = has_any_permutable_node(Schedule);

    if (!has_permutable || has_permutable < 0) {
      Schedule = isl_schedule_free(Schedule);
    } else {
      Schedule = map_to_device(PPCGGen, Schedule);
      PPCGGen->tree = generate_code(PPCGGen, isl_schedule_copy(Schedule));
    }

    if (DumpSchedule) {
      isl_printer *P = isl_printer_to_str(S->getIslCtx());
      P = isl_printer_set_yaml_style(P, ISL_YAML_STYLE_BLOCK);
      P = isl_printer_print_str(P, "Schedule\n");
      P = isl_printer_print_str(P, "========\n");
      if (Schedule)
        P = isl_printer_print_schedule(P, Schedule);
      else
        P = isl_printer_print_str(P, "No schedule found\n");

      printf("%s\n", isl_printer_get_str(P));
      isl_printer_free(P);
    }

    if (DumpCode) {
      printf("Code\n");
      printf("====\n");
      if (PPCGGen->tree)
        printGPUTree(PPCGGen->tree, PPCGProg);
      else
        printf("No code generated\n");
    }

    isl_schedule_free(Schedule);

    return PPCGGen;
  }

  /// Free gpu_gen structure.
  ///
  /// @param PPCGGen The ppcg_gen object to free.
  void freePPCGGen(gpu_gen *PPCGGen) {
    isl_ast_node_free(PPCGGen->tree);
    isl_union_map_free(PPCGGen->sizes);
    isl_union_map_free(PPCGGen->used_sizes);
    free(PPCGGen);
  }

  /// Free the options in the ppcg scop structure.
  ///
  /// ppcg is not freeing these options for us. To avoid leaks we do this
  /// ourselves.
  ///
  /// @param PPCGScop The scop referencing the options to free.
  void freeOptions(ppcg_scop *PPCGScop) {
    free(PPCGScop->options->debug);
    PPCGScop->options->debug = nullptr;
    free(PPCGScop->options);
    PPCGScop->options = nullptr;
  }

  /// Generate code for a given GPU AST described by @p Root.
  ///
  /// @param An isl_ast_node pointing to the root of the GPU AST.
  void generateCode(__isl_take isl_ast_node *Root) {
    ScopAnnotator Annotator;
    Annotator.buildAliasScopes(*S);

    Region *R = &S->getRegion();

    simplifyRegion(R, DT, LI, RI);

    BasicBlock *EnteringBB = R->getEnteringBlock();

    PollyIRBuilder Builder = createPollyIRBuilder(EnteringBB, Annotator);

    GPUNodeBuilder NodeBuilder(Builder, Annotator, this, *DL, *LI, *SE, *DT,
                               *S);

    // Only build the run-time condition and parameters _after_ having
    // introduced the conditional branch. This is important as the conditional
    // branch will guard the original scop from new induction variables that
    // the SCEVExpander may introduce while code generating the parameters and
    // which may introduce scalar dependences that prevent us from correctly
    // code generating this scop.
    BasicBlock *StartBlock =
        executeScopConditionally(*S, this, Builder.getTrue());

    // TODO: Handle LICM
    // TODO: Verify run-time checks
    auto SplitBlock = StartBlock->getSinglePredecessor();
    Builder.SetInsertPoint(SplitBlock->getTerminator());
    NodeBuilder.addParameters(S->getContext());
    Builder.SetInsertPoint(&*StartBlock->begin());
    NodeBuilder.create(Root);
    NodeBuilder.finalizeSCoP(*S);
  }

  bool runOnScop(Scop &CurrentScop) override {
    S = &CurrentScop;
    LI = &getAnalysis<LoopInfoWrapperPass>().getLoopInfo();
    DT = &getAnalysis<DominatorTreeWrapperPass>().getDomTree();
    SE = &getAnalysis<ScalarEvolutionWrapperPass>().getSE();
    DL = &S->getRegion().getEntry()->getParent()->getParent()->getDataLayout();
    RI = &getAnalysis<RegionInfoPass>().getRegionInfo();

    auto PPCGScop = createPPCGScop();
    auto PPCGProg = createPPCGProg(PPCGScop);
    auto PPCGGen = generateGPU(PPCGScop, PPCGProg);

    if (PPCGGen->tree)
      generateCode(isl_ast_node_copy(PPCGGen->tree));

    freeOptions(PPCGScop);
    freePPCGGen(PPCGGen);
    gpu_prog_free(PPCGProg);
    ppcg_scop_free(PPCGScop);

    return true;
  }

  void printScop(raw_ostream &, Scop &) const override {}

  void getAnalysisUsage(AnalysisUsage &AU) const override {
    AU.addRequired<DominatorTreeWrapperPass>();
    AU.addRequired<RegionInfoPass>();
    AU.addRequired<ScalarEvolutionWrapperPass>();
    AU.addRequired<ScopDetection>();
    AU.addRequired<ScopInfoRegionPass>();
    AU.addRequired<LoopInfoWrapperPass>();

    AU.addPreserved<AAResultsWrapperPass>();
    AU.addPreserved<BasicAAWrapperPass>();
    AU.addPreserved<LoopInfoWrapperPass>();
    AU.addPreserved<DominatorTreeWrapperPass>();
    AU.addPreserved<GlobalsAAWrapperPass>();
    AU.addPreserved<PostDominatorTreeWrapperPass>();
    AU.addPreserved<ScopDetection>();
    AU.addPreserved<ScalarEvolutionWrapperPass>();
    AU.addPreserved<SCEVAAWrapperPass>();

    // FIXME: We do not yet add regions for the newly generated code to the
    //        region tree.
    AU.addPreserved<RegionInfoPass>();
    AU.addPreserved<ScopInfoRegionPass>();
  }
};
}

char PPCGCodeGeneration::ID = 1;

Pass *polly::createPPCGCodeGenerationPass() { return new PPCGCodeGeneration(); }

INITIALIZE_PASS_BEGIN(PPCGCodeGeneration, "polly-codegen-ppcg",
                      "Polly - Apply PPCG translation to SCOP", false, false)
INITIALIZE_PASS_DEPENDENCY(DependenceInfo);
INITIALIZE_PASS_DEPENDENCY(DominatorTreeWrapperPass);
INITIALIZE_PASS_DEPENDENCY(LoopInfoWrapperPass);
INITIALIZE_PASS_DEPENDENCY(RegionInfoPass);
INITIALIZE_PASS_DEPENDENCY(ScalarEvolutionWrapperPass);
INITIALIZE_PASS_DEPENDENCY(ScopDetection);
INITIALIZE_PASS_END(PPCGCodeGeneration, "polly-codegen-ppcg",
                    "Polly - Apply PPCG translation to SCOP", false, false)
