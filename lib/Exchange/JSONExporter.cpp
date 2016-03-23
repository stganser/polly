//===-- JSONExporter.cpp  - Export Scops as JSON  -------------------------===//
//
//                     The LLVM Compiler Infrastructure
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//
//
// Export the Scops build by ScopInfo pass as a JSON file.
//
//===----------------------------------------------------------------------===//

#include "polly/DependenceInfo.h"
#include "polly/LinkAllPasses.h"
#include "polly/Options.h"
#include "polly/ScopInfo.h"
#include "polly/ScopPass.h"
#include "polly/Support/ScopLocation.h"
#include "polly/ScheduleOptimizer.h"
#include "llvm/ADT/Statistic.h"
#include "llvm/ADT/DenseSet.h"
#include "llvm/Analysis/RegionInfo.h"
#include "llvm/IR/Module.h"
#include "llvm/Support/FileSystem.h"
#include "llvm/Support/MemoryBuffer.h"
#include "llvm/Support/ToolOutputFile.h"
#include "llvm/Support/ErrorHandling.h"
#include "isl/constraint.h"
#include "isl/map.h"
#include "isl/union_map.h"
#include "isl/printer.h"
#include "isl/set.h"
#include "isl/union_set.h"
#include "isl/schedule.h"
#include "isl/schedule_node.h"
#include "isl/schedule_type.h"
#include "json/reader.h"
#include "json/writer.h"
#include <memory>
#include <string>
#include <system_error>

using namespace llvm;
using namespace polly;

#define DEBUG_TYPE "polly-import-jscop"

STATISTIC(NewAccessMapFound, "Number of updated access functions");

namespace {
static cl::opt<std::string>
    ImportDir("polly-import-jscop-dir",
              cl::desc("The directory to import the .jscop files from."),
              cl::Hidden, cl::value_desc("Directory path"), cl::ValueRequired,
              cl::init("."), cl::cat(PollyCategory));

static cl::opt<std::string>
    ImportPostfix("polly-import-jscop-postfix",
                  cl::desc("Postfix to append to the import .jsop files."),
                  cl::Hidden, cl::value_desc("File postfix"), cl::ValueRequired,
                  cl::init(""), cl::cat(PollyCategory));

static cl::opt<bool> LoadScheduleTree("polly-import-jscop-read-schedule-tree",
                                      cl::desc("Load the new schedule as a schedule tree."),
                                      cl::init(false), cl::ZeroOrMore,
                                      cl::cat(PollyCategory));

struct JSONExporter : public ScopPass {
  static char ID;
  explicit JSONExporter() : ScopPass(ID) {}

  std::string getFileName(Scop &S) const;
  Json::Value getJSON(Scop &S) const;

  /// @brief Export the SCoP @p S to a JSON file.
  bool runOnScop(Scop &S) override;

  /// @brief Print the SCoP @p S as it is exported.
  void printScop(raw_ostream &OS, Scop &S) const override;

  /// @brief Register all analyses and transformation required.
  void getAnalysisUsage(AnalysisUsage &AU) const override;
};

struct JSONImporter : public ScopPass {
  static char ID;
  std::vector<std::string> newAccessStrings;
  explicit JSONImporter() : ScopPass(ID) {}

  std::string getFileName(Scop &S) const;

  /// @brief Import new access functions for SCoP @p S from a JSON file.
  bool runOnScop(Scop &S) override;

  /// @brief Print the SCoP @p S and the imported access functions.
  void printScop(raw_ostream &OS, Scop &S) const override;

  /// @brief Register all analyses and transformation required.
  void getAnalysisUsage(AnalysisUsage &AU) const override;
};
}

char JSONExporter::ID = 0;
std::string JSONExporter::getFileName(Scop &S) const {
  std::string FunctionName = S.getRegion().getEntry()->getParent()->getName();
  std::string FileName = FunctionName + "___" + S.getNameStr() + ".jscop";
  return FileName;
}

void JSONExporter::printScop(raw_ostream &OS, Scop &S) const { S.print(OS); }

Json::Value JSONExporter::getJSON(Scop &S) const {
  Json::Value root;
  unsigned LineBegin, LineEnd;
  std::string FileName;

  getDebugLocation(&S.getRegion(), LineBegin, LineEnd, FileName);
  std::string Location;
  if (LineBegin != (unsigned)-1)
    Location = FileName + ":" + std::to_string(LineBegin) + "-" +
               std::to_string(LineEnd);

  root["name"] = S.getRegion().getNameStr();
  root["context"] = S.getContextStr();
  if (LineBegin != (unsigned)-1)
    root["location"] = Location;
  root["statements"];

  for (ScopStmt &Stmt : S) {
    Json::Value statement;

    statement["name"] = Stmt.getBaseName();
    statement["domain"] = Stmt.getDomainStr();
    statement["schedule"] = Stmt.getScheduleStr();
    statement["accesses"];

    for (MemoryAccess *MA : Stmt) {
      Json::Value access;

      access["kind"] = MA->isRead() ? "read" : "write";
      access["relation"] = MA->getOriginalAccessRelationStr();

      statement["accesses"].append(access);
    }

    root["statements"].append(statement);
  }

  return root;
}

bool JSONExporter::runOnScop(Scop &S) {
  Region &R = S.getRegion();

  std::string FileName = ImportDir + "/" + getFileName(S);

  Json::Value jscop = getJSON(S);
  Json::StyledWriter writer;
  std::string fileContent = writer.write(jscop);

  // Write to file.
  std::error_code EC;
  tool_output_file F(FileName, EC, llvm::sys::fs::F_Text);

  std::string FunctionName = R.getEntry()->getParent()->getName();
  errs() << "Writing JScop '" << R.getNameStr() << "' in function '"
         << FunctionName << "' to '" << FileName << "'.\n";

  if (!EC) {
    F.os() << fileContent;
    F.os().close();
    if (!F.os().has_error()) {
      errs() << "\n";
      F.keep();
      return false;
    }
  }

  errs() << "  error opening file for writing!\n";
  F.os().clear_error();

  return false;
}

void JSONExporter::getAnalysisUsage(AnalysisUsage &AU) const {
  AU.setPreservesAll();
  AU.addRequired<ScopInfo>();
}

Pass *polly::createJSONExporterPass() { return new JSONExporter(); }

char JSONImporter::ID = 0;
std::string JSONImporter::getFileName(Scop &S) const {
  std::string FunctionName = S.getRegion().getEntry()->getParent()->getName();
  std::string FileName = FunctionName + "___" + S.getNameStr() + ".jscop";

  if (ImportPostfix != "")
    FileName += "." + ImportPostfix;

  return FileName;
}

void JSONImporter::printScop(raw_ostream &OS, Scop &S) const {
  S.print(OS);
  for (std::vector<std::string>::const_iterator I = newAccessStrings.begin(),
                                                E = newAccessStrings.end();
       I != E; I++)
    OS << "New access function '" << *I << "'detected in JSCOP file\n";
}

typedef Dependences::StatementToIslMapTy StatementToIslMapTy;

__isl_give isl_set *getDomainForStmtId(const char *stmtId, Scop &S) {
  for (ScopStmt &stmt : S) {
    if (!std::strcmp(stmt.getBaseName(), stmtId)) {
      return stmt.getDomain();
    }
  }
  return nullptr;
}

isl_stat foreachSetHelper(isl_set *s, void *user) {
  std::function<void(isl_set*)> *f = (std::function<void(isl_set*)> *) user;
  (*f)(s);
  return isl_stat_ok;
}

void callLambda(__isl_keep isl_union_set *s, std::function<void(isl_set*)> &f) {
  isl_union_set_foreach_set(s, foreachSetHelper, &f);
}

isl_stat foreachMapHelper(isl_map *m, void *user) {
  std::function<void(isl_map*)> *f = (std::function<void(isl_map*)> *) user;
  (*f)(m);
  return isl_stat_ok;
}

void callLambda(__isl_keep isl_union_map *m, std::function<void(isl_map*)> &f) {
  isl_union_map_foreach_map(m, foreachMapHelper, &f);
}

__isl_give isl_schedule_node *scheduleMapHelper(isl_schedule_node *n, void *user) {
  std::function<isl_schedule_node*(isl_schedule_node*)> *f
      = (std::function<isl_schedule_node*(isl_schedule_node*)> *) user;
  return (*f)(n);
}

__isl_give isl_schedule *callLambda(__isl_take isl_schedule *s,
                std :: function<isl_schedule_node*(isl_schedule_node*)> &f) {
  return isl_schedule_map_schedule_node_bottom_up(s, scheduleMapHelper, &f);
}

DenseSet<const char *> getTupleNamesFromUnionSet(__isl_keep isl_union_set *s) {
  DenseSet<const char *> names;
  std::function<void(isl_set *)> lambda = [&names](isl_set *s) {
    names.insert(isl_set_get_tuple_name(s));
    isl_set_free(s);
  };
  callLambda(s, lambda);
  return names;
}

__isl_give isl_union_set *buildDomainForSttmnts(
    DenseSet<const char*> &stmtNames, Scop &S) {
  isl_union_set *domain = S.getDomains();
  isl_union_set *result
      = isl_union_set_empty(isl_union_set_get_space(domain));
  isl_union_set_free(domain);
  for (const char* stmtName : stmtNames) {
    for (ScopStmt &S : S) {
      if (!std::strcmp(stmtName, S.getBaseName())) {
        result = isl_union_set_add_set(result, S.getDomain());
      }
    }
  }
  return result;
}

ScopStmt& getStmtForName(const char *name, Scop &S) {
  for (ScopStmt &Stmt : S) {
    if (!std::strcmp(name, Stmt.getBaseName())) {
      return Stmt;
    }
  }
  report_fatal_error("Wat?", true);
}

__isl_give isl_union_map *replaceIdsInUnionMap(__isl_take isl_union_map *m,
                                               Scop &S) {
  isl_union_map *result = isl_union_map_empty(S.getParamSpace());
  std::function<void(isl_map *)> lambda = [&result, &S](isl_map *m) {
    const char* tupleName = isl_map_get_tuple_name(m, isl_dim_in);
    ScopStmt &sttmt = getStmtForName(tupleName, S);
    isl_space *stmtSpace = sttmt.getDomainSpace();
    isl_map *newMap = isl_map_set_tuple_id(m, isl_dim_in,
                              isl_space_get_tuple_id(stmtSpace, isl_dim_set));
    for (unsigned i = 0; i < isl_space_dim(stmtSpace, isl_dim_param); ++i) {
      isl_id *id = isl_space_get_dim_id(stmtSpace, isl_dim_param, i);
      newMap = isl_map_set_dim_id(newMap, isl_dim_param, i, id);
    }
    isl_space_free(stmtSpace);
    result = isl_union_map_add_map(result, newMap);
  };
  callLambda(m, lambda);
  isl_union_map_free(m);
  return result;
}

__isl_give isl_schedule *setCoincidence(__isl_take isl_schedule *sched,
                                        __isl_keep isl_schedule_node *oldNode) {
  isl_schedule_node *root = isl_schedule_get_root(sched);
  std::function<isl_schedule_node*(isl_schedule_node*)> lambda
      = [&root, &oldNode](isl_schedule_node *n) {
    isl_schedule_node *res = n;
    if (isl_schedule_node_has_parent(n)) {
      isl_schedule_node *parent = isl_schedule_node_parent(
            isl_schedule_node_copy(n));
      if (isl_schedule_node_is_equal(root, parent)) {
        for (unsigned i = 0; i < isl_schedule_node_band_n_member(res); ++i) {
          int coincidence = isl_schedule_node_band_member_get_coincident(
                oldNode, i);
          res = isl_schedule_node_band_member_set_coincident(res, i,
                                                             coincidence);
        }
      }
      isl_schedule_node_free(parent);
    }
    return res;
  };
  isl_schedule *result = callLambda(sched, lambda);
  isl_schedule_node_free(root);
  return result;
}

__isl_give isl_schedule *setPermutability(__isl_take isl_schedule *sched,
                                      __isl_keep isl_schedule_node *oldNode) {
  isl_schedule_node *root = isl_schedule_get_root(sched);
  std::function<isl_schedule_node*(isl_schedule_node*)> lambda
      = [&root, &oldNode](isl_schedule_node *n) {
    isl_schedule_node *res = n;
    if (isl_schedule_node_get_type(n) == isl_schedule_node_band
        && isl_schedule_node_has_parent(n)) {
      isl_schedule_node *parent = isl_schedule_node_parent(
            isl_schedule_node_copy(n));
      if (isl_schedule_node_is_equal(root, parent)) {
        res = isl_schedule_node_band_set_permutable(n,
                                isl_schedule_node_band_get_permutable(oldNode));
      }
      isl_schedule_node_free(parent);
    }
    return res;
  };
  isl_schedule_node_free(root);
  isl_schedule *result = callLambda(sched, lambda);
  return result;
}

__isl_give isl_schedule *rebuildSchedule(__isl_take isl_schedule_node *n,
                                  Scop &S, __isl_keep isl_union_set *domain) {

  isl_schedule *result = nullptr;

  switch (isl_schedule_node_get_type(n)) {
    case isl_schedule_node_type::isl_schedule_node_filter : {
      isl_union_set *filterSet = isl_schedule_node_filter_get_filter(n);
      DenseSet<const char*> stmtNames = getTupleNamesFromUnionSet(filterSet);
      isl_union_set *newDomain = buildDomainForSttmnts(stmtNames, S);
      isl_union_set_free(filterSet);

      result = rebuildSchedule(isl_schedule_node_get_child(n, 0), S,
                                 newDomain);
      isl_union_set_free(newDomain);
      break;
    }

    case isl_schedule_node_type::isl_schedule_node_domain : {
      result = rebuildSchedule(isl_schedule_node_get_child(n, 0), S, domain);
      break;
    }

    case isl_schedule_node_type::isl_schedule_node_band : {
      isl_union_map *sched =
          isl_schedule_node_band_get_partial_schedule_union_map(n);
      sched = replaceIdsInUnionMap(sched, S);
      result = rebuildSchedule(isl_schedule_node_get_child(n, 0), S, domain);
      result = isl_schedule_insert_partial_schedule(result,
                                isl_multi_union_pw_aff_from_union_map(sched));
      result = setCoincidence(result, n);
      result = setPermutability(result, n);
      break;
    }

    case isl_schedule_node_type::isl_schedule_node_sequence : {
      isl_schedule_node *fstChild = isl_schedule_node_get_child(n, 0);
      result = rebuildSchedule(fstChild, S, domain);
      for (int i = 1; i < isl_schedule_node_n_children(n); ++i) {
        isl_schedule *childSched = rebuildSchedule(
              isl_schedule_node_get_child(n, i), S, domain);
        result = isl_schedule_sequence(result, childSched);
      }
      break;
    }

    case isl_schedule_node_type::isl_schedule_node_set : {
      isl_schedule_node *fstChild = isl_schedule_node_get_child(n, 0);
      result = rebuildSchedule(fstChild, S, domain);
      for (int i = 1; i < isl_schedule_node_n_children(n); ++i) {
        isl_schedule *childSched = rebuildSchedule(
              isl_schedule_node_get_child(n, i), S, domain);
        result = isl_schedule_set(result, childSched);
      }
      break;
    }

    case isl_schedule_node_type::isl_schedule_node_leaf : {
      result = isl_schedule_from_domain(isl_union_set_copy(domain));
      break;
    }

    default : {
      report_fatal_error(
          "schedule rebuilding is unimplemented for the current node type.",
            false);
    }
  }
  isl_schedule_node_free(n);
  return result;
}

__isl_give isl_schedule *rebuildSchedule(__isl_take isl_schedule *sched,
                                         Scop &S) {
  isl_union_set *domain = S.getDomains();
  isl_schedule *result = rebuildSchedule(isl_schedule_get_root(sched), S,
                                         domain);
  isl_printer *pr = isl_printer_to_str(S.getIslCtx());
  pr = isl_printer_print_schedule(pr, sched);
  errs() << "Original schedule: " << isl_printer_get_str(pr) << '\n';
  pr = isl_printer_flush(pr);
  pr = isl_printer_print_schedule(pr, result);
  errs() << "Rebuilt schedule: " << isl_printer_get_str(pr) << '\n';
  isl_printer_free(pr);
  isl_schedule_free(sched);
  isl_union_set_free(domain);
  return result;
}

bool JSONImporter::runOnScop(Scop &S) {
  Region &R = S.getRegion();
  const Dependences &D = getAnalysis<DependenceInfo>().getDependences();
  const DataLayout &DL =
      S.getRegion().getEntry()->getParent()->getParent()->getDataLayout();

  std::string FileName = ImportDir + "/" + getFileName(S);

  std::string FunctionName = R.getEntry()->getParent()->getName();
  errs() << "Reading JScop '" << R.getNameStr() << "' in function '"
         << FunctionName << "' from '" << FileName << "'.\n";
  ErrorOr<std::unique_ptr<MemoryBuffer>> result =
      MemoryBuffer::getFile(FileName);
  std::error_code ec = result.getError();

  if (ec) {
    std::string errMsg = "File could not be read: " + ec.message();
    errs() << errMsg << '\n';
    report_fatal_error(errMsg, false);
  }

  Json::Reader reader;
  Json::Value jscop;

  bool parsingSuccessful = reader.parse(result.get()->getBufferStart(), jscop);

  if (!parsingSuccessful) {
    std::string errMsg = "JSCoP file could not be parsed";
    errs() << errMsg << '\n';
    report_fatal_error(errMsg, false);
  }

  isl_set *OldContext = S.getContext();
  isl_set *NewContext =
      isl_set_read_from_str(S.getIslCtx(), jscop["context"].asCString());

  for (unsigned i = 0; i < isl_set_dim(OldContext, isl_dim_param); i++) {
    isl_id *id = isl_set_get_dim_id(OldContext, isl_dim_param, i);
    NewContext = isl_set_set_dim_id(NewContext, isl_dim_param, i, id);
  }

  isl_set_free(OldContext);
  S.setContext(NewContext);

  StatementToIslMapTy NewSchedule;

  isl_schedule *scheduleTree;

  if (LoadScheduleTree) {
    errs() << "Loading schedule tree.\n";
    Json :: Value schedTreeVal = jscop["schedTree"];
    scheduleTree = isl_schedule_read_from_str(S.getIslCtx(), schedTreeVal.asCString());
    isl_union_map *schedMap = isl_schedule_get_map(scheduleTree);
    schedMap = replaceIdsInUnionMap(schedMap, S);
    std::function<void(isl_map*)> lambda = [&NewSchedule, &S](isl_map *m) {
      for (ScopStmt &Stmt : S) {
        if (!strcmp(Stmt.getBaseName(), isl_map_get_tuple_name(m, isl_dim_in))) {
          NewSchedule[&Stmt] = m;
        }
      }
    };

    callLambda(schedMap, lambda);
    isl_union_map_free(schedMap);
  } else {
    int index = 0;
    for (ScopStmt &Stmt : S) {
      Json::Value schedule = jscop["statements"][index]["schedule"];
      isl_map *m = isl_map_read_from_str(S.getIslCtx(), schedule.asCString());
      isl_space *Space = Stmt.getDomainSpace();

      // Copy the old tuple id. This is necessary to retain the user pointer,
      // that stores the reference to the ScopStmt this schedule belongs to.
      m = isl_map_set_tuple_id(m, isl_dim_in,
                               isl_space_get_tuple_id(Space, isl_dim_set));
      for (unsigned i = 0; i < isl_space_dim(Space, isl_dim_param); i++) {
        isl_id *id = isl_space_get_dim_id(Space, isl_dim_param, i);
        m = isl_map_set_dim_id(m, isl_dim_param, i, id);
      }
      isl_space_free(Space);
      NewSchedule[&Stmt] = m;
      index++;
    }
  }

  if (!D.isValidSchedule(S, &NewSchedule)) {
    std::string errMsg = "JScop file contains a schedule that changes the "
                         "dependences. Use -disable-polly-legality to continue "
                         "anyways";
    errs() << errMsg << '\n';
    for (StatementToIslMapTy::iterator SI = NewSchedule.begin(),
                                       SE = NewSchedule.end();
         SI != SE; ++SI)
      isl_map_free(SI->second);
    report_fatal_error(errMsg, false);
  }

  if (LoadScheduleTree) {
    scheduleTree = rebuildSchedule(scheduleTree, S);
    scheduleTree = ScheduleTreeOptimizer::optimizeSchedule(scheduleTree);
    S.setScheduleTree(scheduleTree);
    S.markAsOptimized();
    for (ScopStmt &Stmt : S) {
      isl_map_free(NewSchedule[&Stmt]);
    }
  } else {
    auto ScheduleMap = isl_union_map_empty(S.getParamSpace());
    for (ScopStmt &Stmt : S) {
      if (NewSchedule.find(&Stmt) != NewSchedule.end())
        ScheduleMap = isl_union_map_add_map(ScheduleMap, NewSchedule[&Stmt]);
      else
        ScheduleMap = isl_union_map_add_map(ScheduleMap, Stmt.getSchedule());
    }
    S.setSchedule(ScheduleMap);
  }

  int statementIdx = 0;
  for (ScopStmt &Stmt : S) {
    int memoryAccessIdx = 0;
    for (MemoryAccess *MA : Stmt) {
      Json::Value accesses = jscop["statements"][statementIdx]["accesses"]
                                  [memoryAccessIdx]["relation"];
      isl_map *newAccessMap =
          isl_map_read_from_str(S.getIslCtx(), accesses.asCString());
      isl_map *currentAccessMap = MA->getAccessRelation();

      if (isl_map_dim(newAccessMap, isl_dim_param) !=
          isl_map_dim(currentAccessMap, isl_dim_param)) {
        std::string errMsg = "JScop file changes the number of parameter "
                             "dimensions";
        errs() << errMsg << '\n';
        isl_map_free(currentAccessMap);
        isl_map_free(newAccessMap);
        report_fatal_error(errMsg, false);
      }

      isl_id *OutId = isl_map_get_tuple_id(currentAccessMap, isl_dim_out);
      newAccessMap = isl_map_set_tuple_id(newAccessMap, isl_dim_out, OutId);

      if (MA->isArrayKind()) {
        // We keep the old alignment, thus we cannot allow accesses to memory
        // locations that were not accessed before if the alignment of the
        // access is not the default alignment.
        bool SpecialAlignment = true;
        if (LoadInst *LoadI = dyn_cast<LoadInst>(MA->getAccessInstruction())) {
          SpecialAlignment =
              DL.getABITypeAlignment(LoadI->getType()) != LoadI->getAlignment();
        } else if (StoreInst *StoreI =
                       dyn_cast<StoreInst>(MA->getAccessInstruction())) {
          SpecialAlignment =
              DL.getABITypeAlignment(StoreI->getValueOperand()->getType()) !=
              StoreI->getAlignment();
        }

        if (SpecialAlignment) {
          isl_set *newAccessSet = isl_map_range(isl_map_copy(newAccessMap));
          isl_set *currentAccessSet =
              isl_map_range(isl_map_copy(currentAccessMap));
          bool isSubset = isl_set_is_subset(newAccessSet, currentAccessSet);
          isl_set_free(newAccessSet);
          isl_set_free(currentAccessSet);

          if (!isSubset) {
            std::string errMsg = "JScop file changes the accessed memory";
            errs() << errMsg << '\n';
            isl_map_free(currentAccessMap);
            isl_map_free(newAccessMap);
            report_fatal_error(errMsg, false);
          }
        }
      }

      // We need to copy the isl_ids for the parameter dimensions to the new
      // map. Without doing this the current map would have different
      // ids then the new one, even though both are named identically.
      for (unsigned i = 0; i < isl_map_dim(currentAccessMap, isl_dim_param);
           i++) {
        isl_id *id = isl_map_get_dim_id(currentAccessMap, isl_dim_param, i);
        newAccessMap = isl_map_set_dim_id(newAccessMap, isl_dim_param, i, id);
      }

      // Copy the old tuple id. This is necessary to retain the user pointer,
      // that stores the reference to the ScopStmt this access belongs to.
      isl_id *Id = isl_map_get_tuple_id(currentAccessMap, isl_dim_in);
      newAccessMap = isl_map_set_tuple_id(newAccessMap, isl_dim_in, Id);

      if (!isl_map_has_equal_space(currentAccessMap, newAccessMap)) {
        std::string errMsg = "JScop file contains access function with incompatible "
                "dimensions";
        errs() << errMsg << '\n';
        isl_map_free(currentAccessMap);
        isl_map_free(newAccessMap);
        report_fatal_error(errMsg, false);
      }

      auto NewAccessDomain = isl_map_domain(isl_map_copy(newAccessMap));
      auto CurrentAccessDomain = isl_map_domain(isl_map_copy(currentAccessMap));

      NewAccessDomain =
          isl_set_intersect_params(NewAccessDomain, S.getContext());
      CurrentAccessDomain =
          isl_set_intersect_params(CurrentAccessDomain, S.getContext());

      if (isl_set_is_subset(CurrentAccessDomain, NewAccessDomain) ==
          isl_bool_false) {
        std::string errMsg = "Mapping not defined for all iteration domain "
                             "elements";
        errs() << errMsg << '\n';
        isl_set_free(CurrentAccessDomain);
        isl_set_free(NewAccessDomain);
        isl_map_free(currentAccessMap);
        isl_map_free(newAccessMap);
        report_fatal_error(errMsg, false);
      }

      isl_set_free(CurrentAccessDomain);
      isl_set_free(NewAccessDomain);

      if (!isl_map_is_equal(newAccessMap, currentAccessMap)) {
        // Statistics.
        ++NewAccessMapFound;
        newAccessStrings.push_back(accesses.asCString());
        MA->setNewAccessRelation(newAccessMap);
      } else {
        isl_map_free(newAccessMap);
      }
      isl_map_free(currentAccessMap);
      memoryAccessIdx++;
    }
    statementIdx++;
  }

  return false;
}

void JSONImporter::getAnalysisUsage(AnalysisUsage &AU) const {
  ScopPass::getAnalysisUsage(AU);
  AU.addRequired<DependenceInfo>();
}

Pass *polly::createJSONImporterPass() { return new JSONImporter(); }

INITIALIZE_PASS_BEGIN(JSONExporter, "polly-export-jscop",
                      "Polly - Export Scops as JSON"
                      " (Writes a .jscop file for each Scop)",
                      false, false);
INITIALIZE_PASS_DEPENDENCY(DependenceInfo)
INITIALIZE_PASS_END(JSONExporter, "polly-export-jscop",
                    "Polly - Export Scops as JSON"
                    " (Writes a .jscop file for each Scop)",
                    false, false)

INITIALIZE_PASS_BEGIN(JSONImporter, "polly-import-jscop",
                      "Polly - Import Scops from JSON"
                      " (Reads a .jscop file for each Scop)",
                      false, false);
INITIALIZE_PASS_DEPENDENCY(DependenceInfo)
INITIALIZE_PASS_END(JSONImporter, "polly-import-jscop",
                    "Polly - Import Scops from JSON"
                    " (Reads a .jscop file for each Scop)",
                    false, false)
