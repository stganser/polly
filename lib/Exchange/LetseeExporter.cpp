//===-- LetseeExporter.cpp  - Export Scops as Candl files -----------===//
//
//                     The LLVM Compiler Infrastructure
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//
//
// Export the Scops build by ScopInfo pass as a Candl file.
//
//===----------------------------------------------------------------------===//

#include "polly/LinkAllPasses.h"
#include "polly/Dependences.h"
#include "polly/Options.h"
#include "polly/ScopInfo.h"
#include "polly/ScopPass.h"

#include "llvm/ADT/Statistic.h"
#include "llvm/Analysis/RegionInfo.h"
#include "llvm/Support/FileSystem.h"
#include "llvm/Support/MemoryBuffer.h"
#include "llvm/Support/ToolOutputFile.h"

#include "isl/set.h"
#include "isl/map.h"
#include "isl/constraint.h"
#include "isl/printer.h"

#include <sstream>
#include <iostream>

using namespace llvm;
using namespace polly;

#define DEBUG_TYPE "polly-export-letsee"

struct LetseeExporter : public ScopPass {
  static char ID;
  explicit LetseeExporter() : ScopPass(ID) {}

  std::string getFileName(Scop &S) const;
  virtual bool runOnScop(Scop &S);
  void printScop(raw_ostream &OS, Scop &S) const;
  void getAnalysisUsage(AnalysisUsage &AU) const;

private:
  std::string getCloogInput(Scop &S) const;
  std::string printSet(isl_printer **p, isl_set *s) const;
  std::string printMap(isl_printer **p, isl_map *m) const;
  std::string getArrayId(MemoryAccess *MA) const;
  int outputAccessRelation(MemoryAccess *MA,
                           std::map<std::string, int> &arrayId2String,
                           isl_printer **p, std::stringstream &result) const;
  void outputMatrix(std::string &matrix, std::stringstream &result) const;
};

char LetseeExporter::ID = 0;

std::string LetseeExporter::getFileName(Scop &S) const {
  std::string FunctionName = S.getRegion().getEntry()->getParent()->getName();
  std::string FileName = FunctionName + "___" + S.getNameStr() + ".candl";
  return FileName;
}

void LetseeExporter::printScop(raw_ostream &OS, Scop &S) const { S.print(OS); }

std::string LetseeExporter::printSet(isl_printer **p, isl_set *s) const {
  *p = isl_printer_print_set(*p, s);
  char *char_str = isl_printer_get_str(*p);
  *p = isl_printer_flush(*p);
  std::string str = std::string(char_str);
  free(char_str);
  return str;
}

std::string LetseeExporter::printMap(isl_printer **p, isl_map *m) const {
  *p = isl_printer_print_map(*p, m);
  char *char_str = isl_printer_get_str(*p);
  *p = isl_printer_flush(*p);
  std::string str = std::string(char_str);
  free(char_str);
  return str;
}

std::string LetseeExporter::getArrayId(MemoryAccess *MA) const {
  isl_id *arrayId = MA->getArrayId();
  std::string arrayIdString = std::string(isl_id_get_name(arrayId));
  isl_id_free(arrayId);
  return arrayIdString;
}

int LetseeExporter::outputAccessRelation(MemoryAccess *MA,
                              std::map<std::string, int> &arrayId2String,
                              isl_printer **p, std::stringstream &result) const {
  int arrayIdNum = arrayId2String.find(getArrayId(MA))->second;
  isl_map *accessRelation = MA->getAccessRelation();
  std::string accessRelStr = printMap(p, accessRelation);
  isl_map_free(accessRelation);

  int numLines = 0;

  std::string accessRelVectStr;
  std::stringstream ss(accessRelStr);

  bool isFirstLine = true;

  for (int i = 0; i < 3; ++i)
    std::getline(ss, accessRelVectStr, '\n');

  std::stringstream ssArrayIdNum("");
  ssArrayIdNum << arrayIdNum;

  while (std::getline(ss, accessRelVectStr, '\n')) {
    int firstSpaceIndex = accessRelVectStr.find_first_of(' ');
    std::cout << "old: " << accessRelVectStr << std::endl;
    //accessRelVectStr
    //    = accessRelVectStr.replace(0, firstSpaceIndex,
    //                               (isFirstLine ? ssArrayIdNum.str() : "0"));
    accessRelVectStr = accessRelVectStr.replace(0, firstSpaceIndex, "");
    int firstNotOfSpaceIndex = accessRelVectStr.find_first_not_of(" ");
    accessRelVectStr = accessRelVectStr.replace(0, firstNotOfSpaceIndex, "");
    accessRelVectStr
        = accessRelVectStr.replace(0, firstSpaceIndex,
                                   (isFirstLine ? ssArrayIdNum.str() : "0"));
    std::cout << "new: " << accessRelVectStr << std::endl;
    result << accessRelVectStr;

    if (isFirstLine)
      result << " # " << getArrayId(MA);
    result << std::endl;

    ++numLines;
    isFirstLine = false;
  }
  return numLines;
}

void LetseeExporter::outputMatrix(
    std::string &matrix, std::stringstream &result) const {
  std::stringstream ss(matrix);
  std::string currentLine;

  for (int i = 0; i < 2; ++i)
    std::getline(ss, currentLine, '\n');

  while (std::getline(ss, currentLine, '\n'))
    result << currentLine << std::endl;
}

std::string LetseeExporter::getCloogInput(Scop &S) const {
    std::stringstream result;
    std::string blockDelimiter = "####################";
  isl_ctx *ctx = S.getIslCtx();
  isl_printer *p = isl_printer_to_str(ctx);
  p = isl_printer_set_output_format(p, ISL_FORMAT_POLYLIB);

  //isl_set *context = S.getAssumedContext();
  isl_set *context = S.getContext();
  result << blockDelimiter << std::endl << "# Context" << std::endl;
  {
    std::string contextMatrixStr = printSet(&p, context);
    outputMatrix(contextMatrixStr, result);
  }
  isl_set_free(context);

  result << std::endl << std::endl << blockDelimiter << std::endl
         << "# Number of statements" << std::endl;
  result << S.getSize() << std::endl << std::endl;

  int statementCount = 0;

  std::map<std::string, int> arrayId2String;
  int arrayCount = 1;

  for (Scop::iterator SI = S.begin(), SE = S.end(); SI != SE; ++SI) {
    ScopStmt *Stmt = *SI;

    for (MemoryAccess *MA : *Stmt) {
      std::string arrayIdString = getArrayId(MA);

      if (arrayId2String.find(arrayIdString) == arrayId2String.end())
        arrayId2String[arrayIdString] = arrayCount++;
    }
  }

  for (Scop::iterator SI = S.begin(), SE = S.end(); SI != SE; ++SI) {
    statementCount++;
    ScopStmt *Stmt = *SI;

    result << blockDelimiter << std::endl << "# Statement " << statementCount
           << std::endl;
    result << "# Statement Type" << std::endl;

    char reductionType = 'A';

    for (MemoryAccess *MA : *Stmt) {
      if (MA->isWrite()) {
        switch (MA->getReductionType()) {
          case polly::MemoryAccess::RT_MUL:
            reductionType = 'T';
            break;

          case polly::MemoryAccess::RT_ADD:
            reductionType = 'P';

          default:
            reductionType = 'A';
        }
      }
    }
    result << reductionType << std::endl << std::endl;

    isl_set *stmtDomain = Stmt->getDomain();
    result << "# Iteration Domain" << std::endl;
    std::string iterationDomainStr = printSet(&p, stmtDomain);
    outputMatrix(iterationDomainStr, result);

    result << std::endl << std::endl << "# Loop Labels" << std::endl;

    {
      uint i;

      for (i = 1; i <= Stmt->getNumIterators() - 1; ++i) {
        result << i << ' ';
      }

      if (Stmt->getNumIterators() > 0)
        result << i << std::endl;
    }
    isl_set_free(stmtDomain);
    result << std::endl;

    std::string firstMARepr;
    bool atLeastOneAccess = false;

    for (MemoryAccess *MA : *Stmt) {
      isl_map *accessRelation = MA->getAccessRelation();
      firstMARepr = printMap(&p, accessRelation);
      isl_map_free(accessRelation);
      atLeastOneAccess = true;
      break;
    }

    if (!atLeastOneAccess) {
      errs() << "Failed to export to Letsee format!\n";
      return NULL;
    }
    std::stringstream ss(firstMARepr);
    std::string accessRelDimStr;

    for (int i = 0; i < 3; ++i)
      std::getline(ss, accessRelDimStr, '\n');
    std::string accessRelNumRowsStr, accessRelNumColsStr;
    std::istringstream(accessRelDimStr) >> accessRelNumRowsStr >> accessRelNumColsStr;
    int accessRelNumCols;
    std::istringstream(accessRelNumColsStr) >> accessRelNumCols;

    result << "# Written items" << std::endl;
    std::stringstream accessRelSS;
    int accessRelNumRows  = 0;

    for (MemoryAccess *MA : *Stmt) {

      if (MA->isWrite()) {
        accessRelNumRows
            += outputAccessRelation(MA, arrayId2String, &p, accessRelSS);
      }
    }
    result << accessRelNumRows  << " " << (accessRelNumCols - 1) << std::endl;
    result << accessRelSS.str() << std::endl;

    accessRelSS.clear();
    accessRelSS.str("");
    accessRelNumRows  = 0;

    result << "# Read items" << std::endl;

    for (MemoryAccess *MA : *Stmt) {

      if (MA->isRead()) {
        accessRelNumRows
            += outputAccessRelation(MA, arrayId2String, &p, accessRelSS);
      }
    }
    result << accessRelNumRows << " " << accessRelNumCols << std::endl;
    result << accessRelSS.str() << std::endl;
  }

  result << blockDelimiter << std::endl
         << "# Transformation Candidate" << std::endl << "0" << std::endl;

  isl_printer_free(p);

  return result.str();
}

bool LetseeExporter::runOnScop(Scop &S) {
  Region &R = S.getRegion();
  std::string fileName = getFileName(S);
  std::string functionName = R.getEntry()->getParent()->getName();
  errs() << "Writing SCoP '" << R.getNameStr() << "' in function '"
         << functionName << "' to '" << fileName << "'.\n";

  std::string fileContent = getCloogInput(S);

  // Write to file.
  std::error_code EC;
  tool_output_file F(fileName, EC, llvm::sys::fs::F_Text);
  if (!EC) {
    F.os() << fileContent;
    F.os().close();

    if (!F.os().has_error()) {
      errs() << "\n";
      F.keep();
      return false;
    }
  }
  errs() << "error opening file for writing.\n";
  F.os().clear_error();

  return false;
}

void LetseeExporter::getAnalysisUsage(AnalysisUsage &AU) const {
  AU.setPreservesAll();
  AU.addRequired<ScopInfo>();
}

Pass *polly::createLetseeExporterPass() { return new LetseeExporter(); }

INITIALIZE_PASS_BEGIN(LetseeExporter, "polly-export-letsee",
                      "Polly - Export Scops as Letsee input files"
                      " (Writes a .candl file for each Scop)",
                      false, false);
INITIALIZE_PASS_DEPENDENCY(ScopInfo)
INITIALIZE_PASS_END(LetseeExporter, "polly-export-letsee",
                    "Polly - Export Scops as Letsee input files"
                    " (Writes a .candl file for each Scop)",
                    false, false)
