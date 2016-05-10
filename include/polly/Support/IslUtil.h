#ifndef ISLUTIL_H
#define ISLUTIL_H

#include "isl/constraint.h"
#include "isl/map.h"
#include "isl/union_map.h"
#include "isl/printer.h"
#include "isl/set.h"
#include "isl/union_set.h"
#include "isl/schedule.h"
#include "isl/schedule_node.h"
#include "isl/schedule_type.h"
#include <memory>
#include <string>
#include <system_error>
#include "llvm/ADT/DenseSet.h"

namespace polly {

isl_stat foreachSetHelper(isl_set *s, void *user);

void callLambda(__isl_keep isl_union_set *s,
                       std::function<void(isl_set*)> &f);

isl_stat foreachMapHelper(isl_map *m, void *user);

void callLambda(__isl_keep isl_union_map *m,
                       std::function<void(isl_map*)> &f);

__isl_give isl_schedule_node *scheduleMapHelper(isl_schedule_node *n,
                                                       void *user);

__isl_give isl_schedule *callLambda(__isl_take isl_schedule *s,
                std :: function<isl_schedule_node*(isl_schedule_node*)> &f);

llvm::DenseSet<const char *> getTupleNamesFromUnionSet(
    __isl_keep isl_union_set *s);

}

#endif // ISLUTIL_H
