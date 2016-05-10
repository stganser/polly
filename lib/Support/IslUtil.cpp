#include "polly/Support/IslUtil.h"

using namespace llvm;
using namespace polly;

isl_stat polly::foreachSetHelper(isl_set *s, void *user) {
  std::function<void(isl_set*)> *f = (std::function<void(isl_set*)> *) user;
  (*f)(s);
  return isl_stat_ok;
}

void polly::callLambda(__isl_keep isl_union_set *s,
                              std::function<void(isl_set*)> &f) {
  isl_union_set_foreach_set(s, foreachSetHelper, &f);
}

isl_stat polly::foreachMapHelper(isl_map *m, void *user) {
  std::function<void(isl_map*)> *f = (std::function<void(isl_map*)> *) user;
  (*f)(m);
  return isl_stat_ok;
}

void polly::callLambda(__isl_keep isl_union_map *m,
                              std::function<void(isl_map*)> &f) {
  isl_union_map_foreach_map(m, foreachMapHelper, &f);
}

__isl_give isl_schedule_node *polly::scheduleMapHelper(
    isl_schedule_node *n, void *user) {
  std::function<isl_schedule_node*(isl_schedule_node*)> *f
      = (std::function<isl_schedule_node*(isl_schedule_node*)> *) user;
  return (*f)(n);
}

__isl_give isl_schedule *polly::callLambda(__isl_take isl_schedule *s,
                std :: function<isl_schedule_node*(isl_schedule_node*)> &f) {
  return isl_schedule_map_schedule_node_bottom_up(s, scheduleMapHelper, &f);
}

DenseSet<const char *> polly::getTupleNamesFromUnionSet(
    __isl_keep isl_union_set *s) {
  DenseSet<const char *> names;
  std::function<void(isl_set *)> lambda = [&names](isl_set *s) {
    names.insert(isl_set_get_tuple_name(s));
    isl_set_free(s);
  };
  callLambda(s, lambda);
  return names;
}
