#ifndef RANKMASK_HXX
#define RANKMASK_HXX

#include <vector>

namespace mota {

struct RankMask {
  const std::vector<int> &node_ids;
  int rank_n;

  RankMask(const std::vector<int> &n, int r) : node_ids(n), rank_n(r) {}
};

}

#endif // RANKMASK_HXX
