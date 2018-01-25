#ifndef RANKWEIGHT_HXX
#define RANKWEIGHT_HXX

namespace mota {

// W is a weight class that provides scalarization and ordering
template <class W>
struct RankWeight {
  int rank;
  W weight;
};

// reset all weights to w
template <class W>
void reset_weights(std::vector<RankWeight<W>> &vec, W w) {
  for (RankWeight<W> &rw : vec) rw.weight = w;
}

// sort by increasing scalar weight
template <class W>
void sort_by_scalar_weight(std::vector<RankWeight<W>> &vec) {
  std::sort(vec.begin(), vec.end(),
    [](const RankWeight<W> &a, const RankWeight<W> &b) {
      return static_cast<double>(a.weight) <
             static_cast<double>(b.weight);
  });
}

template <class W, class U>
bool rw_weight_gt(const RankWeight<W> &a, const RankWeight<W> &b) {
  return static_cast<U>(a.weight) > static_cast<U>(b.weight);
};

template <class W, class U>
bool rw_weight_lt(const RankWeight<W> &a, const RankWeight<W> &b) {
  return static_cast<U>(a.weight) < static_cast<U>(b.weight);
};

template <class W>
struct RankWeightCap {
  int rank;
  W weight;
  W cap;
};

// reset all weights to w
template <class W>
void reset_weights(std::vector<RankWeightCap<W>> &vec, W w) {
  for (RankWeightCap<W> &rwc : vec) rwc.weight = w;
}

}

#endif // RANKWEIGHT_HXX
