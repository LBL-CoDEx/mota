#ifndef NETMODELS_HXX
#define NETMODELS_HXX

#include "amrweight.hxx"
#include "netgraph.hxx"
#include "cartesian.hxx"
#include "morton.hxx"
#include "flags.hxx"
#include "random.hxx"

namespace mota {

const bool flag_verbose_nodes = false;

template <int N>
using RanksByCartesian = std::unordered_map<Cartesian<int,N>, std::vector<int>>;

// helper to make a simple rank map
template <int N>
RanksByCartesian<N> make_ranks_by_car(
    Cartesian<int,N> topo,
    int rank_n,
    int node_rank_n,
    bool flag_rand_place)
{
  auto max_rank_n = prod(topo) * node_rank_n;
  if (rank_n > max_rank_n) {
    Say() << "WARNING: decreasing rank_n from " << rank_n
          << " to " << max_rank_n << " to make job fit!";
    rank_n = max_rank_n;
  }
  RanksByCartesian<N> ranks_by_car{};
  int cur_rank = 0;
  if (flag_rand_place) {
    // assign ranks to nodes
    std::uniform_int_distribution<int> rnd_dist(0, prod(topo)-1);
    while (cur_rank < rank_n) {
      Cartesian<int,N> pt(random::roll(rnd_dist), topo);
      pt[0]++; // add switch
      // check for space (implicit insert of pt into ranks_by_car)
      if (ranks_by_car[pt].size() < static_cast<unsigned>(node_rank_n)) {
        // assign rank to node
        if (flag_verbose_nodes) Say() << "assigning rank " << cur_rank << " to node " << pt;
        ranks_by_car[pt].push_back(cur_rank);
        ++cur_rank;
      }
    }
  } else {
    for (Cartesian<int,N> pt : CartesianProd<int,N>(topo)) {
      pt[0]++; // add switch
      if (cur_rank >= rank_n) break;
      std::vector<int> rank_vec{};
      for (int i = 0; i < node_rank_n; ++i) {
        if (cur_rank >= rank_n) break;
        rank_vec.push_back(cur_rank);
        ++cur_rank;
      }
      assert(rank_vec.size() > 0);
      ranks_by_car[pt] = std::move(rank_vec);
    }
  }
  return ranks_by_car;
}

/**************************
 * ND-Torus Network Graph *
 **************************/

// NOTE: N is the dimension of the torus.
// The first dimension of Coord determines whether the node is a switch or compute node.
//   First dimension: (0: switch, >0: compute node)
// The next N coordinates specify shape of torus,
template <class W, long unsigned N>
class TorusNDGraph : public NetworkGraph<W> {
private:
  using Coord = Cartesian<int,N+1>;
  using RanksByCoord = RanksByCartesian<N+1>;
  static const unsigned kNodeDim = 0; // node dim in Coord
  static const unsigned kTorDimLo = 1;
  static const unsigned kTorDimHi = N+1;

public:

  TorusNDGraph() = default;
  virtual ~TorusNDGraph() = default;

  // Coord shape: first dimension is number of compute nodes per network endpoint
  //              next N dimensions specify shape of torus
  // Coord link_widths: link widths (corresponding to topology_redundant in SST/macro)
  // rank_n: instantiate only enough compute nodes to fit rank_n ranks;
  //         still use full torus of switches specified by shape (for routing)
  // node_rank_n: number of ranks per node
  TorusNDGraph(std::string label, Coord shape,
               const Coord &link_widths,
               const RanksByCoord &ranks_by_coord)
    : label_(std::move(label)),
      shape_(std::move(shape))
  {
    shape_[kNodeDim]  += 1; // add a switch

    // make nodes and switches
    for (Coord pt : CartesianProd<int,N+1>(shape_)) {
      bool has_ranks = ranks_by_coord.count(pt);
      if (is_switch(pt) || has_ranks) {
        if (flag_verbose_network && is_switch(pt)) Say() << "Constructing " << pt << " ...";
        NetworkNode<W> &n = this->get_node_create(get_id(pt));
        if (has_ranks) {
          for (auto r : ranks_by_coord.at(pt)) { this->add_rank(n, r); }
        }
        n.edges = make_edges(pt, link_widths, ranks_by_coord);
      }
    }

    if (!flag_quiet) {
      // omit the switch in description
      auto temp = shape_; temp[kNodeDim]--;
      Say() << "\nBuilding ND Torus:";
      Say() << "  Physical       : " << str(temp);
      Say() << "  Ranks          : " << this->rank_n();
      Say() << "  Compute Nodes  : " << this->compute_n();
      Say() << "  Switches       : " << this->switch_n();
      Say() << "  Edges          : " << this->edge_n();
      Say() << "  Redundant links: " << str(link_widths);
      Say() << "  Links          : " << this->link_n();
    }
  }

  const std::string label() const override { return label_; }

  std::string node_id_str(int node_id) const override {
    return Coord(node_id, shape_).str();
  }

private:
  std::string label_;
  Coord shape_;     // torus network shape
  int node_rank_n_; // ranks per node

  static Coord no_switch(Coord x) {
    assert(x[kNodeDim] > 0);
    Coord result{std::move(x)};
    result[kNodeDim]--; // remove switch
    return result;
  }

  int get_id(const Coord &x) const { return x.shape_id(shape_); }
  static bool is_switch(const Coord &x) { return x[kNodeDim] == 0; }
  bool is_switch(int id) const override { return is_switch(Coord{id, shape_}); }
  static Coord get_switch(Coord x) { x[kNodeDim] = 0; return x; }
  int get_switch(int id) const override { return get_id(get_switch(Coord{id, shape_})); }

  // how to link nodes
  using Edge = typename NetworkGraph<W>::Edge;
  using EdgeMap = std::unordered_map<int, Edge>;
  EdgeMap make_edges(const Coord &x, const Coord &link_widths, const RanksByCoord &ranks_by_coord) {
    EdgeMap result;
    if (is_switch(x)) {
      // connect to neighbor switches
      for (unsigned d = kTorDimLo; d < kTorDimHi; ++d) {
        Coord offset(0);
        for (int s = -1; s <= 1; s += 2) {
          offset[d] = s;
          Coord dst = (x + offset) % shape_;
          if (dst == x) continue;
          if (flag_verbose_network) Say() << "  switch link: " << dst;
          result.insert({get_id(dst), Edge{link_widths[d], false}});
        }
      }
      // connect to compute nodes
      for (int i = 1; i < shape_[kNodeDim]; ++i) {
        Coord dst{x};
        dst[kNodeDim] = i;
        if (ranks_by_coord.count(dst)) {
          if (flag_verbose_network) Say() << "  node link: " << dst;
          result.insert({get_id(dst), Edge{link_widths[kNodeDim], true}});
        }
      }
    } else {
      Coord dst = get_switch(x);
      result.insert({get_id(dst), Edge{link_widths[kNodeDim], true}});
    }
    return result;
  }

  // dimension-ordered routing
  using Path = typename NetworkGraph<W>::Path;
  void switch_path(Path &path, int start_id, int end_id) const override {
    path.push_back(start_id);
    Coord here{start_id, shape_}, end{end_id, shape_};
    for (unsigned i = kTorDimLo; i < kTorDimHi; ++i) {
      while (here[i] != end[i]) {
        bool go_up = wrap(end[i] - here[i], shape_[i]) <= shape_[i] / 2;
        here[i] = wrap(here[i] + (go_up ? +1 : -1), shape_[i]);
        path.push_back(get_id(here));
      }
    }
    assert(here == end);
    return;
  }

  
  // number of hops from start to end
  int hop_dist(int src_node_id, int dst_node_id) const override {
    Coord src(src_node_id, shape_), dst(dst_node_id, shape_);
    int result = 0;
    if (src != dst) {
      if (src[kNodeDim] != 0) { result++; }
      for (unsigned i = kTorDimLo; i < kTorDimHi; ++i) {
        auto up_dist = wrap(dst[i] - src[i], shape_[i]);
        bool go_up = (up_dist <= shape_[i] / 2);
        result += (go_up ? up_dist : shape_[i] - up_dist);
      }
      if (dst[kNodeDim] != 0) { result++; }
    }
    return result;
  }

  template <class C>
  static bool zmorton_lt(int lo, int hi, const C &a, const C &b) {
    int dim = hi-lo;
    switch (dim) {
      case 2:
        return EncodeMorton2(a[lo], a[lo+1]) <
               EncodeMorton2(b[lo], b[lo+1]);
        break;
      case 3:
        return EncodeMorton3(a[lo], a[lo+1], a[lo+2]) <
               EncodeMorton3(b[lo], b[lo+1], b[lo+2]);
        break;
      case 4:
        return EncodeMorton4(a[lo], a[lo+1], a[lo+2], a[lo+3]) <
               EncodeMorton4(b[lo], b[lo+1], b[lo+2], b[lo+3]);
        break;
      case 5:
        return EncodeMorton5(a[lo], a[lo+1], a[lo+2], a[lo+3], a[lo+4]) <
               EncodeMorton5(b[lo], b[lo+1], b[lo+2], b[lo+3], b[lo+4]);
        break;
      default:
        // TODO: implement general case
        abort();
    }
    return false;
  }

  // override for id_weights to return an SFC sorted list
  virtual std::vector<RankWeight<W>> rank_weights() const override {
    auto result = NetworkGraph<W>::rank_weights();
    // sort ranks with space-filling curve
    std::sort(result.begin(), result.end(),
      [=](const RankWeight<W> &a, const RankWeight<W> &b) {
        int a_net_id = this->node_at_rank(a.rank).id;
        int b_net_id = this->node_at_rank(b.rank).id;
        return zmorton_lt(kTorDimLo, kTorDimHi,
                          Coord{a_net_id, shape_},
                          Coord{b_net_id, shape_});
    });
    return result;
  }
};

// Template specializations for ND torus

template <class W>
using Torus2DGraph = TorusNDGraph<W,2>;
template <class W>
using Torus3DGraph = TorusNDGraph<W,3>;
template <class W>
using Torus4DGraph = TorusNDGraph<W,4>;
template <class W>
using Torus5DGraph = TorusNDGraph<W,5>;


/**************************
 * Dragonfly Network Graph *
 **************************/

// NOTE: N is the dimension of each dragonfly group (each is a hypercube).
// The first dimension of Coord determines whether the node is a switch or compute node.
//   First dimension: (0: switch, >0: compute node)
// The next N coordinates specify shape of each hypercube group.
// The final dimension is the group number.
template <class W, long unsigned N>
class DragonflyGraph : public NetworkGraph<W> {
private:
  using Coord = Cartesian<int,N+2>;
  using HCCoord = Cartesian<int,N>;
  using RanksByCoord = RanksByCartesian<N+2>;
  static const unsigned kNodeDim = 0; // node dim in Coord
  static const unsigned kHCDimLo = 1;
  static const unsigned kHCDimHi = N+1;
  static const unsigned kGroupDim = N+1; // group dim in Coord

  using rand_engine_t = std::default_random_engine;

public:
  // Coord s: first dimension is number of compute nodes per network endpoint
  //          next N dimensions specify shape of hypercube group
  //          final dimension is the group number
  // inter_link_n is number of inter-group links per switch
  // link_widths: link widths along each dimension (# of redundant links)
  DragonflyGraph() = default;
  DragonflyGraph(const DragonflyGraph &) = default;
  DragonflyGraph(DragonflyGraph &&) = default;
  virtual ~DragonflyGraph() = default;
  DragonflyGraph &operator=(const DragonflyGraph &) = default;
  DragonflyGraph &operator=(DragonflyGraph &&) = default;

  DragonflyGraph(std::string label, Coord s, Coord link_widths,
                 int iln, const RanksByCoord &ranks_by_coord) {
    init(label, s, link_widths, iln, ranks_by_coord);
  }

  void init(std::string label, Coord s, Coord link_widths,
            int iln, const RanksByCoord &ranks_by_coord) {
    _label = label;
    _shape = s;
    _shape[kNodeDim] += 1; // add a switch to node dimension
    inter_link_n = iln;
    if (HC_size() * inter_link_n < (unsigned) group_n()) {
      Say() << "Dragonfly needs more intergroup links";
      abort();
    }

    // make nodes and switches
    for (Coord pt : CartesianProd<int,N+2>(_shape)) {
      bool has_ranks = ranks_by_coord.count(pt);
      if (is_switch(pt) || has_ranks) {
        if (flag_verbose_network && is_switch(pt)) Say() << "Constructing " << pt << " ...";
        NetworkNode<W> &n = this->get_node_create(get_id(pt));
        if (has_ranks) {
          for (auto r : ranks_by_coord.at(pt)) { this->add_rank(n, r); }
        }
        n.edges = make_edges(pt, link_widths, ranks_by_coord);
      }
    }

    { // omit the switch in description
      auto temp = _shape; temp[kNodeDim]--;
      if (!flag_quiet) {
        Say() << "\nBuilding Dragonfly:";
        Say() << "  Physical       : " << str(temp);
        Say() << "  Ranks          : " << this->rank_n();
        Say() << "  Compute Nodes  : " << this->compute_n();
        Say() << "  Switches       : " << this->switch_n();
        Say() << "  Edges          : " << this->edge_n();
        Say() << "  Redundant links: " << str(link_widths);
        Say() << "  Links          : " << this->link_n();

        if (flag_verbose_nodes) {
          Say() << "Nodes with ranks:";
          this->visit([&](const NetworkNode<W> &n) {
            if (!is_switch(n.id)) {
              Say() << "  " << n;
            }
          });
        }
      }
    }
  }

  const std::string label() const override { return _label; }
  int group_n() const { return _shape[kGroupDim]; }

  std::string node_id_str(int node_id) const override {
    return Coord(node_id, _shape).str();
  }

private:
  std::string _label;
  Coord _shape;
  int inter_link_n = 0;
  mutable unsigned _group_stride = 0;

  static Coord no_switch(Coord x) {
    assert(x[kNodeDim] > 0);
    Coord result{std::move(x)};
    result[kNodeDim]--; // remove switch
    return result;
  }

  int get_id(const Coord &x) const { return x.shape_id(_shape); }
  static bool is_switch(const Coord &x) { return x[kNodeDim] == 0; }
  bool is_switch(int id) const override { return is_switch(Coord{id, _shape}); }
  static Coord get_switch(Coord x) { x[kNodeDim] = 0; return x; }
  int get_switch(int id) const override { return get_id(get_switch(Coord{id, _shape})); }

  unsigned group_stride() const {
    if (_group_stride == 0) { 
      _group_stride = std::max(1, group_n() / inter_link_n);
    }
    return _group_stride;
  }

  HCCoord HC_coords(Coord x) const {
    HCCoord result;
    for (unsigned i = 0, d = kHCDimLo; d < kHCDimHi; ++i, ++d) result[i] = x[d];
    return result;
  }

  size_t HC_size() const {
    size_t result = 1;
    for (unsigned d = kHCDimLo; d < kHCDimHi; ++d) result *= _shape[d];
    return result;
  }

  int group_id(const Coord &x) const {
    return HC_coords(x).shape_id(HC_coords(_shape));
  }
  Coord coord_from_group_id(int group_id, int group) const {
    HCCoord hc_coord{group_id, HC_coords(_shape)};
    Coord result;
    result[kNodeDim] = 0; // switch
    for (unsigned i = 0, d = kHCDimLo; d < kHCDimHi; ++i, ++d) {
      result[d] = hc_coord[i];
    }
    result[kGroupDim] = group;
    return result;
  }

  // how to link nodes
  using Edge = typename NetworkNode<W>::Edge;
  using EdgeMap = std::unordered_map<int, Edge>;
  using Path = typename NetworkGraph<W>::Path;

  EdgeMap make_edges(const Coord &x, const Coord &link_widths, const RanksByCoord &ranks_by_coord) {
    EdgeMap result;
    if (is_switch(x)) {
      // connect within group
      for (unsigned d = kHCDimLo; d < kHCDimHi; ++d) {
        // group hypercube is fully connected along each dimension
        for (int i = 0; i < _shape[d]; ++i) {
          if (i == x[d]) continue; // not to itself
          Coord dst{x};
          dst[d] = i;
          if (flag_verbose_network) Say() << "  switch link: " << dst;
          result.insert({get_id(dst), Edge{link_widths[d], false}});
        }
      }
      // connect across groups
      for (int k = 0; k < inter_link_n; ++k) {
        Coord dst{x};
        dst[kGroupDim] = wrap(x[kGroupDim] + group_id(x) + k*group_stride(),
                              group_n());
        if (dst[kGroupDim] == x[kGroupDim]) continue; // don't self-connect
        if (flag_verbose_network) Say() << "  switch link: " << dst;
        result.insert({get_id(dst), Edge{link_widths[kGroupDim], false}});
      }
      // connect to compute nodes
      for (int i = 1; i < _shape[kNodeDim]; ++i) {
        Coord dst{x};
        dst[kNodeDim] = i;
        if (ranks_by_coord.count(dst)) {
          if (flag_verbose_network) Say() << "  node link: " << dst;
          result.insert({get_id(dst), Edge{link_widths[kNodeDim], true}});
        }
      }
    } else {
      Coord dst = get_switch(x);
      result.insert({get_id(dst), Edge{link_widths[kNodeDim], true}});
    }
    return result;
  }

  bool has_group_connection(Coord x, int g) const {
    int group_offset = wrap(g-(x[kGroupDim]+group_id(x)), group_n());
    return (group_offset % group_stride() == 0 &&
            group_offset / group_stride() < (unsigned) inter_link_n);
  }

  // assumes all groups are directly reachable from some switch in group
  Coord find_group_intermediate(Coord start, Coord end) const {
    // check current switch for inter-group connection
    int start_group = start[kGroupDim], end_group = end[kGroupDim];
    if (has_group_connection(start, end_group)) {
      return start;
    }

    // check distance 1 switches
    for (unsigned i = kHCDimLo; i < kHCDimHi; ++i) {
      std::uniform_int_distribution<int> rnd_dist(0, _shape[i]-1); // inclusive bounds
      int k_start = random::roll(rnd_dist); // pick a random start index
      Coord cur{start};
      for (int k = 0; k < _shape[i]; ++k) {
        cur[i] = (k_start+k) % _shape[i];
        if (has_group_connection(cur, end_group)) {
          return cur;
        }
      }
    }

    // find first (by group id) switch in start group connected to end group
    // this is what SST/macro does
    int group_offset = wrap(end_group-start_group, group_n());
    int itm_group_id = group_offset % group_stride();
    if (group_offset / group_stride() >= (unsigned) inter_link_n) {
      itm_group_id += group_stride();
    }
    Coord itm = coord_from_group_id(itm_group_id, start_group);
    assert(has_group_connection(itm, end_group));
    return itm;
  }

  // number of intragroup hops from start to end switches
  int intragroup_hop_dist(Coord src, Coord dst) const {
    assert(src[kGroupDim] == dst[kGroupDim]);
    int result = 0;
    if (src != dst) {
      if (src[kNodeDim] != 0) { result++; }
      for (unsigned i = kHCDimLo; i < kHCDimHi; ++i) {
        if (src[i] != dst[i]) {
          result++;
        }
      }
      if (dst[kNodeDim] != 0) { result++; }
    }
    return result;
  }

  // number of hops from start to end
  int hop_dist(int src_node_id, int dst_node_id) const override {
    Coord src(src_node_id, _shape), dst(dst_node_id, _shape), here(src);
    int result = 0;
    if (src[kGroupDim] != dst[kGroupDim]) {
      Coord itm = find_group_intermediate(src, dst);
      // go from src to src_itm to dst_itm
      result += intragroup_hop_dist(src, itm) + 1;
      itm[kGroupDim] = dst[kGroupDim];
      here = itm;
    }
    result += intragroup_hop_dist(here, dst);
    return result;
  }

  // dimension-ordered routing
  void push_intragroup_path(Path &path, Coord here, Coord end) const {
    path.push_back(get_id(here));
    if (flag_verbose_routing) Say() << "    Traversing " << here;
    for (unsigned i = kHCDimLo; i < kHCDimHi; ++i) {
      if (here[i] != end[i]) {
        here[i] = end[i];
        path.push_back(get_id(here));
        if (flag_verbose_routing) Say() << "    Traversing " << here;
      }
    }
  }

  // find a path from start to end
  void switch_path(Path &path, int start_id, int end_id) const override {
    Coord here{start_id, _shape}, end{end_id, _shape};
    if (flag_verbose_routing) Say() << "Switch path from " << here << " to " << end;
    if (here[kGroupDim] != end[kGroupDim]) {
      // go to intermediate switch
      if (flag_verbose_routing) Say() << "  Finding intermediate switch ...";
      Coord itm = find_group_intermediate(here, end);
      push_intragroup_path(path, here, itm);
      // go to destination group
      itm[kGroupDim] = end[kGroupDim];
      here = itm;
    }
    // go to end switch
    if (flag_verbose_routing) Say() << "  Going to destination switch ...";
    push_intragroup_path(path, here, end);
    return;
  }
};

#ifdef WITH_ARIES
#include "netmodels-aries.hxx"
#else

struct NetParams {
  // common
  Car4 topo;
  Car4 redundant_link_n;
  // for dragonfly
  int group_n;
  int intergroup_link_n; // generic dfly
};

inline NetParams get_net_params(std::string label) {
  if (label == "3dt-edi") {
    return NetParams {Int4{4, 16, 6, 15},
                      Int4{2,  6, 6,  8},  0, 0};
  } else if (label == "3dt-exa") {
    return NetParams {Int4{4, 16, 6, 50},
                      Int4{2,  8, 8,  8},  0, 0};
  } else if (label == "dfly-edi") {
    return NetParams {Int4{4, 16, 6, 15},
                      Int4{2,  1, 3,  2}, 15, 5};
  } else if (label == "dfly-exa") {
    return NetParams {Int4{4, 16, 6, 50},
                      Int4{4,  1, 3,  2}, 50, 9};
  } else {
    Say() << "unsupported network: " << label;
    std::abort();
  }
}

template <class W>
std::shared_ptr<NetworkGraph_>
net_graph_create_helper(std::string label, const RanksByCartesian<4> &ranks_by_coord) {
  auto params = get_net_params(label);
  if (       label == "3dt-edi" ||
             label == "3dt-exa") {
    return std::shared_ptr<NetworkGraph_>(
      new Torus3DGraph<W>(label, params.topo, params.redundant_link_n,
                          ranks_by_coord));
  } else if (label == "dfly-edi" ||
             label == "dfly-exa") {
    return std::shared_ptr<NetworkGraph_>(
      new DragonflyGraph<W,2>(label, params.topo, params.redundant_link_n,
                              params.intergroup_link_n, ranks_by_coord));
  } else {
    Say() << "unsupported network: " << label;
    std::abort();
  }
}

#endif

// create an application graph with weights for the right number of AMR levels
// TODO: there should be a better way to do this
inline std::shared_ptr<NetworkGraph_>
net_graph_create(std::string label, int amr_level_n,
                 const RanksByCartesian<4> &ranks_by_coord) {
  random::init(0);
  std::shared_ptr<NetworkGraph_> result {nullptr};
  switch (amr_level_n) {
    case 1:
      result = net_graph_create_helper<AMRWeight<1>>(label, ranks_by_coord);
      break;
    case 2:
      result = net_graph_create_helper<AMRWeight<2>>(label, ranks_by_coord);
      break;
    case 3:
      result = net_graph_create_helper<AMRWeight<3>>(label, ranks_by_coord);
      break;
    case 4:
      result = net_graph_create_helper<AMRWeight<4>>(label, ranks_by_coord);
      break;
    case 5:
      result = net_graph_create_helper<AMRWeight<5>>(label, ranks_by_coord);
      break;
    default:
      Say() << "ERROR: Number of AMR levels (" << amr_level_n << ") is over maximum supported (11)!";
      std::abort();
      break;
  }
  return result;
}

inline std::shared_ptr<NetworkGraph_>
net_graph_create(std::string label, int amr_level_n,
                 size_t rank_n, int node_rank_n, bool flag_rand_place) {
  random::init(0);
  auto ranks_by_car = make_ranks_by_car<4>(get_net_params(label).topo,
                                           rank_n, node_rank_n, flag_rand_place);
  return net_graph_create(label, amr_level_n, ranks_by_car);
}

}

#endif // NETMODELS_HXX
