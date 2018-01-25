#ifndef APPGRAPH_HXX
#define APPGRAPH_HXX

#include <algorithm>
#include <vector>
#include <string>
#include <memory>

#include "graph.hxx"
#include "cartesian.hxx"
#include "random.hxx"
#include "rankweight.hxx"
#include "rankmask.hxx"
#include "morton.hxx"
#include "flags.hxx"
#include "rcm.h"
#if KNOB_METIS
#  include <metis.h>
#endif

namespace mota {

// tolerance for compute load imbalance, used for GR and RB mappers
const double imba_ratio = 1.01;

// forward declaration
template <class W> class NetworkGraph;

/*********************
 * Application Graph *
 *********************/

template <class W>
struct AppNode {
  using Weight = W;
  using Ptr = std::shared_ptr<AppNode>;

  struct Edge {
    using Weight = size_t;
    int to_msg_n = 0, from_msg_n = 0;
    size_t to_byte_n = 0, from_byte_n = 0;
    // TODO: take into account message count
    size_t weight() const { return to_byte_n + from_byte_n; }
  };

  int id;
  double load;
  Cartesian<int,3> geom_loc;
  size_t cell_n;
  int level;
  int rank;   // job rank mapping assignment
  std::unordered_map<int, Edge> edges;

  explicit AppNode(int _id = kNone) :
    id{_id}, load{0.0}, geom_loc{kNone}, cell_n{0},
    level{kNone}, rank{kNone}, edges{} {}

  W weight() const;
  double scalar_weight() const { return static_cast<double>(weight()); }

  friend std::ostream &operator<<(std::ostream &os, const AppNode &n) {
    os << "AppNode (" << n.id << ", " << n.geom_loc << ": " << n.cell_n << ", "
       << n.weight() << ") -> " << n.rank;
    return os;
  }

  size_t sum_edge_weight() const {
    size_t result = 0;
    for (const auto &p : edges) result += p.second.weight();
    return result;
  }
};

template <class W>
W AppNode<W>::weight() const {
  return W{load, level};
}
template <>
inline double AppNode<double>::weight() const {
  return load;
}

// generic interface class to AppGraph<W>
class AppGraph_ {
public:
  AppGraph_() = default;
  virtual ~AppGraph_() = default;

  virtual void set_node_loc(int id, Cartesian<int,3> loc) = 0;
  virtual void set_node_lev(int id, int lev) = 0;
  virtual void set_node_cell_n(int id, size_t cell_n) = 0;

  virtual void reset_comp_comm() = 0;
  virtual void import(const std::unordered_map<int, double> &comps,
                      const std::unordered_map<std::pair<int, int>, std::pair<int, size_t>> &comms) = 0;
  virtual void print_stats() const = 0;

  virtual void init_rank_masks(std::shared_ptr<NetworkGraph_> net_g,
                               const std::vector<int> &map_rank_n) = 0;
  virtual void run_mapper(std::string mapper_str,
                          std::shared_ptr<NetworkGraph_> net_g) = 0;

  // vector of available mappers
  virtual std::vector<std::string> available_mappers() const = 0;
  // mapping from app node ID to network node ID
  virtual std::unordered_map<int, int> get_mapping() const = 0;
  virtual void clear_mapping() = 0;
  virtual void write_mapping(std::string filename) const = 0;
  virtual void evaluate_mapping(std::shared_ptr<NetworkGraph_> net_g,
                                const std::unordered_map<int, int> &mapping) = 0;
};

template <class W>
class AppGraph : public Graph<AppNode<W>>, public AppGraph_ {
public:
  using NodePtr = typename AppNode<W>::Ptr;

  AppGraph() = default;
  virtual ~AppGraph() = default;

  virtual void set_node_loc(int id, Cartesian<int,3> loc) override {
    this->get_node_create(id).geom_loc = loc;
  }
  virtual void set_node_lev(int id, int lev) override {
    this->get_node_create(id).level = lev;
  }
  virtual void set_node_cell_n(int id, size_t cell_n) override {
    this->get_node_create(id).cell_n = cell_n;
  }

  virtual void reset_comp_comm() override;
  virtual void import(const std::unordered_map<int, double> &comps,
                      const std::unordered_map<std::pair<int, int>, std::pair<int, size_t>> &comms) override;
  virtual void print_stats() const override;

  virtual void init_rank_masks(std::shared_ptr<NetworkGraph_> net_g,
                               const std::vector<int> &map_rank_n) override;

  virtual void run_mapper(std::string mapper_str,
                          std::shared_ptr<NetworkGraph_> net_g) override; // dispatch

  // vector of available mappers
  std::vector<std::string> available_mappers() const { return mapper_strs(); }
  // mapping from app node ID to network node ID
  std::unordered_map<int, int> get_mapping() const {
    std::unordered_map<int, int> result;
    this->visit([&](const AppNode<W> &n) { result[n.id] = n.rank; });
    return result;
  }
  void clear_mapping() { this->visit([](AppNode<W> &n) { n.rank = kNone; }); }
  void write_mapping(std::string filename) const;
  void evaluate_mapping(std::shared_ptr<NetworkGraph_> net_g,
                        const std::unordered_map<int, int> &mapping);

private:

  std::vector<RankMask> rank_masks_; // one rank mask per level
  int cached_max_node_level_ = std::numeric_limits<int>::min();

  // helpers for sorting nodes

  // returns True if a is heavier than b
  static bool node_scalar_weight_gt(const AppNode<W> &a, const AppNode<W> &b) {
    return a.scalar_weight() > b.scalar_weight();
  }
  // returns True if a < b in Z-Morton ordering
  static bool zmorton_lt(Cartesian<int,3> a, Cartesian<int,3> b) {
    return EncodeMorton3(a[0], a[1], a[2]) < EncodeMorton3(b[0], b[1], b[2]);
  }
  static bool node_zmorton_lt(const AppNode<W> &a, const AppNode<W> &b) {
    return zmorton_lt(a.geom_loc, b.geom_loc);
  }

  // inspectors
  int max_node_level() const {
    int max_level = cached_max_node_level_;
    if (max_level == std::numeric_limits<int>::min()) {
      this->visit([&](const AppNode<W> &n) { max_level = std::max(max_level, n.level); });
    }
    return max_level;
  }

  // filter the nodes by level
  using Graph<AppNode<W>>::filtered_nodes;
  std::vector<NodePtr> nodes_in_level(int level) {
    return this->filtered_nodes([=](const AppNode<W> &n) { return n.level == level; });
  }

  // filter the nodes by level and sort
  using Graph<AppNode<W>>::sorted_nodes;
  template <class CompF>
  std::vector<NodePtr> sorted_nodes_in_level(CompF comp_f, int level) {
    return this->sorted_nodes(comp_f, [=](const AppNode<W> &n) { return n.level == level; });
  }

  // mapper functions
  void map_random(NetworkGraph<W> &net_g);
  void map_round_robin(NetworkGraph<W> &net_g);
  void map_knapsack(NetworkGraph<W> &net_g);
  void map_sfc(NetworkGraph<W> &net_g,
               bool flag_per_lev,     // map levels one at a time
               bool flag_remap_levs,  // remap levels for load balance
               std::string tag);
  void map_sfcs(NetworkGraph<W> &net_g) { map_sfc(net_g, true, true , "SFCS"); }
  void map_pfcs(NetworkGraph<W> &net_g) { map_sfc(net_g, true, false, "PFCS"); }
  void map_pfcm(NetworkGraph<W> &net_g) { map_sfc(net_g, false, true, "PFCM"); }
  void map_greedy(NetworkGraph<W> &net_g);
  void map_rcm(NetworkGraph<W> &net_g, bool flag_pfcm = false);
  void map_rcm_rcm (NetworkGraph<W> &net_g) { map_rcm(net_g, false); }
  void map_rcm_pfcm(NetworkGraph<W> &net_g) { map_rcm(net_g, true ); }
#if KNOB_METIS
  void map_rb(NetworkGraph<W> &net_g, bool flag_multicons, bool flag_pfcm = false);
  void map_rb_rb  (NetworkGraph<W> &net_g) { map_rb(net_g, false, false); }
  void map_rb_rbm (NetworkGraph<W> &net_g) { map_rb(net_g, true, false); }
  void map_rb_pfcm(NetworkGraph<W> &net_g) { map_rb(net_g, false, true ); }
#endif

  std::vector<std::string> mapper_strs() const {
    return std::vector<std::string> {
      "rdm", "rr", "ks", "sfcs", "pfcs", "pfcm", "gr", "rcm", "rcm-pfcm",
#if KNOB_METIS
      "rb", "rbm", "rb-pfcm",
#endif
    };
  }

  using MapperFn = void (AppGraph<W>::*)(NetworkGraph<W> &);
  std::unordered_map<std::string, MapperFn> mapper_fns() const {
    return std::unordered_map<std::string, MapperFn> {
      {"rdm"     , &AppGraph<W>::map_random},
      {"rr"      , &AppGraph<W>::map_round_robin},
      {"ks"      , &AppGraph<W>::map_knapsack},
      {"sfcs"    , &AppGraph<W>::map_sfcs},
      {"pfcs"    , &AppGraph<W>::map_pfcs},
      {"pfcm"    , &AppGraph<W>::map_pfcm},
      {"gr"      , &AppGraph<W>::map_greedy},
      {"rcm"     , &AppGraph<W>::map_rcm_rcm},
      {"rcm-pfcm", &AppGraph<W>::map_rcm_pfcm},
#if KNOB_METIS
      {"rb"      , &AppGraph<W>::map_rb_rb},
      {"rbm"     , &AppGraph<W>::map_rb_rbm},
      {"rb-pfcm" , &AppGraph<W>::map_rb_pfcm},
#endif
    };
  }

  // helpers for mapping
  void init_network(NetworkGraph<W> &net_g, bool flag_bias = true, double cap_ratio = imba_ratio);
  void add_pair_traffic(NetworkGraph<W> &net_g, const AppNode<W> &cur,
                        const AppNode<W> &neighbor, bool flag_bidirect);
  void assign_node(NetworkGraph<W> &net_g, int app_id, int rank,
                   bool flag_add_traffic = false, bool flag_bidirect = true,
                   int app_connect_id = kNone);
  std::vector<RankWeight<W>> get_rank_weight_vector(const NetworkGraph<W> &net_g, int lev);

  static void remap_level(const std::vector<NodePtr> &node_vec,
                          std::vector<RankWeight<W>> &lev_rw_vec,
                          std::vector<RankWeight<W>> &rw_vec);
  static void map_knapsack_level(const std::vector<NodePtr> &node_vec,
                                 std::vector<RankWeight<W>> &rw_vec);
  template <class U>
  static void distribute_nodes(const std::vector<NodePtr> &node_vec,
                               std::vector<RankWeight<W>> &rw_vec);
  static void new_distribute_nodes(const std::vector<NodePtr> &node_vec,
                                   std::vector<RankWeightCap<W>> &rwc_vec);
  static void overwrite_subset_weights(const std::vector<RankWeight<W>> &subset_rw_vec,
                                       std::vector<RankWeight<W>> &rw_vec);

  // helpers for greedy algorithm

  // max heap of app edges
  struct EdgeHeap {
  public:
    struct El {
      int src, dst;
      size_t weight;
      friend std::ostream &operator<<(std::ostream &os, const El &e) {
        os << "Edge (" << e.src << ", " << e.dst << "): " << e.weight;
        return os;
      }
    };
    size_t size() const {
      return vec.size();
    }
    void push(El e) {
      vec.push_back(e);
      std::push_heap(vec.begin(), vec.end(), comp_f);
    }
    El pop() {
      std::pop_heap(vec.begin(), vec.end(), comp_f);
      El result = vec.back();
      vec.pop_back();
      return result;
    }

  private:
    std::vector<El> vec;
    static bool comp_f(const El &a, const El &b) { return a.weight < b.weight; }
  };

  int find_most_connected(std::unordered_set<int> S) const;
  int greedy_find_space(NetworkGraph<W> &net_g, int cur_app, int cur_net);
  int greedy_assign(NetworkGraph<W> &net_g, int cur_app, int prev_app,
                    int rank, std::unordered_set<int> &S, EdgeHeap &Q);

  using Edge = typename AppNode<W>::Edge;
  void print(std::ostream &os) const override {
    this->visit([&](const AppNode<W> &n) {
      os << n << std::endl;
      for (const auto &p : n.edges) {
        int dst = p.first;
        const Edge &e = p.second;
        os << "  Edge (" << n.id << ", " << dst << "): " << e.weight() << std::endl;
      }
    });
  }
};

/*********************
 * Application Graph *
 *********************/

template <class W>
void AppGraph<W>::reset_comp_comm() {
  this->visit([](AppNode<W> &n) {
    n.load = 0;
  });
  this->visit_edges([](int src, int dst, Edge &e) {
    e.to_msg_n    = 0;
    e.to_byte_n   = 0;
    e.from_msg_n  = 0;
    e.from_byte_n = 0;
  });
}

template <class W>
void AppGraph<W>::import(const std::unordered_map<int, double> &comps,
                         const std::unordered_map<std::pair<int, int>,
                                                  std::pair<int, size_t>> &comms) {
  for (const auto &p : comps) {
    auto node_id = p.first;
    auto secs = p.second;
    this->get_node(node_id).load += secs;
  }
  for (const auto &p : comms) {
    auto src_id = p.first.first;
    auto dst_id = p.first.second;
    auto msg_n  = p.second.first;
    auto byte_n = p.second.second;
    Edge &sd = this->get_edge_create(src_id, dst_id);
    Edge &ds = this->get_edge_create(dst_id, src_id);
    sd.to_msg_n    += msg_n;
    sd.to_byte_n   += byte_n;
    ds.from_msg_n  += msg_n;
    ds.from_byte_n += byte_n;
  }
}

template <class W>
void AppGraph<W>::print_stats() const {
  size_t total_msg_n = 0, total_byte_n = 0;
  this->visit_edges([&](int src, int dst, const Edge &e) {
    assert(e.to_msg_n  == this->get_edge(dst, src).from_msg_n);
    assert(e.to_byte_n == this->get_edge(dst, src).from_byte_n);
    total_msg_n += e.to_msg_n;
    total_byte_n += e.to_byte_n;
    if (flag_verbose_stats) {
      Say() << "Edge (" << src << " -> " << dst << "): ("
            << e.to_msg_n << ", " << e.to_byte_n << ")";
    }
  });
  if (!flag_quiet) {
    Say() << "Total messages in app graph: " << total_msg_n;
    Say() << "Total    bytes in app graph: " << total_byte_n;
  }
  if (flag_verbose_stats) { 
    this->visit_edges([&](int src, int dst, const Edge &e) {
      double weighted = e.to_byte_n * 100.0 / total_byte_n;
      Say() << "Normed Edge (" << src << ", " << dst << "): " << weighted;
    });
    this->visit([&](const AppNode<W> &n) {
      Say() << "Node " << n.id << ": " << n.load;
    });
  }
}

/***************************************
 * Application Graph Mapping Functions *
 ***************************************/

// initialize the network
template <class W>
void AppGraph<W>::init_network(NetworkGraph<W> &net_g, bool flag_bias, double cap_ratio) {
  // inspect application graph for total and maximum node weights
  int lev_n = max_node_level() + 1;
  std::vector<W> rank_cap(lev_n);
  for (int lev = 0; lev < lev_n; ++lev) {
    W sum_weight(0.0);
    W max_weight(std::numeric_limits<W>::min());
    auto node_vec = nodes_in_level(lev);
    for (const auto &pn : node_vec) {
      sum_weight += pn->weight();
      // NOTE: not std::max (can't overload std::max b/c of header include definition order)
      //       can overload by defining custom max(W a, W b) function
      max_weight = max(max_weight, pn->weight());
    }

    int lev_rank_n = !rank_masks_.empty() ? rank_masks_[lev].rank_n : net_g.rank_n();
    W avg_weight = sum_weight / static_cast<double>(lev_rank_n);
    rank_cap[lev] = max(cap_ratio * avg_weight, max_weight);
  }

  // initialize network capacities / congestion bias
  double bias = 0;
  if (flag_bias) {
    size_t max_edge_weight = 0;
    this->visit_edges([&](int s, int d, const Edge &e) {
      max_edge_weight = std::max(max_edge_weight, e.weight()); // std:: version fine here
    });
    bias = static_cast<double>(max_edge_weight) * this->node_n() * this->node_n();
  }

  if (!rank_masks_.empty()) {
    net_g.reset(0, bias);
    for (int lev = 0; lev < lev_n; ++lev) {
      net_g.inc_rank_caps(rank_masks_[lev], rank_cap[lev]);
    }
  } else {
    W sum_rank_cap(0);
    for (int lev = 0; lev < lev_n; ++lev) {
      sum_rank_cap += rank_cap[lev];
    }
    net_g.reset(sum_rank_cap, bias);
  }

  if (!flag_quiet) {
    Say() << "Initializing network graph load capacities:";
    Say() << "  Rank capacities: " << seq_to_str(rank_cap);
    Say() << "        Edge bias: " << bias;
    Say() << "  Imbalance ratio: " << cap_ratio;
  }
}

template <class W>
void AppGraph<W>::add_pair_traffic(NetworkGraph<W> &net_g, const AppNode<W> &cur,
                                   const AppNode<W> &nbr, bool flag_bidirect) {
  const Edge &e = cur.edges.at(nbr.id);
  int cur_net_id = net_g.node_at_rank(cur.rank).id;
  int nbr_net_id = net_g.node_at_rank(nbr.rank).id;

  // route traffic from here to neighbor
  if (flag_verbose_assign) {
    Say() << "  Adding traffic: ("
          << cur.id << ", " << cur.rank << ", " << cur_net_id
          << ") <--> ("
          << nbr.id << ", " << nbr.rank << ", " << nbr_net_id
          << "): (" << e.to_msg_n << ", " << e.to_byte_n << ")";
  }
  net_g.add_traffic(cur.rank, nbr.rank, e.to_msg_n, e.to_byte_n);

  if (flag_bidirect) {
    // route traffic from neighbor to here, may take a different route
    if (flag_verbose_assign) {
      Say() << "  Adding traffic: ("
            << nbr.id << ", " << nbr.rank << ", " << nbr_net_id
            << ") <--> ("
            << cur.id << ", " << cur.rank << ", " << cur_net_id
            << "): (" << e.from_msg_n << ", " << e.from_byte_n << ")";
    }
    net_g.add_traffic(nbr.rank, cur.rank, e.from_msg_n, e.from_byte_n);
  }
}

// general node assignment function
// flag_add_traffic: keep track of network traffic?
//                   set to true if we want to see network statistics at the end
//                   or if the mapping algorithm uses this information (e.g. greedy)
// flag_bidirect: when adding traffic on edges of the network, do so in both directions?
//                if this is called after all nodes have been assigned, set to false to avoid double counting
//                if this is called as each node is assigned incrementally, set to true
//                traffic in both directions are routed independently
template <class W>
void AppGraph<W>::assign_node(NetworkGraph<W> &net_g, int app_id, int rank,
                              bool flag_add_traffic, bool flag_bidirect,
                              int app_connect_id) {
  // get node objects
  AppNode<W> &app_node = this->get_node(app_id);
  NetworkNode<W> &net_node = net_g.node_at_rank(rank);
  // assign app node to network node
  if (flag_verbose_assign) Say() << "Assigning " << app_node << "\n  to rank " << rank << "\n    at " << net_node;
  app_node.rank = rank;
  net_g.add_load(rank, app_node.weight());
  // add weights of app edges to network edges traversed
  if (flag_add_traffic) {
    if (flag_fast_add_traffic) {
      if (app_connect_id != kNone) {
        // connect only to the one node for speed, but less accurate
        add_pair_traffic(net_g, app_node, this->get_node(app_connect_id), flag_bidirect);
      } else {
        // do nothing
      }
    } else {
      // connect all
      for (const auto &p : app_node.edges) {
        const AppNode<W> &neighbor = this->get_node(p.first);
        // NOTE: traffic is not added if neighbor hasn't been assigned yet
        if (neighbor.rank != kNone) {
          add_pair_traffic(net_g, app_node, neighbor, flag_bidirect);
        }
      }
    }
  }
}

template <class W>
void AppGraph<W>::write_mapping(std::string outfile) const {
  std::ofstream os(outfile);
  if(!os) { Say() << "Could not open file: " << outfile; abort(); }
  os << "<boxes>\n";
  for (unsigned i = 0; i < this->node_n(); ++i) {
    int rank = this->get_node(i).rank;
    assert(rank != kNone); // check node has assignment
    os << "<box id=\"B" << i << "\" ";
    os << "loc=\"" << rank << "\" ";
    os << "cells=\"" << this->get_node(i).cell_n << "\" ";
    os << "/>\n";
  }
  os << "</boxes>\n";
}

template <class T, class U>
std::shared_ptr<T> polymorphic_downcast(std::shared_ptr<U> p) {
#ifdef NDEBUG
  return std::static_pointer_cast<T>(p);
#else
  auto result = std::dynamic_pointer_cast<T>(p);
  assert(result != nullptr);
  return result;
#endif
}

template <class W>
void AppGraph<W>::evaluate_mapping(std::shared_ptr<NetworkGraph_> _net_g,
                                   const std::unordered_map<int, int> &mapping) {
  auto net_g = polymorphic_downcast<NetworkGraph<W>>(_net_g);
  clear_mapping();
  init_network(*net_g);
  this->visit([&](const AppNode<W> &n) {
    assign_node(*net_g, n.id, mapping.at(n.id), true, true);
  });
  net_g->print_stats();
}

template <class W>
void AppGraph<W>::init_rank_masks(std::shared_ptr<NetworkGraph_> _net_g,
                                  const std::vector<int> &map_rank_n) {
  const auto &net_g = *polymorphic_downcast<NetworkGraph<W>>(_net_g);
  int lev_n = map_rank_n.size();
  rank_masks_.clear();
  rank_masks_.reserve(lev_n);
  for (int i = 0; i < lev_n; ++i) {
    rank_masks_.emplace_back(net_g.get_hood(), map_rank_n[i]);
  }
}

// dispatch
template <class W>
void AppGraph<W>::run_mapper(std::string mapper_str,
                             std::shared_ptr<NetworkGraph_> _net_g) {
  auto net_g = polymorphic_downcast<NetworkGraph<W>>(_net_g);
  auto fns = mapper_fns();
  if (fns.count(mapper_str)) {
    auto mf = fns.at(mapper_str);
    (this->*mf)(*net_g);
  } else {
    Say() << "ERROR: unknown mapper: " << mapper_str;
    abort();
  }
}

template <class W>
std::vector<RankWeight<W>>
AppGraph<W>::get_rank_weight_vector(const NetworkGraph<W> &net_g, int lev) {
  return !rank_masks_.empty() ? net_g.rank_weights(rank_masks_[lev])
                              : net_g.rank_weights();
}

/**********
 * Random *
 **********/

template <class W>
void AppGraph<W>::map_random(NetworkGraph<W> &net_g) {
  // initialize app and network graphs
  clear_mapping();
  init_network(net_g);

  for (int lev = 0; lev <= max_node_level(); ++lev) {
    std::vector<NodePtr> node_vec = nodes_in_level(lev);
    auto rw_vec = get_rank_weight_vector(net_g, lev);
    assert(rw_vec.size());

    // random assignment
    std::uniform_int_distribution<int> udist(0,rw_vec.size()-1);
    for (auto &pn : node_vec) {
      const RankWeight<W> &rw = rw_vec[random::roll(udist)];
      assign_node(net_g, pn->id, rw.rank, flag_force_traffic);
    }
  }
}

/***************
 * Round-robin *
 ***************/

template <class W>
void AppGraph<W>::map_round_robin(NetworkGraph<W> &net_g) {
  // initialize app and network graphs
  clear_mapping();
  init_network(net_g);

  for (int lev = 0; lev <= max_node_level(); ++lev) {
    // sort nodes by decreasing weight
    auto node_vec = sorted_nodes_in_level(node_scalar_weight_gt, lev);
    auto rw_vec = get_rank_weight_vector(net_g, lev);
    assert(rw_vec.size());
    // sort ranks by increasing weight
    sort_by_scalar_weight(rw_vec);

    // round robin assignment
    for (unsigned i = 0; i < node_vec.size(); ++i) {
      AppNode<W> &an = *node_vec[i];
      RankWeight<W> &rw = rw_vec[i % rw_vec.size()];
      assign_node(net_g, an.id, rw.rank, flag_force_traffic);
    }
  }
}

/*******************
 * Greedy Knapsack *
 *******************/

// remap nodes that were mapped to lev_rw_vec to rw_vec_orig
// sort both vecs and put most weighted ranks onto least weighted ranks
// also reassign ranks saved in node_vec nodes
template <class W>
void AppGraph<W>::remap_level(const std::vector<NodePtr> &node_vec,
                              std::vector<RankWeight<W>> &lev_rw_vec,
                              std::vector<RankWeight<W>> &rw_vec_orig)
{
  bool flag_subset_remap = lev_rw_vec.size() != rw_vec_orig.size();
  unsigned rank_n = lev_rw_vec.size();

  std::unordered_map<int, int> rank_to_idx; // map from ID to index in rw_vec
  std::vector<RankWeight<W>> rw_vec_subset;
  std::vector<RankWeight<W>> *rw_vec = &rw_vec_orig;

  // handle case where lev_rw_vec has a proper subset of ranks in rw_vec_orig
  if (flag_subset_remap) {
    // rank_to_idx map only contains ranks in lev_rw_vec
    for (const auto &rw : lev_rw_vec) {
      rank_to_idx[rw.rank] = -1;
    }
    rw_vec_subset.reserve(rank_n);
    for (unsigned i = 0; i < rw_vec_orig.size(); ++i) {
      const auto &rw = rw_vec_orig[i];
      if (rank_to_idx.count(rw.rank)) {
        rw_vec_subset.push_back(rw);
        rank_to_idx[rw.rank] = i;
      }
    }
    rw_vec = &rw_vec_subset;
  }
  assert(rw_vec->size() == lev_rw_vec.size());

  sort_by_scalar_weight(lev_rw_vec);
  sort_by_scalar_weight(*rw_vec);

  std::unordered_map<int, int> mapping;
  for (unsigned i = 0; i < rank_n; ++i) {
    int target_rank = (*rw_vec)[rank_n-i-1].rank;
    mapping[lev_rw_vec[i].rank] = target_rank;
    int orig_idx = flag_subset_remap ? rank_to_idx[target_rank] : rank_n-i-1;
    rw_vec_orig[orig_idx].weight += lev_rw_vec[i].weight;
  }
  for (const NodePtr &pn : node_vec) {
    pn->rank = mapping[pn->rank];
  }
}

// NOTE: assumes node_vec is in the order we want to assign them, and
//       rw_vec is already a min heap
template <class W>
void AppGraph<W>::map_knapsack_level(const std::vector<NodePtr> &node_vec,
                                     std::vector<RankWeight<W>> &rw_vec) {
  // for each node in node_vec
  for (const NodePtr &pn : node_vec) {
    // map it to the rank with the lowest weight
    std::pop_heap(rw_vec.begin(), rw_vec.end(), rw_weight_gt<W, double>);
    pn->rank = rw_vec.back().rank;
    rw_vec.back().weight += pn->weight();
    std::push_heap(rw_vec.begin(), rw_vec.end(), rw_weight_gt<W, double>);
  }
}

template <class W>
void AppGraph<W>::map_knapsack(NetworkGraph<W> &net_g)
{
  // initialize app and network graphs
  clear_mapping();
  init_network(net_g);

  auto rw_vec = net_g.rank_weights();
  // order and map each level separately
  for (int lev = 0; lev <= max_node_level(); ++lev) {
    // sort nodes by decreasing weight
    auto node_vec = sorted_nodes_in_level(node_scalar_weight_gt, lev);
    auto lev_rw_vec = get_rank_weight_vector(net_g, lev);
    reset_weights(lev_rw_vec, W(0.0));
    map_knapsack_level(node_vec, lev_rw_vec);
    remap_level(node_vec, lev_rw_vec, rw_vec);
  }
  // make assignments
  this->visit([&](const AppNode<W> &n) {
    assign_node(net_g, n.id, n.rank, flag_force_traffic, false);
  });
}

/*******
 * SFC *
 *******/

// This version of distribute_nodes is based on BoxLib distribute()
// Instead of a per-node weight capacity limit, it uses an increasing cumulative limit
//   as it scans from the beginning of the list, spacing large nodes out.
// Single pass only: assign remaining nodes to last list element (flag_last).
//   This doesn't work well for multi-constraint weights b/c nodes get dumped
//   on the last rank.
//
// W is the class of the weights in the App and Network Graphs
// U is the class to static_cast the weights to for purposes of distribution
//   this is useful if e.g. we want to use a scalarized version of the weights
template <class W>
template <class U>
void AppGraph<W>::distribute_nodes(const std::vector<NodePtr> &node_vec,
                                   std::vector<RankWeight<W>> &rw_vec)
{
  if (!flag_quiet && flag_verbose_distribute) Say() << "Distribute begin:";
  // compute average weight
  U avg_weight; {
    U sum_weight = 0;
    for (const NodePtr &pn : node_vec) { sum_weight += static_cast<U>(pn->weight()); }
    for (const RankWeight<W> &rw : rw_vec) { sum_weight += static_cast<U>(rw.weight); }
    avg_weight = sum_weight / rw_vec.size();
    if (!flag_quiet && flag_verbose_distribute) {
      Say() << "  sum_weight = " << sum_weight << ", avg_weight = " << avg_weight;
    }
  }

  unsigned node_idx = 0; // app node idx to be mapped
  U weight_so_far = 0;
  for (unsigned rw_idx = 0; rw_idx < rw_vec.size(); ++rw_idx) {
    RankWeight<W> &rw = rw_vec[rw_idx];
    U rw_weight = static_cast<U>(rw.weight); // handy alias
    weight_so_far += rw_weight;
    if (!flag_quiet && flag_verbose_distribute) {
      Say() << "  bin " << rw_idx << ": weight = " << rw_weight;
      Say() << "  total weight_so_far = " << weight_so_far;
    }
    bool flag_last = (rw_idx == rw_vec.size()-1);
    while (node_idx < node_vec.size()) {
      AppNode<W> &n = *node_vec[node_idx];
      U n_weight = static_cast<U>(n.weight()); // handy alias
      if (!flag_quiet && flag_verbose_distribute) {
        Say() << "  node " << node_idx << ": weight = " << n_weight;
        Say() << "    bin check: " << rw_weight << " ?<= " << avg_weight << " = " << (rw_weight <= avg_weight);
        Say() << "  cumul. check: " << n_weight << " ?<= " << avg_weight*(rw_idx+1) - weight_so_far
              << " = " << (n_weight <= avg_weight*(rw_idx+1) - weight_so_far);
      }
      if (flag_last ||
          (rw_weight <= avg_weight &&
           weight_so_far+n_weight <= avg_weight*(rw_idx+1))) {
        // assign node
        if (!flag_quiet && flag_verbose_distribute) {
          Say() << "  assigning node " << node_idx << " to bin " << rw_idx;
        }
        n.rank = rw.rank;
        rw.weight += n.weight();
        rw_weight = static_cast<U>(rw.weight);
        weight_so_far += n_weight;
        ++node_idx;
      } else {
        break; // go to next
      }
    }
  }
}

// New version of distribute that uses a capacity based on the mean and max weights.
//   Uses multiple passes until all nodes have been assigned.
//   This capacity is raised until a valid assignment is found.
template <class W>
void AppGraph<W>::new_distribute_nodes(const std::vector<NodePtr> &node_vec,
                                       std::vector<RankWeightCap<W>> &rwc_vec)
{
  // compute capacity ratio to multiply average by
  double cap_ratio = 1.01;
  const double cap_multiplier = std::sqrt(2.0);

  bool flag_done = false;
  while (!flag_done) {
    // reset rank assigned weights
    reset_weights(rwc_vec, W(0.0));

    // assign nodes to bins in rwc_vec
    unsigned rw_idx = 0;
    int direction = +1; // start searching in positive direction
    flag_done = true; // unless set to false inside loop
    for (auto np : node_vec) {
      AppNode<W> &n = *np;
      // find a bin where this node fits
      unsigned counter = 0;
      while (!(rwc_vec[rw_idx].weight + n.weight() <= rwc_vec[rw_idx].cap * cap_ratio)) {
        ++counter;
        if (counter >= 2*rwc_vec.size()) { flag_done = false; break; }
        // swap directions if needed
        if (rw_idx >= rwc_vec.size()-1) {
          direction = -1;
        } else if (rw_idx <= 0) {
          direction = +1;
        }
        rw_idx += direction;
      }
      if (!flag_done) {
        cap_ratio = 1+(cap_multiplier*(cap_ratio-1));
        if (!flag_quiet && flag_verbose_distribute) {
          Say() << "Distribute: retrying with cap_ratio = " << cap_ratio;
        }
        break;
      } else {
        // assign node to bin
        n.rank = rwc_vec[rw_idx].rank;
        auto old_space = rwc_vec[rw_idx].cap * cap_ratio - rwc_vec[rw_idx].weight;
        rwc_vec[rw_idx].weight += n.weight();
        auto new_space = old_space - n.weight();
        if (!flag_quiet && flag_verbose_distribute) {
          Say() << "Assigning node ("
                << n.id   << ", " << n.weight() << ") to ("
                << rwc_vec[rw_idx].rank << ", " << rw_idx << ", " << old_space  << ") -> ("
                << rwc_vec[rw_idx].rank << ", " << rw_idx << ", " << new_space  << ")";
        }
      }
    }
  }
}

// overwrite the weights in rw_vec with each rank's weight in subset_rw_vec
template <class W>
void AppGraph<W>::overwrite_subset_weights(const std::vector<RankWeight<W>> &subset_rw_vec,
                                           std::vector<RankWeight<W>> &rw_vec) {
  std::unordered_map<int, int> rank_to_idx; // map from ID to index in rw_vec
  for (unsigned i = 0; i < rw_vec.size(); ++i) {
    rank_to_idx[rw_vec[i].rank] = i;
  }
  for (const auto &rw : subset_rw_vec) {
    rw_vec[rank_to_idx[rw.rank]].weight = rw.weight;
  }
}

template <class W>
void AppGraph<W>::map_sfc(NetworkGraph<W> &net_g, bool flag_per_lev, bool flag_remap_levs, std::string tag)
{
  // initialize app and network graphs
  clear_mapping();
  init_network(net_g);

  auto rw_vec = net_g.rank_weights();
  if (flag_per_lev) {
    // order and map each level separately
    for (int lev = 0; lev <= max_node_level(); ++lev) {
      // sort nodes by Z-Morton ordering
      auto node_vec = sorted_nodes_in_level(node_zmorton_lt, lev);
      if (flag_remap_levs) {
        // SFCS: // map onto temporary rw_vec, then merge for better load balance
        // remap destroys spatial locality between bins since they are sorted by weight
        auto lev_rw_vec = get_rank_weight_vector(net_g, lev);
        reset_weights(lev_rw_vec, W(0.0));
        distribute_nodes<double>(node_vec, lev_rw_vec);
        remap_level(node_vec, lev_rw_vec, rw_vec);
      } else {
        // PFCS: // map each level onto bins incrementally without remapping
        // preserves spatial locality between bins
        // since rw_vec is sorted in a spatial-aware way (if possible)
        auto lev_rw_vec = get_rank_weight_vector(net_g, lev);
        distribute_nodes<double>(node_vec, lev_rw_vec);
        overwrite_subset_weights(lev_rw_vec, rw_vec);
      }
    }
  } else {
    // PFCM: // order and map all levels simultaneously
    // may result in worse level load balance since all levels are together
    // also preserves spatial locality between bins
    if (!rank_masks_.empty()) {
      Say() << "WARNING: mota::pfcm mapper ignores rank masks (maps to all ranks)!";
    }
    auto node_vec = this->sorted_nodes(node_zmorton_lt);
    distribute_nodes<double>(node_vec, rw_vec);
  }
  // make assignments
  this->visit([&](const AppNode<W> &n) {
    assign_node(net_g, n.id, n.rank, flag_force_traffic, false);
  });
}

/**********************
 * Greedy Topological *
 **********************/

template <class W>
int AppGraph<W>::find_most_connected(std::unordered_set<int> S) const {
  int result = kNone;
  size_t max_weight = 0;
  for (int id : S) {
    size_t w = this->get_node(id).sum_edge_weight();
    if (w > max_weight) {
      max_weight = w;
      result = id;
    }
  }
  return result;
}

// find a nearby network rank with space for the current application node
template <class W>
int AppGraph<W>::greedy_find_space(NetworkGraph<W> &net_g, int cur_app_id, int cur_net_id) {
  AppNode<W> &app_node = this->get_node(cur_app_id);
  return net_g.find_space(app_node.weight(), cur_net_id);
}

// assign an app node to a network rank
template <class W>
int AppGraph<W>::greedy_assign(NetworkGraph<W> &net_g, int cur_app_id, int prev_app_id,
                               int rank, std::unordered_set<int> &S, EdgeHeap &Q) {

  AppNode<W> &app_node = this->get_node(cur_app_id);
  int cur_net_id = net_g.node_at_rank(rank).id;

  // assign the app node to the network node
  assign_node(net_g, cur_app_id, rank, true, true, prev_app_id);
  S.erase(cur_app_id);

  // add cur_app_id node's edges to Q
  for (const auto &p : app_node.edges) {
    const int src = cur_app_id, dst = p.first;
    if (S.count(dst) > 0) { // only add edge if dst is in S
      const Edge &e = p.second;
      typename EdgeHeap::El el {src, dst, e.weight()};
      Q.push(el);
    }
  }

  return cur_net_id; // current network node
}

// Greedy algorithm based on Hoefler SC'12 paper
//   cur_app = m, prev_app = u, cur_net = s
template <class W>
void AppGraph<W>::map_greedy(NetworkGraph<W> &net_g) {
  // compute capacity ratio
  double cap_ratio = 1.01;
  const double cap_multiplier = std::sqrt(2.0);

  bool flag_fail;
  do {
    // initialize app and network graphs
    clear_mapping();
    init_network(net_g, true, cap_ratio);

    // S is the set of unassigned app node ids
    std::unordered_set<int> S;
    this->visit([&](const AppNode<W> &n) { S.insert(n.id); });

    // Q is a heap of edges from assigned to unassigned app nodes
    EdgeHeap Q;
    int cur_net_id = kNone, cur_app_id = kNone, prev_app_id = kNone, rank = kNone;

    // while there are unassigned app nodes
    flag_fail = false;
    while (!flag_fail && S.size() > 0) {
      cur_app_id = find_most_connected(S);
      if (flag_verbose_greedy) Say() << "\nFound most connected node (scan): " << this->get_node(cur_app_id);
      rank = greedy_find_space(net_g, cur_app_id, cur_net_id);
      if (!(flag_fail = (rank == kNone))) {
        cur_net_id = greedy_assign(net_g, cur_app_id, prev_app_id, rank, S, Q);
      }
      while (!flag_fail && S.size() > 0 && Q.size() > 0) {
        typename EdgeHeap::El el = Q.pop();
        if (S.count(el.dst) > 0) { // only traverse edge if node is in S
          if (flag_verbose_greedy) Say() << "\nFound most connected node (heap): " << el;
          cur_app_id = el.dst;
          prev_app_id = el.src;
          rank = greedy_find_space(net_g, cur_app_id, cur_net_id);
          if (!(flag_fail = (rank == kNone))) {
            cur_net_id = greedy_assign(net_g, cur_app_id, prev_app_id, rank, S, Q);
          }
        }
      }
    }
    if (flag_fail) {
      cap_ratio = 1+(cap_multiplier*(cap_ratio-1));
      if (flag_verbose_greedy) {
        Say() << "Greedy: retrying with cap_ratio = " << cap_ratio;
      }
    }
  } while (flag_fail);
}


/*******
 * RCM *
 *******/

namespace _ {
  template <class NodeW, class EdgeW>
  std::vector<int> rcm_perm(const CSRMat<NodeW, EdgeW> &mat) {
    const int n = mat.size(), flags = 0;
    std::vector<int> perm(n,0), deg(n,0);
    std::vector<signed char> mask(n,0);
    genrcmi(n, flags, mat.xadj().data(), mat.adj().data(),
            perm.data(), mask.data(), deg.data());
    return mat.idxs_to_node_ids(perm);
  }
}

template <class W>
void AppGraph<W>::map_rcm(NetworkGraph<W> &net_g, bool flag_pfcm_override) {
  // initialize app and network graphs
  clear_mapping();
  init_network(net_g);
  
  // order app nodes
  std::vector<NodePtr> node_vec;
  if (flag_pfcm_override) {
    // just use the PFCM ordering for the app nodes
    node_vec = this->sorted_nodes(node_zmorton_lt);
  } else {
    // get app permutation
    auto app_perm = _::rcm_perm(this->get_csr());
    node_vec = this->make_node_vec(app_perm);
  }

  // get net permutation
  auto net_perm = _::rcm_perm(net_g.get_csr_no_switch());
  auto rwc_vec = net_g.rank_weight_caps(net_perm);

  // distribute app nodes to network nodes
  new_distribute_nodes(node_vec, rwc_vec);

  // make assignments
  this->visit([&](const AppNode<W> &n) {
    assign_node(net_g, n.id, n.rank, flag_force_traffic, false);
  });
}


/***********************
 * Recursive Bisection *
 ***********************/

#if KNOB_METIS
namespace _ {
  template <class T>
  std::vector<int64_t> to_int64(const std::vector<T> &vec) {
    std::vector<int64_t> result;
    result.reserve(vec.size());
    for (auto val : vec) {
      result.push_back(static_cast<int64_t>(val));
    }
    return result;
  }

  template <class T>
  std::vector<int64_t> scale_to_int64(const std::vector<T> &vec) {
    double max_aval = 0;
    for (auto val : vec) {
      auto dval = static_cast<double>(val);
      if (fabs(dval) > max_aval) max_aval = fabs(dval);
    }
    double target = std::sqrt(double(std::numeric_limits<int64_t>::max()));
    double scale_factor = target / max_aval;
    std::vector<int64_t> result;
    result.reserve(vec.size());
    for (auto val : vec) {
      auto dval = static_cast<double >(val);
      auto ival = static_cast<int64_t>(scale_factor * dval);
      result.push_back(ival);
    }
    return result;
  }

  template <class T>
  void summarize_vec(std::vector<T> vec) {
    Say() << "size=" << vec.size();
    auto min_v = std::numeric_limits<T>::max();
    auto max_v = std::numeric_limits<T>::min();
    for (const auto &v : vec) {
      min_v = std::min(v, min_v);
      max_v = std::max(v, max_v);
    }
    Say() << ", min=" << min_v;
    Say() << ", max=" << max_v;
    Say() << ", (";
    for (size_t i = 0; i < std::min(4ul, vec.size()); ++i) {
      Say() << vec[i] << ", ";
    }
    Say() << "..., ";
    for (size_t i = std::max(0ul, vec.size()-4); i < vec.size(); ++i) {
      Say() << vec[i] << ", ";
    }
    Say() << "\b\b)";
  }

  template <class NodeW, class EdgeW>
  std::vector<double>
  get_constraints(const CSRMat<NodeW, EdgeW> &mat,
                  const std::vector<std::function<double(NodeW)>> &cons) {
    std::vector<double> result;
    result.reserve(mat.size() * cons.size());
    for (auto id : mat.node_ids()) {
      for (auto lam : cons) {
        result.push_back(lam(mat.node_weight(id)));
      }
    }
    return result;
  }

  template <class NodeW, class EdgeW>
  std::vector<int64_t>
  metis_part(const CSRMat<NodeW,EdgeW> &mat,
             int64_t nparts,
             std::vector<std::function<double(NodeW)>> cons) {
    static_assert(sizeof(idx_t) == sizeof(int64_t), "Requires 64-bit METIS idx_t");
    static_assert(sizeof(real_t) == sizeof(double), "Requires 64-bit METIS real_t");
    int64_t nvtxs = mat.size();
    // fast, easy case
    if (nparts == 1 || nvtxs == 1) { return std::vector<int64_t>(nvtxs,0); }
    // set up metis parameters
    int64_t ncon = cons.size(), objval;
    std::vector<int64_t> options(METIS_NOPTIONS), part(nvtxs);
    std::vector<double> tpwgts(nparts*ncon, 1.0 / nparts),
                        ubvec(ncon, imba_ratio);
    // set up weights vectors
    std::vector<int64_t> constraints = scale_to_int64(get_constraints(mat, cons)),
                         edge_wgts   = scale_to_int64(mat.edge_wgts());
    // check parameters
    if (flag_verbose_metis) {
      Say() << "nparts: " << nparts;
      Say() << "nvtx: " << nvtxs;
      Say() << "xadj: ";
      summarize_vec(mat.xadj());
      Say() << " adj: ";
      summarize_vec(mat.adj());
      Say() << "ncon: " << ncon;
    }
    // do partitioning
    METIS_SetDefaultOptions(options.data());
    METIS_PartGraphRecursive( &nvtxs, &ncon,
      to_int64(mat.xadj()).data(),
      to_int64(mat.adj()).data(),
      constraints.data(), nullptr,
      edge_wgts.data(),
      &nparts, tpwgts.data(), ubvec.data(),
      options.data(), &objval, part.data()
    );
    return part;
  }

  inline std::vector<int> part_to_perm(std::vector<int64_t> part) {
    // build (partition, index) pairs out of part vector
    std::vector<std::array<int, 2>> pp; {
      for (unsigned idx = 0; idx < part.size(); ++idx) {
        pp.push_back({static_cast<int>(part[idx]),
                      static_cast<int>(idx)});
    }}
    // sort indexes by which partition they were assigned to
    sort(pp.begin(), pp.end());
    // build permutation vector
    std::vector<int> result; {
      for (const auto &p : pp) result.push_back(p[1]);
    }
    return result;
  }
}

template <class W>
struct MetisCons;

// specialize to use different set of constraints for your node weight
template <class W>
struct MetisCons {
  // default single constraint: static cast weight to double
  static std::vector<std::function<double(W)>> cons(bool flag_multicons) {
    return std::vector<std::function<double(W)>> {
      [] (W x) { return static_cast<double>(x); }
    };
  }
};

// NOTE: this sometimes generates a (cannot bisect a graph with 0 vertices)
//       warning from METIS.  Seems harmless...
template <class W>
void AppGraph<W>::map_rb(NetworkGraph<W> &net_g, bool flag_multicons, bool flag_pfcm_override) {
  // initialize app and network graphs
  clear_mapping();
  init_network(net_g);

  // get sizes
  auto app_n = this->node_n();
  auto rank_n = net_g.rank_n();

  // order app nodes
  std::vector<NodePtr> node_vec;
  if (flag_pfcm_override) {
    // just use the PFCM ordering for the app nodes
    node_vec = this->sorted_nodes(node_zmorton_lt);
  } else {
    // get app partition ordering
    //std::cout << "Partitioning app graph" << std::endl;
    auto app_csr = this->get_csr();
    auto app_part = _::metis_part(app_csr, std::min(rank_n, app_n), MetisCons<W>::cons(flag_multicons));
    auto app_perm = _::part_to_perm(app_part);
    app_perm = app_csr.idxs_to_node_ids(app_perm);
    node_vec = this->make_node_vec(app_perm);
  }

  // get network partition ordering
  //std::cout << "Partitioning net graph" << std::endl;
  auto net_csr = net_g.get_csr_no_switch();
  auto net_part = _::metis_part(net_csr, rank_n, MetisCons<W>::cons(flag_multicons));
  auto net_perm = _::part_to_perm(net_part);
  net_perm = net_csr.idxs_to_node_ids(net_perm);
  auto rwc_vec = net_g.rank_weight_caps(net_perm);

  // distribute app nodes to network nodes
  new_distribute_nodes(node_vec, rwc_vec);

  // make assignments
  this->visit([&](const AppNode<W> &n) {
    assign_node(net_g, n.id, n.rank, flag_force_traffic, false);
  });
}
#endif // KNOB_METIS

// create an application graph with weights for the right number of AMR levels
// TODO: there should be a better way to do this
inline std::shared_ptr<AppGraph_> app_graph_create(int amr_level_n) {
  std::shared_ptr<AppGraph_> result {nullptr};
  switch (amr_level_n) {
    case 1:
      result = std::shared_ptr<AppGraph_>(new AppGraph<AMRWeight<1>>());
      break;
    case 2:
      result = std::shared_ptr<AppGraph_>(new AppGraph<AMRWeight<2>>());
      break;
    case 3:
      result = std::shared_ptr<AppGraph_>(new AppGraph<AMRWeight<3>>());
      break;
    case 4:
      result = std::shared_ptr<AppGraph_>(new AppGraph<AMRWeight<4>>());
      break;
    case 5:
      result = std::shared_ptr<AppGraph_>(new AppGraph<AMRWeight<5>>());
      break;
    default:
      Say() << "ERROR: Number of AMR levels (" << amr_level_n << ") is over maximum supported (11)!";
      std::abort();
      break;
  }
  return result;
}

}

#endif // APPGRAPH_HXX
