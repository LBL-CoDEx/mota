#ifndef NETGRAPH_HXX
#define NETGRAPH_HXX

#include <memory>
#include <fstream>
#include <sstream>
#include <tuple>
#include "graph.hxx"
#include "rankweight.hxx"
#include "rankmask.hxx"
#include "flags.hxx"

namespace mota {

// forward declaration
template <class W> class AppGraph;

/*************************
 * Generic Network Graph *
 *************************/

// W is the type of capacity and utilization (weight) assigned to a rank
// W supports: construct/assign to/from double, compare, add, subtract,
//             multiply, divide by double, operator<<(ostream)
template <class W>
struct NetworkNode {
  using Weight = W;
  using Ptr = std::shared_ptr<NetworkNode>;

  struct Edge {
    using Weight = double;
    size_t byte_n = 0;
    int link_n = 1.0;
    bool node_link = false; // true iff links node & switch (not switch-to-switch)
    double congest_bias = 0.0; // congestion bias for mappers and routing
    Edge(int l, bool nl, double bias=0.0) : link_n{l}, node_link{nl}, congest_bias{bias} {}
    double weight() const { return (double) link_n; }
    double congest() const { return (double) byte_n / link_n; }
    double congest_biased() const { return congest() + congest_bias; }
    size_t edge_byte_n() const { return byte_n; }
    size_t link_byte_n() const { return byte_n / link_n; }
  };

  // each network node can have multiple ranks assigned to it
  struct Rank {
    int rank;          // job rank
    NetworkNode &node; // node it resides on
    W cap = 1.0;       // capacity
    W util = 0.0;      // utilization
    size_t comp_byte_n = 0;
    size_t on_rank_comm_msg_n = 0;
    size_t on_rank_comm_byte_n = 0;
    size_t on_node_comm_msg_n = 0;
    size_t on_node_comm_byte_n = 0;
    size_t off_node_comm_msg_n = 0;
    size_t off_node_comm_byte_n = 0;

    Rank(int _rank, NetworkNode &_node) : rank{_rank}, node{_node} {}
    W space() const { return cap - util; } // free space

    friend std::ostream &operator<<(std::ostream &os, const Rank &r) {
      os << "  Rank (rank, cap, util) = ("
         << r.rank << ", " << r.cap << ", " << r.util << ")";
      return os;
    }
  };
  using RankPtr = std::shared_ptr<Rank>;

  int id = kNone;
  std::vector<RankPtr> ranks{};          // ranks on this node
  std::unordered_map<int, Edge> edges{}; // keyed on neighbor node ID
  size_t router_msg_n;
  size_t router_byte_n;

  explicit NetworkNode(int _id) : id{_id} {}

  size_t rank_n() const { return ranks.size(); }

  // return sum of the space over ranks on the node
  W weight() const {
    W sum_space {0};
    for (const auto &r : ranks) { sum_space += r->space(); }
    return sum_space;
  }
  double scalar_weight() const { return static_cast<double>(weight()); }

  std::string ranks_to_str() const {
    std::stringstream oss;
    for (const auto &r : ranks) oss << *r << '\n';
    return oss.str();
  }
  friend std::ostream &operator<<(std::ostream &os, const NetworkNode &n) {
    os << "NetworkNode " << n.id << " ranks:\n" << n.ranks_to_str();
    return os;
  }
};

// collected stats
struct NetworkStats {
  bool flag_switch_hops_only = false;
  std::unordered_set<int> active_nodes;
  size_t total_msg_n = 0,  // incl. on-rank/on-node "messages"
         total_byte_n = 0, // incl. on-rank/on-node "messages"
         on_rank_msg_n = 0,
         on_rank_byte_n = 0,
         on_node_msg_n = 0,
         on_node_byte_n = 0,
         off_node_msg_n = 0,
         off_node_byte_n = 0,
         off_node_msg_hop_n = 0,
         off_node_byte_hop_n = 0,
         active_edge_n = 0,
         active_link_n = 0,
         max_edge_byte_n = 0,
         max_link_byte_n = 0;
  void reset() {
    active_nodes.clear();
    total_msg_n = total_byte_n = 0;
    on_rank_msg_n = on_rank_byte_n = 0;
    on_node_msg_n = on_node_byte_n = 0;
    off_node_msg_n = off_node_byte_n = 0;
    off_node_msg_hop_n = off_node_byte_hop_n = 0;
    active_edge_n = active_link_n = 0;
    max_edge_byte_n = max_link_byte_n = 0;
  }
  double avg_off_node_msg_hop_n() const { return off_node_msg_n ? (double) off_node_msg_hop_n / off_node_msg_n : 0; }
  double avg_off_node_byte_hop_n() const { return off_node_byte_n ? (double) off_node_byte_hop_n / off_node_byte_n : 0; }
  double avg_active_edge_byte_n() const { return active_edge_n ? (double) off_node_byte_hop_n / active_edge_n : 0; }
  double avg_active_link_byte_n() const { return active_link_n ? (double) off_node_byte_hop_n / active_link_n : 0; }

  static void write_stats_header(std::string filename);
};


// abstract class to a network graph
class NetworkGraph_
{
public:
  NetworkGraph_() = default;
  virtual ~NetworkGraph_() = default;

  virtual const std::string label() const = 0;
  virtual size_t rank_n() const = 0;
  virtual size_t switch_n() const = 0;
  virtual size_t compute_n() const = 0;
  virtual size_t link_n() const = 0;

  // stats functions
  virtual void append_stats(std::string filename, std::string mapper,
                            size_t map_ms) const = 0;
  virtual void print_stats() const = 0;
  virtual NetworkStats get_stats() = 0;

  // cost estimation
  virtual double estimate_cost() const = 0;

  // find a neighborhood of nodes around a rank
  virtual std::pair<std::vector<int>, double>
  neighborhood(int src_rank_id, std::vector<std::pair<int, int>> rank_n_weights) const = 0;
  virtual void set_hood(const std::vector<int> *hood_ptr) = 0;
  virtual const std::vector<int> & get_hood() const = 0; 
};


// concrete network graph using W for weight type
template <class W>
class NetworkGraph : public Graph<NetworkNode<W>>, public NetworkGraph_ {
  friend AppGraph<W>;
public:
  using Node    = NetworkNode<W>;
  using NodeW   = typename NetworkNode<W>::Weight;
  using NodePtr = typename NetworkNode<W>::Ptr;
  using Edge    = typename NetworkNode<W>::Edge;
  using EdgeW   = typename NetworkNode<W>::Edge::Weight;
  using Rank    = typename NetworkNode<W>::Rank;
  using RankPtr = typename NetworkNode<W>::RankPtr;
  using Path    = std::vector<int>;

  NetworkGraph() = default;
  NetworkGraph(const NetworkGraph &) = default;
  NetworkGraph(NetworkGraph &&) = default;
  virtual ~NetworkGraph() = default;

  NetworkGraph &operator=(const NetworkGraph &) = default;
  NetworkGraph &operator=(NetworkGraph &&) = default;

  virtual const std::string label() const = 0;
  virtual size_t rank_n() const override { return ranks.size(); }
  virtual size_t switch_n() const override {
    size_t switch_n = 0;
    this->visit([&](const NetworkNode<W> &n) {
      if (is_switch(n)) ++switch_n;
    });
    return switch_n;
  }
  virtual size_t compute_n() const override { return this->node_n() - switch_n(); }
  virtual size_t link_n() const override {
    size_t result = 0;
    this->visit_edges([&](int src, int dst, const Edge &e) { result += e.link_n; });
    return result;
  }

  // number of hops between two nodes
  virtual int hop_dist(int src_node_id, int dst_node_id) const {
    // default to inefficient implementation
    // override for better performance
    return get_path(src_node_id, dst_node_id).size() - 1; // path has inclusive bounds
  };
  // find a neighborhood of nodes around a rank
  virtual std::pair<std::vector<int>, double>
  neighborhood(int src_rank_id, std::vector<std::pair<int, int>> rank_n_weights) const override;
  virtual void set_hood(const std::vector<int> *hp) override { hood_ptr = hp; }
  virtual const std::vector<int> & get_hood() const override { return *hood_ptr; }

  // stats functions
  virtual void append_stats(std::string filename, std::string mapper, size_t map_ms) const override;
  virtual void print_stats() const override;
  virtual NetworkStats get_stats() override { return stats; }

protected:

  void add_rank(NetworkNode<W> &n, int rank) {
    RankPtr rp{new Rank{rank, n}};
    n.ranks.push_back(rp);
    this->ranks[rank] = std::move(rp);
  }
  const NetworkNode<W> &node_at_rank(int rank) const {
    return ranks.at(rank)->node;
  }
  NetworkNode<W> &node_at_rank(int rank) {
    return ranks.at(rank)->node;
  }

  // return vector of RankWeights with utilizations
  // override if we want to return the nodes in a locality-aware custom order
  virtual std::vector<RankWeight<W>> rank_weights() const {
    std::vector<RankWeight<W>> rw_vec;
    if (flag_debug_rank_order) {
      this->visit([&](const NetworkNode<W> &n) {
        for (const auto &r : n.ranks) {
          rw_vec.push_back(RankWeight<W>{r->rank, r->util});
        }
      });
    } else {
      for (const auto &p : ranks) {
        int rank = p.first;
        const auto &r = p.second;
        rw_vec.push_back(RankWeight<W>{rank, r->util});
      }
    }
    return rw_vec;
  }

private:
  mutable NetworkStats stats;
  std::unordered_map<int, RankPtr> ranks;
  const std::vector<int> *hood_ptr;

  virtual bool is_switch(int id) const = 0;
  virtual int get_switch(int id) const = 0;
  virtual bool is_switch(const NetworkNode<W> &n) const { return is_switch(n.id); }
  virtual int get_switch(const NetworkNode<W> &n) const { return get_switch(n.id); }

  // bias is the edge congestion bias
  void reset(W rank_cap, double bias) {
    stats.reset();
    for (auto &p: ranks) {
      p.second->cap = rank_cap;
      p.second->util = 0.0;
      p.second->comp_byte_n = 0;
      p.second->on_rank_comm_msg_n = 0;
      p.second->on_rank_comm_byte_n = 0;
      p.second->on_node_comm_msg_n = 0;
      p.second->on_node_comm_byte_n = 0;
      p.second->off_node_comm_msg_n = 0;
      p.second->off_node_comm_byte_n = 0;
    }
    this->visit([=](NetworkNode<W> &n) {
      n.router_msg_n = 0;
      n.router_byte_n = 0;
    });
    this->visit_edges([=](int src, int dst, Edge &e) {
      e.byte_n = 0;
      e.congest_bias = bias;
    });
    _space_heap.clear();
  }

  // increment the rank_cap of all ranks indicated in rank_mask by cap_inc
  void inc_rank_caps(RankMask rank_mask, W cap_inc) {
    int num_set = 0;
    for (auto id : rank_mask.node_ids) {
      if (num_set == rank_mask.rank_n) break;
      for (auto &r : this->get_node(id).ranks) {
        r->cap += cap_inc;
        if (++num_set == rank_mask.rank_n) break;
      }
    }
  }

  int find_space(W load, int start_id); // find rank with space for load
  void add_load(int rank, W load); // add load to a rank
  // add traffic between two ranks
  void add_traffic(int src_rank, int dst_rank, size_t msg_n, size_t byte_n);

  // return vector of RankWeights using order from vector of node IDs
  std::vector<RankWeight<W>> rank_weights(std::vector<int> node_ids) const {
    std::vector<RankWeight<W>> rw_vec;
    for (auto id : node_ids) {
      const NetworkNode<W> &n = this->get_node(id);
      for (const auto &r : n.ranks) {
        rw_vec.push_back(RankWeight<W>{r->rank, r->util});
      }
    }
    return rw_vec;
  }

  // return vector of RankWeightCaps using order from vector of node IDs
  std::vector<RankWeightCap<W>> rank_weight_caps(std::vector<int> node_ids) const {
    std::vector<RankWeightCap<W>> rwc_vec;
    for (auto id : node_ids) {
      const NetworkNode<W> &n = this->get_node(id);
      for (const auto &r : n.ranks) {
        rwc_vec.push_back(RankWeightCap<W>{r->rank, r->util, r->cap});
      }
    }
    return rwc_vec;
  }

  // return vector of RankWeights using a RankMask
  std::vector<RankWeight<W>> rank_weights(RankMask mask) const {
    std::vector<RankWeight<W>> rw_vec;
    bool done = (mask.rank_n == 0);
    for (auto id : mask.node_ids) {
      if (done) break;
      const NetworkNode<W> &n = this->get_node(id);
      for (const auto &r : n.ranks) {
        rw_vec.push_back(RankWeight<W>{r->rank, r->util});
        if ((int) rw_vec.size() == mask.rank_n) {
          done = true;
          break;
        }
      }
    }
    assert((int) rw_vec.size() == mask.rank_n);
    return rw_vec;
  }

  CSRMat<NodeW, EdgeW> get_csr_no_switch() const {
    // filter and project onto switches
    auto is_compute = [&](const NetworkNode<W> &n) { return !is_switch(n); };
    auto projection = [&](int node_id) { return get_switch(node_id); };
    return this->get_csr(is_compute, projection);
  }

  // max space heap
  struct SpaceHeap {
    std::vector<RankPtr> heap{};
    static bool rank_space_lt(const RankPtr &a, const RankPtr &b) {
      return static_cast<double>(a->space()) < static_cast<double>(b->space());
    };
    void clear() { heap.clear(); }
    void stage(RankPtr p) { heap.push_back(p); }
    void reheapify() { std::make_heap(heap.begin(), heap.end(), rank_space_lt); }
    bool is_init() const { return heap.size() > 0; }
    void add_load_to_head(W load) {
      pop_heap(heap.begin(), heap.end(), rank_space_lt);
      heap.back()->util += load;
      push_heap(heap.begin(), heap.end(), rank_space_lt);
    }
    void sort() { std::sort(heap.begin(), heap.end(), rank_space_lt); }
  };
  SpaceHeap _space_heap;

  void init_space_heap() {
    _space_heap.clear();
    if (flag_debug_rank_order) {
      this->visit([&](const NetworkNode<W> &n) {
        for (const auto &r : n.ranks) {
          _space_heap.stage(r);
        }
      });
      _space_heap.sort();
    } else {
      for (const auto &p : ranks) { _space_heap.stage(p.second); }
      _space_heap.reheapify();
    }
  }
  bool use_space_heap() const { return _space_heap.is_init(); }
  int space_heap_get_max() const { return _space_heap.heap.front()->rank; }

  // node struct for dijkstra heap
  struct HE {
    int id;
    double dist;
    int prev; // pointer backwards
    static bool comp(const HE &a, const HE &b) {
      return a.dist > b.dist; // min-heap
    }
  };

  std::string he_str(const HE &he) const {
    std::ostringstream oss;
    oss << "(" << this->node_id_str(he.id) << ", " << he.dist
        << ", " << this->node_id_str(he.prev) << ")";
    return oss.str();
  }

  void push_path(Path &path, int start_id, int end_id,
                 const std::unordered_map<int, HE> &visited) const {
    if (end_id != start_id) {
      push_path(path, start_id, visited.at(end_id).prev, visited);
    }
    path.push_back(end_id);
    return;
  }

  // find the closest node that satisfies the condition
  // append result onto path passed in
  template <class CondF, class EdgeFiltF>
  void dijkstra(Path &path, int start_id, CondF cond_f, EdgeFiltF edge_filt_f) const;

  // convenience function: create and return the path found
  template <class CondF, class EdgeFiltF>
  Path dijkstra(int start_id, CondF cond_f, EdgeFiltF edge_filt_f) const {
    Path path;
    dijkstra(path, start_id, cond_f, edge_filt_f);
    return path;
  }

  Path get_path(int start_id, int end_id) const;
  virtual void switch_path(Path &path, int start_id, int end_id) const;

  void print(std::ostream &os) const override {
    this->visit([&](const NetworkNode<W> &n) {
      os << n << "\n";
      for (const auto &p2 : n.edges) {
        int src = n.id, dst = p2.first;
        const Edge &e = p2.second;
        os << "  Edge (" << src << ", " << dst << "): ("
           << e.link_n << ", " << e.edge_byte_n() << ")\n";
      }
    });
  }

  // stats helpers
  std::tuple<W,W,W> weight_min_avg_max() const;
  template <class T, class HistF, class FiltF>
  std::vector<std::size_t> make_edge_hist(size_t bin_n, T min, T max, HistF hf, FiltF ff) const;
  void compute_edge_byte_stats() const;

  void print_edges() const;
  void print_hist() const;

  // cost estimate helpers
  double roofline_node_cost_fn(const NetworkNode<W> &n) const;
  double edge_cost_fn(const Edge &e) const;
  double router_cost_fn(const Node &n) const;
  double estimate_cost() const;
};

// find the closest node that satisfies the condition
// append result onto path passed in
// TODO: randomize between all nodes with the same shortest path distance
template <class W>
template <class CondF, class EdgeFiltF>
void NetworkGraph<W>::dijkstra(Path &path, int start_id, CondF cond_f, EdgeFiltF edge_filt_f) const {

  // quick check if start satisfies the condition
  if (cond_f(this->get_node(start_id))) {
    path.push_back(start_id);
    return;
  }

  std::unordered_set<int> unvisited;
  std::unordered_map<int, HE> visited;

  // q contains min-heap of distances to nodes from start
  // and may contain duplicates and already visited nodes
  std::vector<HE> q;
  q.push_back(HE{start_id, 0, kNone});

  // run Dijkstra's algorithm variant with redundant heap entries
  // (generally faster than decrease-key algorithms for sparse graphs)
  unvisited.insert(start_id);
  while (unvisited.size()) {
    // get (and remove) closest node from q
    assert(q.size());
    pop_heap(q.begin(), q.end(), HE::comp);
    HE cur_node(std::move(q.back()));
    q.pop_back();
    if (flag_verbose_dijkstra) Say() << "  Popped " << he_str(cur_node);

    // if cur_node has not been visited
    if (unvisited.count(cur_node.id)) {
      if (flag_verbose_dijkstra) Say() << "  Visiting " << he_str(cur_node);
      // mark visited
      unvisited.erase(cur_node.id);
      visited[cur_node.id] = cur_node;
      // check to see if it satisfies our condition
      if (cond_f(this->get_node(cur_node.id))) {
        push_path(path, start_id, cur_node.id, visited);
        return;
      }
      // update distances to neighbors
      for (const auto &p : this->get_node(cur_node.id).edges) {
        if (!edge_filt_f(p.second)) continue;
        int neighbor_id = p.first;
        if (visited.count(neighbor_id) == 0) {
          unvisited.insert(neighbor_id);
          const Edge &e = p.second;
          double new_dist = cur_node.dist + e.congest_biased();
          HE neighbor{neighbor_id, new_dist, cur_node.id};
          if (flag_verbose_dijkstra) Say() << "  Pushing " << he_str(neighbor);
          q.push_back(std::move(neighbor));
          push_heap(q.begin(), q.end(), HE::comp);
        }
      }
    }
  }

  // condition never satisfied
  return;
}

template <class W>
int NetworkGraph<W>::find_space(W load, int start_id) {
  // check node for a rank with space and save rank to captured variable
  int rank = kNone;
  auto has_space = [&](const NetworkNode<W> &n) {
    for (const auto &r : n.ranks) {
      if (load <= r->space()) {
        rank = r->rank;
        return true;
      }
    }
    return false;
  };
  auto f_true = [](const Edge &e) { return true; };

  // check for easy case
  if (start_id != kNone && has_space(this->get_node(start_id))) {
    return rank; // this value is set from inside has_space
  }

  // use dijkstra from starting point, if one is given
  if (start_id != kNone) {
    if (flag_verbose_greedy) Say() << "Running Dijkstra's SSSP algorithm ...";
    Path path = dijkstra(start_id, has_space, f_true);
    if (rank != kNone) return rank;
  }

  // fallback if no start_id given or dijkstra failed
  if (flag_verbose_greedy) Say() << "Running linear search ...";
  this->visit_until(has_space);
  if (rank != kNone) return rank;

  // did not find a rank with sufficient space
  if (flag_verbose_greedy) {
    Say() << "WARNING: could not find rank with sufficient capacity!";
  }
  return kNone;
}

template <class W>
void NetworkGraph<W>::add_load(int rank, W load) {
  stats.active_nodes.insert(node_at_rank(rank).id);
  if (use_space_heap()) {
    if (rank == space_heap_get_max()) {
      // does the load add for us
      _space_heap.add_load_to_head(load);
    } else {
      ranks.at(rank)->util += load;
      _space_heap.reheapify();
    }
  } else {
    // do it ourselves
    ranks.at(rank)->util += load;
  }
}

template <class W>
typename NetworkGraph<W>::Path
NetworkGraph<W>::get_path(int start_id, int end_id) const {
  int start_switch = get_switch(start_id), end_switch = get_switch(end_id);
  if (flag_verbose_routing) Say() << "Path from " << this->node_id_str(start_id)
                          << " to " << this->node_id_str(end_id);
  Path path;
  if (start_id != start_switch) {
    path.push_back(start_id);
  }
  switch_path(path, start_switch, end_switch);
  if (flag_verbose_routing) {
    for (int switch_id : path) Say() << "  Traversing to " << this->node_id_str(switch_id);
  }
  if (end_id != end_switch) {
    if (flag_verbose_routing) Say() << "  Traversing to " << this->node_id_str(end_id);
    path.push_back(end_id);
  }
  return path;
}

// default to using dijkstra to find shortest path
// can override for custom routing protocol
template <class W>
void NetworkGraph<W>::switch_path(Path &path, int start_id, int end_id) const {
  dijkstra(path, start_id,
           [=](const NetworkNode<W> &n) { return n.id == end_id; },
           [=](const Edge &e) { return !e.node_link; });
}

template <class W>
void NetworkGraph<W>::add_traffic(int src_rank_id, int dst_rank_id, size_t msg_n, size_t byte_n) {
  
  // all traffic
  stats.total_msg_n += msg_n;
  stats.total_byte_n += byte_n;

  Rank &src_rank = *ranks.at(src_rank_id);
  Rank &dst_rank = *ranks.at(dst_rank_id);
  Node &src_node = src_rank.node;
  Node &dst_node = dst_rank.node;
  int src_node_id = src_node.id;
  int dst_node_id = dst_node.id;

  if (src_node_id == dst_node_id) {
    if (src_rank_id == dst_rank_id) { // same rank
      src_rank.on_rank_comm_msg_n += msg_n;
      src_rank.on_rank_comm_byte_n += byte_n;
      stats.on_rank_msg_n += msg_n;
      stats.on_rank_byte_n += byte_n;
    } else { // same node, different ranks
      src_rank.on_node_comm_msg_n += msg_n;
      src_rank.on_node_comm_byte_n += byte_n;
      stats.on_node_msg_n += msg_n;
      stats.on_node_byte_n += byte_n;
    }
  } else { // different nodes
    src_rank.off_node_comm_msg_n += msg_n;
    src_rank.off_node_comm_byte_n += byte_n;
    stats.off_node_msg_n += msg_n;
    stats.off_node_byte_n += byte_n;

    auto path = get_path(src_node_id, dst_node_id);
    unsigned hop_n = stats.flag_switch_hops_only ? path.size()-3 : path.size()-1;
    stats.off_node_msg_hop_n += hop_n;
    stats.off_node_byte_hop_n += byte_n * hop_n;

    // add byte_n to network links along path
    int prev = src_node_id;
    for (auto cur : path) {
      this->get_node(cur).router_msg_n += msg_n;
      this->get_node(cur).router_byte_n += byte_n;
      if (cur == src_node_id) continue;
      this->get_edge(prev, cur).byte_n += byte_n;
      prev = cur;
    }
  }
}

inline void NetworkStats::write_stats_header(std::string filename) {
  std::ofstream of(filename);
  of << "Topology" << '\t';
  of << "Ranks" << '\t';
  of << "Mapper" << '\t';

  of << "Mapping Time (ms)" << '\t';

  of << "Num. Active Nodes" << '\t';

  of << "Total Bytes Sent" << '\t';
  of << "Total Bytes On-rank" << '\t';
  of << "Total Bytes On-node" << '\t';
  of << "Total Bytes Off-node" << '\t';
  of << "Total Byte-Hops" << '\t';
  of << "Avg. Hops / Byte" << '\t';

  of << "Total Messages Sent" << '\t';
  of << "Total Messages On-rank" << '\t';
  of << "Total Messages On-node" << '\t';
  of << "Total Messages Off-node" << '\t';
  of << "Total Message-Hops" << '\t';
  of << "Avg. Hops / Msg" << '\t';
  
  of << "Max/Avg. Load" << '\t';

  of << "Active Edges" << '\t';
  of << "Avg. Bytes / Act. Edge" << '\t';
  of << "Max Edge Bytes" << '\t';

  of << "Active Links" << '\t';
  of << "Avg. Bytes / Act. Link" << '\t';
  of << "Max Link Bytes" << '\t';

  of << '\n';
}

template <class W>
void NetworkGraph<W>::append_stats(std::string filename, std::string mapper, size_t map_ms) const {
  std::ofstream of(filename, std::ofstream::app);
  of << label() << '\t';
  of << rank_n() << '\t';
  of << mapper << '\t';

  of << map_ms << '\t';

  of << stats.active_nodes.size() << '\t';

  of << stats.total_byte_n << '\t';
  of << stats.on_rank_byte_n << '\t';
  of << stats.on_node_byte_n << '\t';
  of << stats.off_node_byte_n << '\t';
  of << stats.off_node_byte_hop_n << '\t';
  of << stats.avg_off_node_byte_hop_n() << '\t';

  of << stats.total_msg_n << '\t';
  of << stats.on_rank_msg_n << '\t';
  of << stats.on_node_msg_n << '\t';
  of << stats.off_node_msg_n << '\t';
  of << stats.off_node_msg_hop_n << '\t';
  of << stats.avg_off_node_msg_hop_n() << '\t';

  W min_wt, avg_wt, max_wt;
  std::tie(min_wt, avg_wt, max_wt) = weight_min_avg_max();
  of << max_wt / avg_wt << '\t';

  of << stats.active_edge_n << '\t';
  of << stats.avg_active_edge_byte_n() << '\t';
  of << stats.max_edge_byte_n << '\t';

  of << stats.active_link_n << '\t';
  of << stats.avg_active_link_byte_n() << '\t';
  of << stats.max_link_byte_n << '\t';

  of << '\n';
}

template <class W>
std::tuple<W,W,W> NetworkGraph<W>::weight_min_avg_max() const {
  auto vec = this->rank_weights();
  W sum = 0;
  W cur_min = std::numeric_limits<W>::max();
  W cur_max = std::numeric_limits<W>::min();
  for (const RankWeight<W> &rw: vec) {
    cur_min = min(cur_min, rw.weight);
    cur_max = max(cur_max, rw.weight);
    sum += rw.weight;
  }
  W avg = sum / W{static_cast<double>(rank_n())};
  return std::tie(cur_min, avg, cur_max);
}

template <class W>
template <class T, class HistF, class FiltF>
std::vector<std::size_t>
NetworkGraph<W>::make_edge_hist(size_t bin_n, T min, T max, HistF hf, FiltF ff) const {
  std::vector<std::size_t> result;
  result.resize(bin_n); // bins between min and max
  T bin_sz = (max-min)/bin_n;
  // return the bin corresponding to value
  auto bin = [=](T val)->size_t {
    size_t raw_bin = val < min ? 0 : (val-min)/bin_sz;
    return std::min(bin_n-1, raw_bin);
  };
  this->visit_edges([&](int src, int dst, const Edge &e) {
    if (ff(src, dst, e)) { ++result[bin(hf(src, dst, e))]; }
  });
  return result;
}

template <class W>
void NetworkGraph<W>::print_edges() const {
  std::vector<std::tuple<int,int,size_t>> sorted_edges;
  this->visit_edges([&](int src, int dst, const Edge &e) {
    if (!is_switch(dst) && e.edge_byte_n()) {
      sorted_edges.push_back(std::tuple<int,int,size_t>{src, dst, e.edge_byte_n()});
    }
  });
  std::sort(sorted_edges.begin(), sorted_edges.end());
  for (auto &e : sorted_edges) {
    int src, dst; size_t edge_byte_n; std::tie(src, dst, edge_byte_n) = e;
    Say() << "Switch-node (" << src << "," << dst << "): " << edge_byte_n;
  }
}
template <class W>
void NetworkGraph<W>::print_hist() const {
  // filter out edges that have no traffic
  auto filter = [&](int src, int dst, const Edge &e) {
    return e.edge_byte_n() > 0;
  };
  auto edge_byte_hist = make_edge_hist<size_t>(100, 0, stats.max_edge_byte_n,
    [=](int src, int dst, const Edge &e) { return e.edge_byte_n(); },
    filter);
  Say() << "Edge-Bytes Histogram (min, max)=(0, " << stats.max_edge_byte_n << ")";
  for (auto val : edge_byte_hist) { Say() << "  " << val; }

  auto link_byte_hist = make_edge_hist<size_t>(100, 0, stats.max_link_byte_n,
    [=](int src, int dst, const Edge &e) { return e.link_byte_n(); },
    filter);
  Say() << "Link-Bytes Histogram (min, max)=(0, " << stats.max_link_byte_n << ")";
  for (auto val : link_byte_hist) { Say() << "  " << val; }
}

template <class W>
void NetworkGraph<W>::compute_edge_byte_stats() const {
  stats.max_edge_byte_n = stats.max_link_byte_n = 0;
  stats.active_edge_n = stats.active_link_n = 0;
  this->visit_edges([&](int src, int dst, const Edge &e) {
    if (!stats.flag_switch_hops_only || (is_switch(src) && is_switch(dst))) {
      stats.max_link_byte_n = std::max(e.link_byte_n(), stats.max_link_byte_n);
      stats.max_edge_byte_n = std::max(e.edge_byte_n(), stats.max_edge_byte_n);
      if (e.edge_byte_n()) { stats.active_edge_n++; stats.active_link_n += e.link_n; }
    }
  });
}

template <class W>
void NetworkGraph<W>::print_stats() const {
  Say() << "Mapped to " << rank_n() << " ranks";
  W min_wt, avg_wt, max_wt;
  std::tie(min_wt, avg_wt, max_wt) = weight_min_avg_max();
  Say() << "Rank load max          : " << max_wt;
  Say() << "Rank load avg          : " << avg_wt;
  Say() << "Rank load max/avg      : " << max_wt / avg_wt;

  Say() << "Total bytes sent       : " << stats.total_byte_n;
  Say() << "Total bytes on-rank    : " << stats.on_rank_byte_n;
  Say() << "Total bytes on-node    : " << stats.on_node_byte_n;
  Say() << "Total bytes off-node   : " << stats.off_node_byte_n;
  Say() << "Total byte-hops        : " << stats.off_node_byte_hop_n;
  Say() << "Avg.  byte-hops        : " << stats.avg_off_node_byte_hop_n();

  Say() << "Total messages sent    : " << stats.total_msg_n;
  Say() << "Total messages on-rank : " << stats.on_rank_msg_n;
  Say() << "Total messages on-node : " << stats.on_node_msg_n;
  Say() << "Total messages off-node: " << stats.off_node_msg_n;
  Say() << "Total message-hops     : " << stats.off_node_msg_hop_n;
  Say() << "Avg.  message-hops     : " << stats.avg_off_node_msg_hop_n();

  compute_edge_byte_stats();

  Say() << "Total active edges     : " << stats.active_edge_n;
  Say() << "Avg. active edge-bytes : " << stats.avg_active_edge_byte_n();
  Say() << "Max edge-bytes         : " << stats.max_edge_byte_n;

  Say() << "Total active links     : " << stats.active_link_n;
  Say() << "Avg. active link-bytes : " << stats.avg_active_link_byte_n();
  Say() << "Max link-bytes         : " << stats.max_link_byte_n;

  auto cost = estimate_cost();
  Say() << "Estimated cost: " << cost;

  if (flag_verbose_stats) {
    Say() << "All Weights (rank, weight):";
    auto vec = this->rank_weights();
    // sort by rank
    std::sort(vec.begin(), vec.end(),
      [](const RankWeight<W> &a, const RankWeight<W> &b) {
        return a.rank < b.rank;
    });
    size_t zero_n{0};
    for (const RankWeight<W> &rw: vec) {
      if (rw.weight == 0) {
        ++zero_n;
      }
      Say() << rw.rank << '\t' << rw.weight;
    }
    Say() << "  Num zero weight: " << zero_n;
  }
  Say() << ""; // blank line
}

// compute roofline cost using flops and DRAM traffic
template <class W>
double NetworkGraph<W>::roofline_node_cost_fn(const NetworkNode<W> &n) const {

  double max_rank_comp    = 0; // Flops cost due to comp
  double dram_comp_byte_n = 0; // DRAM traffic due to compute
  double dram_msg_byte_n  = 0; // DRAM traffic due to on-rank/on-node messages
  double dram_msg_n       = 0; // number of on-node messages sent

  for (RankPtr rp: n.ranks) {

    // flops due to compute
    max_rank_comp = std::max(max_rank_comp, (double) rp->util);
    dram_comp_byte_n += rp->comp_byte_n;

    // DRAM traffic due to on-rank and on-node messages
    dram_msg_byte_n += rp->on_rank_comm_byte_n;
    dram_msg_n      += rp->on_rank_comm_msg_n;
    dram_msg_byte_n += rp->on_node_comm_byte_n;
    dram_msg_n      += rp->on_node_comm_msg_n;
  }

  const double comp_coeff = 1.0,    // convert to seconds
               dram_bw = 500e9,     // 500 GB/s aggregate node DRAM BW
               dram_msg_oh = 100e-6; // 100 us msg send overhead

  auto cpu_cost = max_rank_comp * comp_coeff;
  auto dram_cost = (dram_comp_byte_n + dram_msg_byte_n) / dram_bw +
                   dram_msg_n * dram_msg_oh;
  return std::max(cpu_cost, dram_cost);
}

template <class W>
double NetworkGraph<W>::edge_cost_fn(const Edge &e) const {
  const double link_bw = 120e9; // 100 GB/s network link bandwidth
  return e.link_byte_n() / link_bw;
}

template <class W>
double NetworkGraph<W>::router_cost_fn(const Node &n) const {
  const double router_bw = 80e9; // 500 GB/s router bandwidth
  return n.router_byte_n / router_bw;
}

template <class W>
double NetworkGraph<W>::estimate_cost() const {

  double max_node_cost = 0;
  double max_router_cost = 0;
  double max_edge_cost = 0;

  this->visit([&](const Node &n) {
    // amount of network traffic handled by this node
    // injection or routing if node is compute or switch, resp.
    auto cur_router_cost = router_cost_fn(n);
    max_router_cost = std::max(max_router_cost, cur_router_cost);

    // handle compute nodes
    if (!is_switch(n.id)) {
      // node's local cost (compute cost + on-node messaging + off-node messaging)
      auto cur_node_cost = roofline_node_cost_fn(n) + cur_router_cost;
      max_node_cost = std::max(max_node_cost, cur_node_cost);
    }
  });

  this->visit_edges([&](int src, int dst, const Edge &e) {
    auto cur_edge_cost = edge_cost_fn(e);
    max_edge_cost = std::max(max_edge_cost, cur_edge_cost);
  });

  Say() << "Max node   cost: " << max_node_cost;
  Say() << "Max edge   cost: " << max_edge_cost;
  Say() << "Max router cost: " << max_router_cost;

  return std::max(std::max(max_node_cost, max_router_cost), max_edge_cost);
}

////////////////////////////////////////////////////////////
// Find a neighborhood of close nodes around a given rank //
////////////////////////////////////////////////////////////

namespace _ {
  inline int pair_n(int x) {
    return x*(x-1)/2;
  }
}

// find close neighborhood of nodes around src rank id with at least want_rank_n ranks
// src_rank_id: where to start searching
// rank_n_weights: vector of pairs: (rank_n, weight) used to calculate the overall score
//                 for each rank count, multiply its score by the corresponding weight
// return value: (node ids, score)
template <class W>
std::pair<std::vector<int>, double>
NetworkGraph<W>::neighborhood(int src_rank_id, std::vector<std::pair<int, int>> rank_n_weights) const
{
  std::vector<int> result; // node ids
  double score = 0;

  std::sort(rank_n_weights.begin(), rank_n_weights.end());
  assert(rank_n_weights.size() == 0 ||
         (1 < rank_n_weights.front().first &&
          rank_n_weights.back().first <= (int) rank_n()));

  int want_rank_n = rank_n();
  const auto &src_node = node_at_rank(src_rank_id);
  int cur_node_id = src_node.id;
  int cur_node_rank_n = src_node.rank_n();

  // add source node
  result.push_back(cur_node_id);
  int total_rank_n = cur_node_rank_n;
  int total_pairs_dist = 0;
  if (flag_verbose_neighborhood) {
    Say() << "  Added " << this->node_id_str(cur_node_id)
          << ", total ranks: " << total_rank_n
          << ", avg dist: " << 0;
  }
  if (total_rank_n >= want_rank_n) {
    return {result, 0};
  }

  // set score_idx to point to first entry with more ranks than total so far
  unsigned score_idx = 0;
  while (score_idx < rank_n_weights.size() &&
         rank_n_weights[score_idx].first <= total_rank_n) {
    score_idx++;
  }

  // find set of node IDs that have at least 1 rank
  // create a map : node_id -> (node_rank_n, sum_dist)
  // where sum_dist is sum of pairwise rank distances from the candidate node
  //   to already mapped nodes.
  std::unordered_map<int, std::pair<int, int>> candidates;
  for (const auto &p : ranks) {
    const auto &node = p.second->node;
    if (candidates.count(node.id) == 0) {
      candidates[node.id] = {node.rank_n(), 0};
    }
  }

  if (flag_verbose_neighborhood) {
    Say() << "  Candidates:";
    for (const auto &p : candidates) {
      Say() << "    " << this->node_id_str(p.first) << " : " << p.second.first << " ranks";
    }
  }

  while (total_rank_n < want_rank_n) {
    // cur_node_id was just added
    candidates.erase(cur_node_id);
    double min_avg_dist = std::numeric_limits<double>::max();
    int next_node_id = 0;
    // update candidates with their pairwise rank distances to cur_node
    for (auto &p : candidates) {
      const int cand_node_id = p.first;
      const int cand_node_rank_n = p.second.first;
      int &sum_dist = p.second.second;
      // multiply hop distance by number of rank pairs across the two nodes
      sum_dist += hop_dist(cand_node_id, cur_node_id) * (cand_node_rank_n * cur_node_rank_n);
      double avg_dist = static_cast<double>(sum_dist + total_pairs_dist) /
                        _::pair_n(total_rank_n + cand_node_rank_n);
      if (flag_verbose_neighborhood) {
        Say() << "    Hop Distance from " << this->node_id_str(cand_node_id)
              << " to " << this->node_id_str(cur_node_id)
              << ": " << hop_dist(cand_node_id, cur_node_id)
              << ", candidate avg: " << avg_dist;
      }
      // keep track of what should be the next node to add
      if (avg_dist < min_avg_dist) {
        next_node_id = cand_node_id;
        min_avg_dist = avg_dist;
      }
    }
    cur_node_id = next_node_id;
    cur_node_rank_n = candidates.at(cur_node_id).first;
    total_pairs_dist += candidates.at(cur_node_id).second;
    // add cur_node_id to result
    result.push_back(cur_node_id);
    total_rank_n += cur_node_rank_n;
    if (flag_verbose_neighborhood) {
      Say() << "  Added " << this->node_id_str(cur_node_id)
            << ", ranks: " << cur_node_rank_n
            << ", total ranks: " << total_rank_n
            << ", avg dist: " << min_avg_dist;
    }
    // contribute to score
    while (score_idx < rank_n_weights.size() &&
           rank_n_weights[score_idx].first <= total_rank_n) {
      score += rank_n_weights[score_idx].second * min_avg_dist;
      if (flag_verbose_neighborhood) {
        Say() << "  Contributed to score (" << rank_n_weights[score_idx].first
              << ", " << rank_n_weights[score_idx].second
              << " * " << min_avg_dist << ") -> " << score;
      }
      score_idx++;
    }
  }
  assert(score_idx == rank_n_weights.size());

  return {result, score};
}

}

#endif // NETGRAPH_HXX
