#ifndef GRAPH_HXX
#define GRAPH_HXX

#include <functional>
#include <string>

#include "flags.hxx"
#include "util.hxx"
#include "cartesian.hxx"
#include "csr.hxx"

/*****************
 * Generic Graph *
 *****************/

namespace mota {

// requires Node objects have members:
//   id: node ID
//   weight(): returns the weight of the node (castable to double)
//   edges: map from destination node ID to Edge object, which has a weight() member
template <class Node>
class Graph {
public:
  using NodeW   = typename Node::Weight;
  using Edge    = typename Node::Edge;
  using EdgeW   = typename Node::Edge::Weight;
  using NodePtr = typename Node::Ptr;

  Graph() = default;
  Graph(const Graph &) = default;
  Graph(Graph &&) = default;
  virtual ~Graph() {}

  Graph &operator=(const Graph &) = default;
  Graph &operator=(Graph &&) = default;

  size_t node_n() const { return nodes.size(); }
  size_t edge_n() const {
    size_t result = 0;
    visit([&](const Node &n) { result += n.edges.size(); });
    return result;
  }

  friend std::ostream &operator<<(std::ostream &os, const Graph &g) {
    g.print(os); // dynamic dispatch to print() member
    return os;
  }

  // return sparse matrix representation of this graph
  // only includes nodes for which filter returns true
  // for projection function P: A connects B iff P(A) connects P(B)
  // connect sibling nodes (P(A)==P(B)) with large weight
  template <class NodeFiltF = std::function<bool(const Node &n)>,
            class ProjectF = std::function<int(int)>>
  CSRMat<NodeW, EdgeW>
  get_csr(NodeFiltF filter = [](const Node &n) { return true; },
          ProjectF project = [](int x) {return x;},
          double sibling_weight_factor = 10.0) const {

    std::unordered_map<int, NodeW> node_wgts;
    std::unordered_map<int, std::vector<int>> unproject; // invert projection

    visit([&](const Node &n) {
      if (filter(n)) {
        node_wgts[n.id] = n.weight();
        if (flag_verbose_csr) Say() << "Adding node " << n.id;
        unproject[project(n.id)].push_back(n.id);
        if (flag_verbose_csr) Say() << "  projects to " << project(n.id);
      }
    });

    std::unordered_map<std::pair<int, int>, EdgeW> edge_wgts;
    EdgeW max_edge_weight(0);
    for (const auto &p1 : node_wgts) {
      int node_id = p1.first;
      int project_id = project(node_id);
      // connect to nodes projected to neighbors of projection
      for (const auto &p2 : get_node(project_id).edges) {
        const int dst = p2.first;
        const auto &edge = p2.second;
        EdgeW ew = edge.weight();
        if (unproject.count(dst)) { // if dst is a projection target (e.g. a switch with nodes attached)
          for (int neighbor : unproject.at(dst)) {
            if (flag_verbose_csr) Say() << "Connecting " << node_id << " to " << neighbor << " via (" << project_id << ", " << dst << ")";
            edge_wgts[std::make_pair(node_id, neighbor)] = ew;
            if (ew > max_edge_weight) max_edge_weight = ew;
          }
        }
      }
    }

    // connect nodes projected to the same place with large edge weight
    double sib_edge_weight = (max_edge_weight == 0) ? 1.0 :
                              max_edge_weight*sibling_weight_factor;
    for (const auto &p1 : node_wgts) {
      int node_id = p1.first;
      int project_id = project(node_id);
      // connect to other nodes projected to same place
      for (int neighbor : unproject[project_id]) {
        if (neighbor != node_id) {
          if (flag_verbose_csr) Say() << "Connecting " << node_id << " to " << neighbor << " via " << project_id;
          edge_wgts[std::make_pair(node_id, neighbor)] = sib_edge_weight;
        }
      }
    }

    return CSRMat<NodeW, EdgeW>{std::move(node_wgts), std::move(edge_wgts)};
  }


protected:
  
  // return reference to Node, create if needed
  Node & get_node_create(int id) {
    NodePtr &pn = nodes[id];
    if (!pn) pn = NodePtr(new Node(id));
    return *pn;
  }

  // return reference to existing Node
  Node & get_node(int id) { return *nodes.at(id); }
  const Node & get_node(int id) const { return *nodes.at(id); }

  // return reference to Edge, create if needed
  Edge & get_edge_create(int src_id, int dst_id) {
    return get_node(src_id).edges[dst_id];
  }

  // return reference to existing Edge
  Edge & get_edge(int src_id, int dst_id) {
    return get_node(src_id).edges.at(dst_id);
  }
  const Edge & get_edge(int src_id, int dst_id) const {
    return get_node(src_id).edges.at(dst_id);
  }

  // node visitors
  template <class F>
  void visit(F f) const {
    for (const auto &p : nodes) f(*p.second);
  }
  template <class F>
  void visit(F f) {
    for (auto &p : nodes) f(*p.second);
  }

  // early exit node visitors
  template <class F>
  void visit_until(F f) const {
    for (const auto &p : nodes) { if (f(*p.second)) break; }
  }
  template <class F>
  void visit_until(F f) {
    for (auto &p : nodes) { if (f(*p.second)) break; }
  }

  // edge visitors assume Node has subtype "Node::Edge"
  // and has members "int id" and "unordered_map<int, Edge> edges"
  template <class F>
  void visit_edges(F f) const {
    visit([=](const Node &n) { for (const auto &p : n.edges) { f(n.id, p.first, p.second); }});
  }
  template <class F>
  void visit_edges(F f) {
    visit([=](Node &n) { for (auto &p : n.edges) { f(n.id, p.first, p.second); }});
  }

  // get vector of node pointers
  std::vector<NodePtr> make_node_vec(std::vector<int> node_ids) {
    std::vector<NodePtr> result;
    for (auto id : node_ids) result.push_back(nodes.at(id));
    return result;
  }
  template <class FiltF>
  std::vector<NodePtr> filtered_nodes(FiltF filt_f) {
    std::vector<NodePtr> node_vec;
    for (const auto &p : nodes) {
      const NodePtr &pn = p.second;
      if (filt_f(*pn)) {
        node_vec.push_back(NodePtr{pn});
      }
    }
    return node_vec;
  }
  template <class CompF, class FiltF = std::function<bool(const Node &)>>
  std::vector<NodePtr> sorted_nodes(CompF comp_f,
                                    FiltF filt_f = [](const Node &){return true;}) {
    auto node_vec = filtered_nodes(filt_f);
    std::sort(node_vec.begin(), node_vec.end(),
      [&](const NodePtr &pa, const NodePtr &pb) { return comp_f(*pa, *pb); });
    return node_vec;
  }

  virtual std::string node_id_str(int node_id) const {
    return std::to_string(node_id);
  }

private:
  // map from node ID to NodePtr
  std::unordered_map<int, NodePtr> nodes;

  virtual void print(std::ostream &os) const = 0;
};

}

#endif // GRAPH_HXX
