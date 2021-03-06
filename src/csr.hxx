#ifndef _1b833e26_2871_4cce_8c46_dd736f8581e0
#define _1b833e26_2871_4cce_8c46_dd736f8581e0

#include <vector>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>

#include "util.hxx"

namespace mota {

template <class NodeW=double, class EdgeW=double>
class CSRMat {
public:

  // construct from components
  CSRMat(std::unordered_map<int, NodeW> nw,
         std::unordered_map<std::pair<int, int>, EdgeW> ew) :
    _nodes{}, _edges{},
    _node_wgts_map(std::move(nw)),
    _edge_wgts_map(std::move(ew))
  {
    for (const auto &p : _node_wgts_map) {
      _nodes.push_back(p.first);
    }
    std::sort(_nodes.begin(), _nodes.end());

    // use of intermediate edge_sets normalizes all edges to be bi-directional
    std::unordered_map<int, std::unordered_set<int>> edge_sets;
    for (const auto &p : _edge_wgts_map) {
      int src = p.first.first, dst = p.first.second;
      edge_sets[src].insert(dst);
      edge_sets[dst].insert(src);
    }
    for (const auto &p : edge_sets) {
      int src = p.first;
      for (int dst : p.second) {
        _edges[src].push_back(dst);
      }
      std::sort(_edges[src].begin(), _edges[src].end());
    }
  }

  size_t size() const { return _nodes.size(); }
  const std::vector<int> &node_ids() const { return _nodes; }
  NodeW node_weight(int id) const { return _node_wgts_map.at(id); }

  std::vector<int> xadj() const {
    std::vector<int> result;
    result.push_back(0);
    for (auto id : _nodes) {
      int edge_n = _edges.count(id) ? _edges.at(id).size() : 0;
      int idx = result.back() + edge_n;
      result.push_back(idx);
    }
    return result;
  }

  std::vector<int> adj() const {
    std::vector<int> result;
    // mapping from node IDs to index [0, num_nodes)
    std::unordered_map<int, int> mapping;
    for (unsigned i = 0; i < _nodes.size(); ++i) {
      mapping[_nodes[i]] = i;
    }
    for (auto src : _nodes) {
      if (_edges.count(src)) {
        for (auto dst : _edges.at(src)) {
          result.push_back(mapping.at(dst));
        }
      }
    }
    return result;
  }

  // vector of node weights in order
  std::vector<NodeW> node_wgts() const {
    std::vector<NodeW> result;
    result.reserve(size()); // pre-allocate
    for (auto id : _nodes) { // _nodes is sorted
      result.push_back(_node_wgts_map.at(id));
    }
    return result;
  }

  // returns a vector of edge weights in order
  std::vector<EdgeW> edge_wgts() const {
    std::vector<EdgeW> result;
    for (auto src : _nodes) { // _nodes is sorted
      if (_edges.count(src)) { // node might not be connected to anything
        for (auto dst : _edges.at(src)) { // _edges.at(src) is sorted
          EdgeW w{0};
          auto id = std::make_pair(src, dst);
          if (_edge_wgts_map.count(id)) {
            w += _edge_wgts_map.at(id);
          }
          id = std::make_pair(dst, src);
          if (_edge_wgts_map.count(id)) {
            w += _edge_wgts_map.at(id);
          }
          result.push_back(w);
        }
      }
    }
    return result;
  }

  std::vector<int> idxs_to_node_ids(std::vector<int> idxs) const {
    std::vector<int> result;
    for (auto idx : idxs) {
      result.push_back(_nodes[idx]);
    }
    return result;
  }

  friend std::ostream &operator<<(std::ostream &os, const CSRMat &mat) {
    os << "CSRMat (" << mat.size() << ")\n";
    auto nw = mat.node_wgts();
    auto ew = mat.edge_wgts();
    int nidx = 0, eidx = 0;
    for (auto nid : mat._nodes) {
      os << "(" << nid << "," << nw[nidx++] << ") : ";
      for (auto eid : mat._edges.at(nid)) { 
        os << "(" << eid << "," << ew[eidx++] << "),";
      }
      os << '\n';
    }
    return os;
  }

private:
  std::vector<int> _nodes; // sorted node IDs
  std::unordered_map<int, std::vector<int>> _edges;
  std::unordered_map<int, NodeW> _node_wgts_map;
  std::unordered_map<std::pair<int, int>, EdgeW> _edge_wgts_map;
};

}

#endif
