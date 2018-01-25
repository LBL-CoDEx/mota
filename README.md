## Mota Mapper: A Multi-Objective Topology-Aware Mapper Library

### Prerequisites: 
- C++11 compliant compiler
- External codes:
  - Get morton ordering algorithms as in Fabian Giesen's blog: http://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/

### Network Graph:
- To create a network graph, call factory funcction `net_graph_create(std::string label, int cons_n, const RanksByCartesian<4> &ranks_by_coord)` to create a network graph of type `label`, whose nodes are assigned ranks by `ranks_by_coord`, and have capacities for `cons_n` constraints.
- See `src/netmodels.hxx` for more details

### Application Graph: 
- To create an application graph, call factory function `app_graph_create(int cons_n)`, with `cons_n`   set to the number of constraints desired (e.g. number of AMR levels that require load balancing)
- To specify the nodes and edges in the application graph, call `AppGraph_::import()`
- To map the application graph to a network graph, call `AppGraph_::run_mapper()`
- See `src/appgraph.hxx` for more details

## Copyright

"Mota Mapper (Mota)" Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved.

If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare derivative works, and perform publicly and display publicly, and to permit other to do so.
