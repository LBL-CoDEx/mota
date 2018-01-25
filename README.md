## Mota Mapper: A Multi-Objective Topology-Aware Mapper Library

### Description:

Mota is a library that provides several heuristics for mapping tasks to compute ranks (processes) on a networked computer cluster or supercomputer.  It was originally designed for task placement of adaptive mesh refinement codes (AMR).  It is multi-objective in the sense that it simultaneously balances (potentially multiple types of) computational load on each rank as well as the communication traffic between the boxes.  The heuristics used for mapping include algorithms such as list assignment and space-filling curves, as well as algorithms from graph analysis such as adjacency matrix bandwidth reduction, recursive bisection, and greedy placement.

Mota has been used in conjunction with ProgrAMR and SST/Macro for the purpose of simulating AMR performance on future network interconnection topologies and is being integrated into the AMReX framework for box placement during dynamic regridding.  We utilize switch-link topology models for current and future supercomputer network interconnects to do task placement and make a detailed evaluation of the network utilization and performance.

### Prerequisites: 
- C++11 compliant compiler
- External codes:
  - [METIS graph partitioning library](http://glaros.dtc.umn.edu/gkhome/metis/metis/download)
  - [Morton ordering algorithms](http://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/)

### Network Graph:
- To create a network graph, call factory funcction `net_graph_create(std::string label, int cons_n, const RanksByCartesian<4> &ranks_by_coord)` to create a network graph of type `label`, whose nodes are assigned ranks by `ranks_by_coord`, and have capacities for `cons_n` constraints.
- See `src/netmodels.hxx` for more details

### Application Graph: 
- To create an application graph, call factory function `app_graph_create(int cons_n)`, with `cons_n`   set to the number of constraints desired (e.g. number of AMR levels that require load balancing)
- To specify the nodes and edges in the application graph, call `AppGraph_::import()`
- To map the application graph to a network graph, call `AppGraph_::run_mapper()`
- See `src/appgraph.hxx` for more details

## Copyright

"Mota Mapper (Mota)" Copyright (c) 2018, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved.

If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare derivative works, and perform publicly and display publicly, and to permit other to do so.
