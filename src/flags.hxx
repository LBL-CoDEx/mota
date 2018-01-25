#ifndef FLAGS_HXX
#define FLAGS_HXX

namespace mota {

// can set to false to speed things up
const bool flag_force_traffic = false;
// can set to true to speed things up
const bool flag_fast_add_traffic = false;

const bool flag_verbose_assign = false;
const bool flag_verbose_dijkstra = false;
const bool flag_verbose_neighborhood = false;
const bool flag_verbose_greedy = false;
const bool flag_verbose_stats = false;
const bool flag_verbose_distribute = false;
const bool flag_verbose_csr = false;
const bool flag_verbose_metis = false;
const bool flag_verbose_network = false;
const bool flag_verbose_routing = false;

const bool flag_debug_rank_order = false;

extern bool flag_quiet;

}

#endif // FLAGS_HXX
