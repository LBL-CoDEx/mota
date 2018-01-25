#include "random.hxx"
#include "util.hxx"

namespace mota {
namespace random {

rand_engine_t rnd_eng;
bool rand_is_init = false;

void init(int rank_here) {
  if (!rand_is_init) {
    // set up random number generator
    auto rnd_def_seed = rand_engine_t::default_seed;
    auto seed = env<rand_result_t>("rnd_seed", rnd_def_seed);
    rnd_eng.seed(seed + rank_here);
    rnd_eng();
    rnd_eng();
    rand_is_init = true;
  }
}

}}
