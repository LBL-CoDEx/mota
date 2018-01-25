#ifndef RANDOM_HXX
#define RANDOM_HXX

#include <random>
#include <cassert>

namespace mota {
namespace random {

using rand_engine_t = std::default_random_engine;
using rand_result_t = rand_engine_t::result_type;

extern rand_engine_t rnd_eng;
extern bool rand_is_init;

void init(int rank_here);

template <class T>
typename T::result_type roll(T dist) {
  assert(rand_is_init);
  return dist(rnd_eng);
}

}}

#endif // RANDOM_HXX
