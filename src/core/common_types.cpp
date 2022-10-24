/* 
 * File:   common_types.cpp
 * Author: kirill
 * 
 * Created on July 29, 2013, 9:56 AM
 */

#include "common_types.h"

constexpr charge_balance cb_names::__cbt[__psize];
constexpr const char * cb_names::__cbn[__psize];

template <class R = std::minstd_rand>
std::random_device::result_type new_seed(std::random_device::result_type seed,
                                         int seed_num) {
  std::random_device::result_type result = seed;
  R gen(seed);
  std::uniform_int_distribution<std::random_device::result_type>
      u(std::random_device::min(), std::random_device::max());
  for(int i = 0; i < seed_num; i++) {
    result = u(gen);
  }
  return result;
}

void rnd_indexer_t::set_properties(int random_count, int64_t total_comb, int symm_op) {
  this->random_count = random_count;
  this->total_comb = total_comb;
  if( random_count == 0) {
    method = sampling_method_t::DISABLED;
    return;
  }
  total_samples = 4 * symm_op * (random_count + 1000);
  if( (total_comb < 20000) || (2 * total_samples > total_comb) ) {
    method = sampling_method_t::ALL;
    return;
  }
  method = sampling_method_t::BY_INDEXES;
  base_index = 0;
  index_step = 2 * total_comb / symm_op + 1;
  create_indexes();
}

void rnd_indexer_t::create_indexes() {
  rnd_engine_t rnd(initial_seed);
  std::uniform_int_distribution<int64_t> ds(0, total_comb - 1);
  for(int64_t i = 0; i < total_samples; i++) {
    int64_t x = ds(rnd);
    if (x >= base_index && x < base_index + index_step) {
      indexes.push_back(x);
    }
  }
  std::sort(indexes.begin(), indexes.end());
  indexes.erase(std::unique(indexes.begin(), indexes.end()), indexes.end());
}

void rnd_indexer_t::postprocess_rnd_container(std::vector<struct_info> &rnd_container) {
  if( rnd_container.size() <= random_count )
    return;
  rnd_engine_t rnd(new_seed(initial_seed, 1));
  std::shuffle(rnd_container.begin(), rnd_container.end(), rnd);
  rnd_container.resize(random_count);
  std::sort(rnd_container.begin(), rnd_container.end(),
            [](const struct_info a, const struct_info b) -> bool {
              return a.index < b.index;
            });
}
