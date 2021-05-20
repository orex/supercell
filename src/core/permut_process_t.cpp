//
// Created by kirill on 4/26/21.
//

#include "permut_process_t.h"
#include "science/combinatorics.h"

#define XXH_INLINE_ALL
#define XXH_STATIC_LINKING_ONLY
#include <xxhash.h>

t_vec_comb permut_process_t::setmm_values(const t_vec_comb &vc, int start_index, bool max_value) const {
  auto it = std::upper_bound(permi.cbegin(), permi.cend(), start_index);
  if( it == permi.cend() ) {
    return vc;
  }
  int pmipos = std::distance(permi.cbegin(), it);
  int sti = start_index;
  t_vec_comb result = vc;
  int mv = result[start_index];
  for (int i = sti, j = pmipos; i < vc.size(); i++) {
    if (i == permi[j]) {
      std::fill(result.begin() + sti, result.begin() + i, mv);
      mv = result[i];
      j++;
      sti = i;
    }
    if( max_value == (result[i] > mv) )
      mv = result[i];
  }
  std::fill(result.begin() + sti, result.end(), mv);
  return result;
}

bool permut_process_t::next_comb(t_vec_comb &vc) const {
  return next_complex_permutation(vc, permi);
}

void permut_process_t::process_no_merge() {
  ps_size = 0;
  auto vc = first_cmb;
  while(true) {
    ps[ps_size].weight = 1;
    ps[ps_size].cmb = vc;
    ps_size++;
    assert(ps_size <= ps.size());
    if( vc != last_cmb )
      next_comb(vc);
    else
      break;
  }
}

void permut_process_t::process_merge() {
  ps_size = 0;
  int blength = first_cmb.size() * sizeof(t_vec_comb::value_type);
  int syms_num = syms.rows();
  int start_index = 0;

  for(start_index = 0; start_index < first_cmb.size() && first_cmb[start_index] == last_cmb[start_index] ; start_index++) {};

  auto fmcomp = [blength](const auto & first1, const auto &first2) -> int {
    return memcmp(&(*first1), &(*first2), blength);
  };
  auto hash_comb = [blength](const t_vec_comb::value_type * p) -> std::size_t {
    return XXH3_64bits(p, blength);
  };
  auto eq_comb = [blength](const t_vec_comb::value_type * p1, const t_vec_comb::value_type * p2) -> bool {
    return p1 == p2 || memcmp(p1, p2, blength) == 0;
  };
  bool has_lower_syms = false;
  t_vec_comb vmax = setmm_values(first_cmb, start_index, true);
  for(int i = 0; i < syms_num && !has_lower_syms; i++) {
    auto itc = allcmb.get_iterator(i, 0);
    auto its = syms.get_iterator(i, 0);
    for(int j = 0; j < vmax.size(); j++) {
      *(itc + *(its + j)) = vmax[j];
    }
    has_lower_syms = fmcomp(itc, first_cmb.cbegin()) < 0;
  }
  if( has_lower_syms ) {
    // std::cout << "Lower syms " << std::endl;
    return;
  }
  t_vec_comb vmin = setmm_values(first_cmb, start_index, false);
  std::vector<bool> always_higher(syms_num);
  for(int i = 0; i < syms_num; i++) {
    auto itc = allcmb.get_iterator(i, 0);
    auto its = syms.get_iterator(i, 0);
    for(int j = start_index; j < vmin.size(); j++) {
      *(itc + *(its + j)) = vmin[j];
    }
    always_higher[i] = fmcomp(itc, last_cmb.cbegin()) > 0;
  }

  hash_set<t_vec_comb::value_type *, decltype(hash_comb), decltype(eq_comb)> hs(syms_num, hash_comb, eq_comb);

  auto vc = first_cmb;

  while(true) {
    has_lower_syms = false;
    for(int i = 0; i < syms_num && !has_lower_syms; i++) {
      if( always_higher[i] )
        continue;
      auto itc = allcmb.get_iterator(i, 0);
      auto its = syms.get_iterator(i, 0);
      for(int j = start_index; j < vc.size(); j++) {
        *(itc + *(its + j)) = vc[j];
      }
      has_lower_syms = fmcomp(itc, vc.cbegin()) < 0;
    }
    // std::cout << (!has_lower_syms ? "*" : "" )  << std::endl;
    if( !has_lower_syms ) {
      hs.clear();
      for(int i = 0; i < syms_num; i++) {
        auto itc = allcmb.get_iterator(i, 0);
        if( always_higher[i] ) {
          auto its = syms.get_iterator(i, 0);
          for (int j = start_index; j < vc.size(); j++) {
            *(itc + *(its + j)) = vc[j];
          }
        }
        hs.put(&(*itc));
      }
      ps[ps_size].weight = hs.size();
      ps[ps_size].cmb = vc;
      ps_size++;
      assert(ps_size <= ps.size());
    }
    if( vc != last_cmb )
      next_comb(vc);
    else
      break;
  }
}

