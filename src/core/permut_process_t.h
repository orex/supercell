//
// Created by kirill on 4/26/21.
//

#ifndef SUPERCELL_SRC_CORE_PERMUT_PROCESS_T_H_
#define SUPERCELL_SRC_CORE_PERMUT_PROCESS_T_H_

#include <vector>
#include <boost/optional.hpp>
#include <Eigen/Core>
#include <Eigen/Geometry>

#define XXH_INLINE_ALL
#define XXH_STATIC_LINKING_ONLY
#define XXH_NO_STREAM
#define XXH_NO_STDLIB
#include <xxhash.h>

#include "containers/hash_unique.h"
#include "science/combinatorics.h"

#include "common_types.h"

template <typename T, std::size_t Aligment = sizeof(T)>
class t_proc_prm {
 private:
  typedef std::vector<T, boost::alignment::aligned_allocator<T, Aligment>> arr2d_t;
  int _cols;
  int _rows;
  int row_size;
  arr2d_t arr2d;
 public:
  typedef typename arr2d_t::const_iterator const_iterator;
  typedef typename arr2d_t::iterator iterator;

  t_proc_prm() : _cols(0), _rows(0), row_size(0) {};
  t_proc_prm(int rows, int cols, T def = T()) : _cols(cols), _rows(rows) {
    for(row_size = cols; row_size * sizeof(T) % Aligment != 0; row_size++) {};
    arr2d.resize(row_size * _rows, def);
  }

  inline int cols() const {
    return _cols;
  }

  inline int rows() const {
    return _rows;
  }

  inline T & at(int r, int c) {
    return arr2d.at(r * row_size + c);
  }

  inline const T & at(int r, int c) const {
    return arr2d.at(r * row_size + c);
  }

  inline T & operator()(int r, int c) {
    return arr2d[r * row_size + c];
  }

  inline const T & operator()(int r, int c) const {
    return arr2d[r * row_size + c];
  }

  const_iterator get_iterator(int r, int c) const {
    return arr2d.cbegin() + r * row_size + c;
  }

  iterator get_iterator(int r, int c) {
    return arr2d.begin() + r * row_size + c;
  }

  const_iterator next_row(const_iterator it) const {
    return it + row_size;
  }

  iterator next_row(iterator it) {
    return it + row_size;
  }

};

typedef t_proc_prm<base_prm_t, 32> t_all_comb;
typedef t_proc_prm<int, 32> t_symm_set;

class permut_process_t {
 private:
  const t_symm_set &syms;
  const std::vector<int> &permi;
  t_all_comb allcmb;
  t_vec_comb first_cmb;
  t_vec_comb last_cmb;
  t_vec_comb setmm_values(const t_vec_comb &vc, int start_index, bool max_value) const;
  bool next_comb(t_vec_comb &vc) const;
 private:
   class hash_comb {
   private:
     const int b_length;
   public:
     hash_comb(int cmb_size) : b_length(cmb_size * sizeof(t_vec_comb::value_type)) {};
     std::size_t operator()(const t_vec_comb::value_type * p) const {
       return XXH3_64bits(p, b_length);
     }
   };
   class eq_comb {
   private:
     const int b_length;
   public:
     eq_comb(int cmb_size) : b_length(cmb_size * sizeof(t_vec_comb::value_type)) {};
     inline bool operator()(const t_vec_comb::value_type * p1, const t_vec_comb::value_type * p2) const {
       return p1 == p2 || std::memcmp(p1, p2, b_length) == 0;
     };
   };
   hash_set<t_vec_comb::value_type *, hash_comb, eq_comb> hs;
 public:
  int ps_size;
  std::vector<struct_info_t> ps;
  std::int64_t total_combination_chunk;
  std::pair<double, double> minmax_energy;
  int min_weight;

  permut_process_t(const t_symm_set &symmetries, const std::vector<int> &perm_indexes,
                   int max_proc_struct) :
      syms(symmetries), permi(perm_indexes),
      allcmb(symmetries.rows(), symmetries.cols()),
      ps_size(0), ps(max_proc_struct), min_weight(0),
      hs(symmetries.rows(), hash_comb(symmetries.cols()), eq_comb(symmetries.cols()))
  {
    assert(perm_indexes.back() == symmetries.cols());
  };
  inline void set_proc_range(const t_vec_comb &first, const t_vec_comb &last) {
    first_cmb = first;
    last_cmb = last;
  }

  void process_merge() noexcept;
  void process_no_merge();
};

template<class Container>
bool next_complex_permutation(Container &vc, const std::vector<int> &permi) {
  for (auto it = --permi.cend(); it != permi.cbegin(); --it) {
    if (std::next_permutation(vc.begin() + *(it - 1), vc.begin() + *it)) {
      return true;
    }
  }
  return false;
}

template <class Container, typename Index>
Index next_k_complex_permutation(Container &vc, Container &vc_last, Index k,
                                 const std::vector<int> &permi) {
  if (k == 0)
    return k;

  vc_last = vc;
  Index km = k - 1;
  for (auto it = --permi.cend(); it != permi.cbegin() && km != 0; --it) {
    km = next_k_permutations(vc.begin() + *(it - 1), vc.begin() + *it, km);
  }

  if (km == 0) {
    vc_last = vc;
    next_complex_permutation(vc, permi);
    return k;
  } else {
    bool has_next_permutation = true;
    vc = vc_last;
    Index i;
    for (i = 0; i < k && has_next_permutation; i++) {
      vc_last = vc;
      has_next_permutation = next_complex_permutation(vc, permi);
    }
    return i;
  }
}

#endif //SUPERCELL_SRC_CORE_PERMUT_PROCESS_T_H_
