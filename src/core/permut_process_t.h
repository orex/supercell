//
// Created by kirill on 4/26/21.
//

#ifndef SUPERCELL_SRC_CORE_PERMUT_PROCESS_T_H_
#define SUPERCELL_SRC_CORE_PERMUT_PROCESS_T_H_

#include <vector>
#include <boost/optional.hpp>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "containers/hash_unique.h"

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
 public:
  int ps_size;
  std::vector<struct_info> ps;

  permut_process_t(const t_symm_set &symmetries, const std::vector<int> &perm_indexes,
                   int max_proc_struct) :
      syms(symmetries), permi(perm_indexes),
      allcmb(symmetries.rows(), symmetries.cols()),
      ps_size(0), ps(max_proc_struct) {
    assert(perm_indexes.back() == symmetries.cols());
  };
  inline void set_proc_range(const t_vec_comb &first, const t_vec_comb &last) {
    first_cmb = first;
    last_cmb = last;
  }

  void process_merge();
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

template<class Container, typename Index>
Index next_k_complex_permutation(Container &vc, Container &vc_last, Index k, const std::vector<int> &permi) {
  for(Index i = 0; i < k; i++) {
    vc_last = vc;
    if( !next_complex_permutation(vc, permi) )
      return i + 1;
  }
  return k;
}



#endif //SUPERCELL_SRC_CORE_PERMUT_PROCESS_T_H_
