/* 
 * File:   combinatorics.h
 * Author: kirill
 *
 * Created on July 17, 2013, 11:32 AM
 */

#ifndef COMBINATORICS_H
#define	COMBINATORICS_H

#include <vector>
#include <map>
#include <assert.h>
#include <climits>
#include <algorithm>

template <typename T>
int get_sign(const T &v) {
  if( v == 0)
  {
    return 0;
  }
  else if (v > 0)
  {
    return 1;
  }
  else
  {
    return -1;
  }

}

template <typename To, typename Ti>
bool safe_multiplication(const std::vector<Ti> &v, To &result)
{
  result = 0;
  if(v.empty())
    return false;

  int sign = get_sign(v[0]);
  for(int i = 1; i < v.size(); i++) {
    sign *= get_sign(v[i]);
  }

  if(sign == 0) {
    return true;
  }
  else if (sign < 0 && std::numeric_limits<To>::min() == 0)
  {
    return false;
  }

  result = 1;

  bool overflow = false;
  for(int i = 0; i < v.size(); i++)
  {
    To vm = static_cast<To>(std::abs(v[i]));
    if( vm <= std::numeric_limits<To>::max() / result )
    {
      result = result * vm;
    }
    else
    {
      overflow = true;
      result = 0;
      break;
    }
  }
  if(!overflow)
    result = result * static_cast<To>(sign);

  return !overflow;
}

template <typename T>
bool factorial(int k, T &result)
{
  std::vector<int> v;
  v.reserve(k);
  for(int i = 1; i <= k; i++)
    v.push_back(i);
    
  return safe_multiplication(v, result);
}

inline int __delete_vector_by_value(std::vector<int> &nm, int value, int max_operation = -1)
{
  assert(value > 0);
  int result = 0;
  if(max_operation == 0)
    return 0;

  for(int i = 0; i < nm.size(); i++)
  {
    while(nm[i] % value == 0)
    {
      nm[i] /= value;
      result++;
      if(result == max_operation)
        break;
    }
    if(result == max_operation)
      break;
  }
  return result;
}

template <typename T>
bool num_combinations(const std::vector<int> &nm, T &result)
{
  if(nm.empty()) {
    result = 0;
    return true;
  }
  std::vector<int> top_fract;
  std::vector<int> bottom_fract;

  int max_value = 0;
  int counter = 1;
  for(std::vector<int>::const_iterator it = nm.begin(); it != nm.end(); it++)
  {
    for(int i = 1; i <= *it; i++)
    {
      max_value = std::max(max_value, i);
      bottom_fract.push_back(i);
      top_fract.push_back(counter);
      counter++;
    }
  }

  if(max_value == 0) {
    result = 0;
    return true;
  }


  for(int i = max_value; i > 1; i--)
  {
    std::vector<int> bt1, top1;
    top1 = top_fract;
    bt1 = bottom_fract;
    int num_delete_bottom = __delete_vector_by_value(bt1, i);
    int num_delete_top = __delete_vector_by_value(top1, i);
    int num_delete = std::min(num_delete_bottom, num_delete_top);
    num_delete_bottom = __delete_vector_by_value(bottom_fract, i, num_delete);
    num_delete_top = __delete_vector_by_value(top_fract, i, num_delete);

    assert(num_delete_bottom == num_delete_top && num_delete_top == num_delete);
  }
  int bottom_result;
  bool sm_bottom = safe_multiplication(bottom_fract, bottom_result);
  assert(sm_bottom && (bottom_result == 1) );

  return safe_multiplication(top_fract, result);
}

template < typename Tm, typename To>
bool num_combinations(const std::map<Tm, int> &nm, To & result)
{
  typedef std::map<Tm, int> map_type;
  
  std::vector<int> nmv;
  for(typename map_type::const_iterator it = nm.begin(); it != nm.end(); it++)
    nmv.push_back(it->second);
  
  return num_combinations(nmv, result);
}

template < typename T>
std::vector<T> create_start_combination(const std::map<T, int> &nm)
{
  std::vector<T> result;
  for(typename std::map<T, int>::const_iterator it = nm.begin(); it != nm.end(); it++)
  {
    for(int j = 0; j < (*it).second; j++)
      result.push_back((*it).first);
  }
  
  bool prev_perm_exist = std::prev_permutation(result.begin(), result.end());
  bool next_perm_exist = std::next_permutation(result.begin(), result.end());
  
  assert( (!prev_perm_exist) && (!next_perm_exist) );
  
  return result;
}

template <class _Compare, class _ConstBidirectionalIterator, typename IndexType>
inline IndexType __permutation_count(_ConstBidirectionalIterator __first,
                              _ConstBidirectionalIterator __last, _Compare __comp) {
  static_assert(std::is_arithmetic<IndexType>::value, "Result type should be arithmetic.");
  _ConstBidirectionalIterator __it = __last;
  if( __first == __last )
    return IndexType(0);

  auto eq = [&__comp](const typename _ConstBidirectionalIterator::value_type &__a,
      const typename _ConstBidirectionalIterator::value_type &__b) -> bool {
    return !__comp(__a, __b) && !__comp(__b, __a);
  };

  std::vector<IndexType> __mx;
  __mx.reserve(16);

  IndexType result(1);
  IndexType len(1);
  for(auto it = __first; it != __last; ++it) {
    auto __itm = std::find_if(__mx.crbegin(), __mx.crend(), eq);
    IndexType l = (__itm == __mx.crend()) ? IndexType(0) : (*__itm + 1);
    IndexType ld = result / l, lr = result % l;
    result = ld * len + (lr * len) / l;
    len++;
  }

  return result;
}

template <class _Compare, class _BidirectionalIterator>
inline
_BidirectionalIterator
__next_permutation_it(_BidirectionalIterator __first, _BidirectionalIterator __last, _Compare __comp)
{
  _BidirectionalIterator __i = __last;
  if (__first == __last || __first == --__i)
    return __last;
  while (true)
  {
    _BidirectionalIterator __ip1 = __i;
    if (__comp(*--__i, *__ip1))
    {
      _BidirectionalIterator __j = __last;
      while (!__comp(*__i, *--__j))
        ;
      std::swap(*__i, *__j);
      std::reverse(__ip1, __last);
      return __i;
    }
    if (__i == __first)
    {
      std::reverse(__first, __last);
      return __last;
    }
  }
}

template <class _BidirectionalIterator, class _Compare>
inline
_BidirectionalIterator
next_permutation_it(_BidirectionalIterator __first, _BidirectionalIterator __last, _Compare __comp)
{
  return __next_permutation_it(__first, __last, __comp);
}

template <class _BidirectionalIterator>
inline
_BidirectionalIterator
next_permutation_it(_BidirectionalIterator __first, _BidirectionalIterator __last)
{
  return __next_permutation_it(__first, __last,
                          std::less<typename std::iterator_traits<_BidirectionalIterator>::value_type>());
}

template <class _Compare, class _ConstBidirectionalIterator, typename IndexType>
inline
IndexType __permutation_index(_ConstBidirectionalIterator __first, _ConstBidirectionalIterator __last, _Compare __comp) {
  static_assert(std::is_arithmetic<IndexType>::value, "Result type should be arithmetic.");
  IndexType result(0);
  _ConstBidirectionalIterator __it = __last;
  if( __first == __last )
    return result;

  typedef std::pair<_ConstBidirectionalIterator, IndexType> pr_t;
  std::vector<pr_t> __mx;
  __mx.reserve(16);

  auto cmp = [&__comp](const pr_t &__a, const typename _ConstBidirectionalIterator::value_type &__b) -> bool {
    return __comp(*(__a.first), __b);
  };
  int len = 0;
  IndexType tot_prm = 1;
  do {
    --__it;
    auto __itmx = std::lower_bound(__mx.begin(), __mx.end(), *(__it), cmp);
    if( __itmx == __mx.end() || __comp(*(__it), *(__itmx->first)) ) {
      __itmx = __mx.emplace(__itmx, __it, 0);
    }
    __itmx->second++;
    IndexType ms(0);
    for(auto __itx = __mx.begin(); __itx != __itmx; ++__itx) {
      ms += __itx->second;
    }
    result += (tot_prm / __itmx->second) * ms + ((tot_prm % __itmx->second) * ms) / __itmx->second;
    len++;
    tot_prm = (tot_prm / __itmx->second) * len + ((tot_prm % __itmx->second) * len) / __itmx->second;
  } while( __it != __first );
  return result;
}


template < typename IndexType, class _ConstBidirectionalIterator, class _Compare>
inline
IndexType
permutation_index(_ConstBidirectionalIterator __first, _ConstBidirectionalIterator __last, _Compare __comp)
{
  return __permutation_index<_Compare, _ConstBidirectionalIterator, IndexType>(__first, __last, __comp);
}

template <typename IndexType, class _ConstBidirectionalIterator>
inline
IndexType
permutation_index(_ConstBidirectionalIterator __first, _ConstBidirectionalIterator __last)
{
  typedef std::less<typename std::iterator_traits<_ConstBidirectionalIterator>::value_type> less_t;
  return __permutation_index<less_t, _ConstBidirectionalIterator, IndexType>(__first, __last, less_t());
}

template <typename IndexType, class _BidirectionalIterator, class _Compare>
inline IndexType next_k_permutations(_BidirectionalIterator __first,
                                     _BidirectionalIterator __last, IndexType k,
                                     _Compare __comp) {
  if (__first == __last)
    return k;

  _BidirectionalIterator p_last = __last;
  p_last--;
  if (__first == p_last)
    return k;

  IndexType result = 0;
  while (k > 0) {
    auto itc = p_last;
    auto itp = itc;
    for (itp = itc; itc != __first; --itc) {
      itp--;
      if (__comp(*itp, *itc) || __comp(*itc, *itp))
        break;
    }
    if (itc == __first)
      return k;

    bool acc = __comp(*itp, *itc);
    if (acc) {
      IndexType total_count = std::distance(itc, __last);
      IndexType total_cur_count = 0;
      IndexType km = 1;
      for (itp = itc; itc != __first; --itc) {
        itp--;
        if (__comp(*itc, *itp))
          break;

        total_count++;
        if (!__comp(*itp, *itc)) {
          total_cur_count++;
        } else {
          total_cur_count = 1;
        }
        IndexType kt = (km * total_count) / total_cur_count;
        if (kt > k)
          break;
        km = kt;
      }
      // All permutations run from fist permutation until last.
      // km = total_number of permutations for the sequence.
      if (itc == __first) {
        result += k / km;
        k %= km;
      } else if (km == 1) {
        k--;
        if (!std::next_permutation(__first, __last, __comp))
          result++;
      } else {
        k -= km - 1;
        std::reverse(itc, __last);
      }
    } else {
      k--;
      if (!std::next_permutation(__first, __last, __comp))
        result++;
    }
  }
  return result;
}

template <typename IndexType, class _BidirectionalIterator>
inline IndexType next_k_permutations(_BidirectionalIterator __first,
                                     _BidirectionalIterator __last,
                                     IndexType k) {
  typedef std::less<
      typename std::iterator_traits<_BidirectionalIterator>::value_type>
      less_t;
  return next_k_permutations(__first, __last, k, less_t());
}

#endif	/* COMBINATORICS_H */

