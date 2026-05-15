//
// Created by kirill on 4/18/21.
//

#ifndef HASH_UNIQUE_H_
#define HASH_UNIQUE_H_

#include <vector>
#include <limits>

template<class T, class Hash = std::hash<T>, class Eq = std::equal_to<T>, typename hashType = std::size_t >
class hash_set {
 private:
  static constexpr hashType c_ocp_flg = std::numeric_limits<hashType>::max() - (std::numeric_limits<hashType>::max() >> 1);
  typedef std::pair<hashType, T> bucket_t;
  typedef std::vector<bucket_t> buckets_t;
  buckets_t buckets;
  const Hash _h;
  const Eq _eq;
  std::size_t mask;
  std::size_t _size;
 public:
  typedef T value_type;
  hash_set(int num_elements, const Hash hash, const Eq eq ): _h(hash), _eq(eq), _size(0) {
    int cnt;
    for(cnt = std::max(16, num_elements * 4); (cnt & (cnt - 1)) != 0; cnt = (cnt & (cnt - 1))) {};
    buckets.resize(cnt, {0, T()});
    mask = cnt - 1;
  };

  bool put(const T & t) {
    hashType hs = _h(t) | c_ocp_flg;
    std::size_t i = hs & mask;
    for(i = i; true ; i = (i + 1) & mask ) {
      auto &bi = buckets[i];
      if (( bi.first & c_ocp_flg) == 0 ) {
        bi = {hs, t};
        _size++;
        return true;
      } else if( buckets[i].first == hs && _eq(t, bi.second) )
        return false;
    }
  };

  std::size_t size() const {
    return _size;
  }

  void clear() {
    _size = 0;
    std::fill(buckets.begin(), buckets.end(), bucket_t(0, T()));
  }
};



#endif //HASH_UNIQUE_H_
