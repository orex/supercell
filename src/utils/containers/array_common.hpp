/*
 * File:   array_common.hpp
 * Author: kirill
 *
 * Created on January 13, 2014, 3:13 PM
 */

#ifndef ARRAY_COMMON_HPP
#define	ARRAY_COMMON_HPP

#include <vector>
#include <algorithm>

namespace array_common
{
  template <typename T>
  void delete_singles(std::vector<T> &data)
  {
    std::sort(data.begin(), data.end());
    int pos_cmb = 0;
    for(int i = 0; i < int(data.size()) - 1; i++)
    {
      if( data[i] == data[i + 1])
      {
        data[pos_cmb] = data[i];
        pos_cmb++;
      }
    }
    data.resize(pos_cmb);
  }
}

#endif	/* ARRAY_COMMON_HPP */
