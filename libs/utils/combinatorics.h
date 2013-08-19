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

inline
int64_t factorial(int k)
{
  int64_t result = 1;
  
  for(int i = 1; i <= k; i++)
    result *= i;
    
  return result;
}

inline
int64_t num_combinations(const std::vector<int> &nm)
{

  int64_t down_fact = 1;
  int64_t up_fact = 1;
  int counter = 1;
  
  bool limit_achived = false;
  
  for(std::vector<int>::const_iterator it = nm.begin(); it != nm.end(); it++)
  {
    for(int i = 1; i <= *it; i++)
    {
      if( up_fact > LLONG_MAX / counter )
      {
        limit_achived = true;
        break;
      }  
      
      up_fact *= counter;
      
      if(up_fact % i == 0)
        up_fact /= i;
      else          
        down_fact *= i;
      
      counter++;
    }
    if( limit_achived )
      break;
  }

  if(!limit_achived)
  {  
    assert(up_fact % down_fact == 0);
    return up_fact / down_fact;
  }
  else
    return LLONG_MAX;
}

template < typename T>
int64_t num_combinations(const std::map<T, int> &nm)
{
  typedef std::map<T, int> map_type;
  
  std::vector<int> nmv;
  for(typename map_type::const_iterator it = nm.begin(); it != nm.end(); it++)
    nmv.push_back(it->second);
  
  return num_combinations(nmv);  
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


#endif	/* COMBINATORICS_H */

