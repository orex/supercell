/* 
 * File:   c_assigned_value_serialize.h
 * Author: kirill
 *
 * Created on September 17, 2013, 10:48 AM
 */

#ifndef C_ASSIGNED_VALUE_SERIALIZE_H
#define	C_ASSIGNED_VALUE_SERIALIZE_H

#include "av_base.hpp"
#include <boost/serialization/split_free.hpp>

namespace boost {
namespace serialization {

  template<class Archive, typename T>
  inline void serialize(Archive & ar, c_assigned_value<T> & g, const unsigned int version)
  {  
    boost::serialization::split_free(ar, g, version); 
  }
  
  template<class Archive, typename T>
  inline void save(Archive & ar, const c_assigned_value<T> & g, unsigned int version)
  {
    bool assigned = g.assigned();
    ar & assigned;
    if(assigned)
    {  
      T value = g.value();
      ar & value;
    }  
  }
  
  template<class Archive, typename T>
  inline void load(Archive & ar, c_assigned_value<T> & g, unsigned int version)
  {
    g.clear();
    bool assigned;
    ar & assigned;
    if(assigned)
    {  
      T value;
      ar & value;
      g = value;
    }  
  }  

} // namespace serialization
} // namespace boost

#endif	/* C_ASSIGNED_VALUE_SERIALIZE_H */

