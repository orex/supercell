/* 
 * File:   TAssignedValue.h
 * Author: kirill
 *
 * Created on August 14, 2013, 10:55 AM
 */

#ifndef TASSIGNEDVALUE_H
#define	TASSIGNEDVALUE_H

#include <string>
#include <assert.h>

#include <boost/lexical_cast.hpp>

template <typename T>
class c_assigned_value 
{
  bool _assigned;
  T    _value;
public:
  c_assigned_value(): _assigned(false) {};
  
  c_assigned_value<T> & operator=(const T &p)
  { 
    _assigned = true;
    _value = p;
    return *this;
  };
  
  bool assigned() const
  {return _assigned; }
  
  T value() const
  {
    assert(_assigned);
    return _value; 
  }  
  
  std::string to_string() const
  {
    if( assigned() )
      return boost::lexical_cast<std::string>(_value);
    else
      return "N/A";
  }
  
  bool assign(const c_assigned_value<T> &p)
  {
    if( p.assigned() )
      *this = p;
  }
  
  virtual ~c_assigned_value() {};
};

#endif	/* TASSIGNEDVALUE_H */

