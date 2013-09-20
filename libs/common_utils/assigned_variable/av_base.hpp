/* 
 * File:   TAssignedValue.h
 * Author: kirill
 *
 * Created on August 14, 2013, 10:55 AM
 */

#ifndef TASSIGNEDVALUE_H
#define	TASSIGNEDVALUE_H

#include <string>
#include <cassert>
#include <ostream>

#include <boost/lexical_cast.hpp>

template <typename T>
class c_assigned_value;

template <typename T>
std::ostream& operator<<(std::ostream& os, const c_assigned_value<T> &cav);

template <typename T>
class c_assigned_value 
{
protected:  
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

  T value_def(const T def_value) const
  {
    if( assigned() )
      return value(); 
    else
      return def_value;
  }  
  
  std::string to_string() const
  {
    if( assigned() )
      return boost::lexical_cast<std::string>(_value);
    else
      return "N/A";
  }
  
  void assign(const c_assigned_value<T> &p)
  {
    if( p.assigned() )
      *this = p;
  }
  
  void clear()
  {
    _assigned = false;
  }
  
  friend std::ostream& operator<< <T>(std::ostream& os, const c_assigned_value<T> &cav);
  
  virtual ~c_assigned_value() {};
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const c_assigned_value<T> &cav)
{
  if(cav.assigned())
    os << cav.value();
  else
    os << "N/A";
    
  return os;
}

#endif	/* TASSIGNEDVALUE_H */

