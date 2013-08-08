/* 
 * File:   common_types.h
 * Author: kirill
 *
 * Created on July 29, 2013, 9:56 AM
 */

#ifndef COMMON_TYPES_H
#define	COMMON_TYPES_H

#include <string>
#include <map>

enum charge_balance {cb_no, cb_formal, cb_input, cb_try};

class cb_names
{
private:
  class cb_init
  {
  public:
    cb_init();
  };  
  
  friend class cb_init;
private:
  static cb_init _init;  
  static std::map<charge_balance, std::string> mp_cb;
public:  
  static std::string get_name(const charge_balance cb);
  static bool get_cb(const std::string str, charge_balance &cb);
};

#endif	/* COMMON_TYPES_H */

