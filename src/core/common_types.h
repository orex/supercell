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
#include <vector>
#include <set>

#include "assigned_variable/av_base.hpp"

enum charge_balance {cb_no, cb_input, cb_try};

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
  static bool get_cb(const std::string &str, charge_balance &cb);
};

class c_man_atom_prop_item
{
public:
  c_assigned_value<double> charge;
  c_assigned_value<bool> fixed;  
  c_assigned_value<int> population;
 
  c_man_atom_prop_item() {};
  void assign(const c_man_atom_prop_item &orig)
  { 
    charge.assign(orig.charge);
    fixed.assign(orig.fixed);
    population.assign(orig.population);
  };
};        

class c_man_atom_prop
{
public:
  typedef std::map<std::string, c_man_atom_prop_item> data_type;
protected:
  data_type data_map;
public:
  c_man_atom_prop();
  const data_type &data()
  { return data_map; };
  c_man_atom_prop_item operator[](const std::string &lbl)
  { return data_map[lbl]; };
  bool exist(const std::string &lbl)
  { return data_map.count(lbl) != 0; };
  virtual void convert_properties(const std::set<std::string> &labels) = 0;
  ~c_man_atom_prop();
};        

#endif	/* COMMON_TYPES_H */

