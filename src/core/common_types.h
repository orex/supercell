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

#include <boost/optional.hpp>

enum charge_balance {cb_no, cb_yes, cb_try};

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
  boost::optional<double> charge;
  boost::optional<bool> fixed;
  boost::optional<int> population;
 
  c_man_atom_prop_item() {};
  void assign(const c_man_atom_prop_item &orig)
  { 
    charge = orig.charge;
    fixed = orig.fixed;
    population= orig.population;
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

class c_struct_sel
{
protected:
  typedef std::map<std::string, int> sel_prop;
  sel_prop sp;
  int str_count(const std::string &st) const
  {
    if(sp.count(st) > 0)
      return sp.at(st);
    return 0;        
  };
public:
  bool save_all() const 
  { return sp.empty(); };
  int str_first_count() const 
  { return str_count("f"); };
  int str_last_count() const 
  { return str_count("a"); };
  int str_low_count() const 
  { return str_count("l"); };
  int str_high_count() const 
  { return str_count("h"); };
  int str_random_count() const 
  { return str_count("r"); };
  void assign_base(const c_struct_sel &orig)
  { sp = orig.sp; };
};        


#endif	/* COMMON_TYPES_H */

