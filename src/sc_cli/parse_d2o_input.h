/* 
 * File:   parse_d2o_input.h
 * Author: kirill
 *
 * Created on July 19, 2013, 11:27 AM
 */

#ifndef PARSE_D2O_INPUT_H
#define	PARSE_D2O_INPUT_H

#include <string>
#include <vector>
#include <map>
#include <eigen3/Eigen/Core>

#include <regex>

#include "common_types.h"

class c_man_atom_prop_item_cli : public c_man_atom_prop_item
{
public:
  enum class lbl_type {ltPlain, ltWC, ltRegex};
  
  std::string label;
  lbl_type l_type;
};

class c_man_atom_prop_cli : public c_man_atom_prop
{
protected:
  int verbose_level;
  std::vector< c_man_atom_prop_item_cli > vc_raw;
public:
  void regex_test(std::string test_str);
  int  search_count(std::string &str, 
                    std::vector<std::string> &match,
                    const std::regex &rx);
  bool get_param(std::string &right_str, const std::regex &rx, int param_num, std::string &param);
  bool get_params(std::string right_str, c_man_atom_prop_item_cli &c_prop);
  
  bool get_labels(std::string left_str, 
                  std::vector< std::string > &lables_pattern,
                  c_man_atom_prop_item_cli::lbl_type &lt);
  
  bool parse_input_item(std::string inp_str);
public:
  void set_verbose(int level)
  { verbose_level = level; }
  bool parse_input(const std::vector<std::string> &inp, std::string &param_error);
  virtual void convert_properties(const std::set<std::string> &labels);
};

class parse_d2o_input 
{
public:
  static bool get_supercell_size(const std::string &sc_size_str, std::vector<int> &scs);
  static bool get_charge_balance(std::string cb_str, charge_balance &cb);
};

bool parse_sel_input(const std::vector<std::string> &inp, c_struct_sel & out, std::string &param_error);

#endif	/* PARSE_D2O_INPUT_H */

