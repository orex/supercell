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

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

#include "common_types.h"

class parse_d2o_input 
{
public:
  static bool get_supercell_size(const std::string &sc_size_str, std::vector<int> &scs);
  static bool get_charge_balance(std::string cb_str, charge_balance &cb);

  template <typename T>
  static bool get_manual_populations_map(const std::vector<std::string> &mp_str, 
                                         std::map<std::string,T> &mp_map, 
                                         const std::string data_rx = "[-+]?[0-9]*\\.?[0-9]")
  { 
    bool result = true;
    mp_map.clear();
    boost::smatch what;
    const boost::regex e("\\s*(.+):(" + data_rx + ")\\s*");

    for(int i = 0; i < mp_str.size(); i++)
    {
      result = result && boost::regex_match(mp_str[i], what, e);
      if(!result) break;
      mp_map[what[1]] = boost::lexical_cast<T>(what[2]);
    }  
  
    if(!result)  
      mp_map.clear();
    
    return result;
  }  
};

#endif	/* PARSE_D2O_INPUT_H */

