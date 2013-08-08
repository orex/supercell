/* 
 * File:   parse_d2o_input.cpp
 * Author: kirill
 * 
 * Created on July 19, 2013, 11:27 AM
 */

#include "parse_d2o_input.h"

#include "string_utils.h"

#include <assert.h>

using namespace std;
using namespace boost;

bool parse_d2o_input::get_supercell_size(const std::string &sc_size_str, std::vector<int> &scs)
{
  assert(scs.size() == 3);
  bool result;
  
  smatch what;
  const regex e("\\s*(\\d+)x(\\d+)x(\\d+)\\s*");
  result = regex_match(sc_size_str, what, e);
  
  if(result)
  {  
    scs[0] = lexical_cast<int>(what[1]);
    scs[1] = lexical_cast<int>(what[2]);
    scs[2] = lexical_cast<int>(what[3]);
  }
  
  return result;
}


bool parse_d2o_input::get_charge_balance(std::string cb_str, charge_balance &cb)
{
  bool result = false;
  trim(cb_str);
  
  result = cb_names::get_cb(cb_str, cb);

  return result;
}
