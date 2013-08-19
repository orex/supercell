/* 
 * File:   parse_d2o_input.cpp
 * Author: kirill
 * 
 * Created on July 19, 2013, 11:27 AM
 */

#include "parse_d2o_input.h"

#include "string_utils.h"
#include "d2o_main_class.h"

#include <assert.h>

#include <boost/regex.hpp>

#include <iostream>

using namespace std;
using namespace boost;  

int c_man_atom_prop_cli::search_count(std::string &str, 
                                      std::vector<std::string> &match,
                                      const boost::regex &rx)
{
  int result = 0;
  
  smatch r_match;      

  while(regex_search(str, r_match, rx))
  {
    string replaced_str = "";
    
    match.clear();    
    for(smatch::const_iterator it = r_match.begin(); it != r_match.end(); it++)
      match.push_back(*it);
    
    replaced_str = str.substr(0, r_match.position()) + " " + 
                   str.substr(r_match.position() + r_match.length());
    
    str = replaced_str;    
    
    result++;
  }
  
  return result;
}

void c_man_atom_prop_cli::regex_test(std::string test_str)
{
  regex rx_charge("(^|\\s)(fixed|notfixed)(\\s|$)");  
  vector<string> rx_out;
  int find = search_count(test_str, rx_out, rx_charge);
  cout << "find: " << find << endl;
  if(find == 1)
  {
    cout << rx_out[0] << endl;
    cout << rx_out[1] << endl;    
    cout << rx_out[2] << endl;    
    cout << rx_out[3] << endl;
  }  
}

bool c_man_atom_prop_cli::get_param(std::string &right_str, const boost::regex &rx, 
                                    int param_num, std::string &param)
{
  bool result;
  
  vector<string> match;
  int sc = search_count(right_str, match, rx);
  if(sc == 0)
  {
    param = "";
    result = true;
  }
  else if( sc == 1 )
  {
    param = match[param_num];
    trim(param);
    result = param != "";
  } 
  else if( sc > 1 )
  {
    result = false;
  }  
  
  return result;
}

bool c_man_atom_prop_cli::get_params(std::string right_str, c_man_atom_prop_item_cli &c_prop)
{
  bool result = true;
  string param;
  bool good_param;
  
  trim(right_str);
  
  if(right_str.length() < 2)
    return false;
  
  if( (right_str[0] == '{' ) && (right_str[right_str.length() - 1] == '}' ) )
    right_str = right_str.substr(1, right_str.length() - 2);

  trim(right_str);
  if(right_str.length() < 2)
    return false;
  
  //extract charge
  regex rx_charge("(^|\\s)c(harge){0,1}\\s*=\\s*(\\S*)(\\s|$)");
  good_param = get_param(right_str, rx_charge, 3, param);
  result = result && good_param;
  if(good_param && (param != "") )
  {
    double tmp;
    if( is_double(param, tmp) )
      c_prop.charge = tmp;
    else
      result = false;
  }  
  
  //extract population
  regex rx_pop("(^|\\s)p(opulation){0,1}\\s*=\\s*(\\d+)(\\s|$)");
  good_param = get_param(right_str, rx_pop, 3, param);
  result = result && good_param;
  if(good_param && (param != "") )
    c_prop.population = lexical_cast<int>(param);


  //extract fixation
  regex rx_fix("(^|\\s)(fixed|notfixed)(\\s|$)");
  good_param = get_param(right_str, rx_fix, 2, param);
  result = result && good_param;
  if(good_param && (param != "") )
    c_prop.fixed = (param == "fixed");
  
  trim(right_str);
  
  result = result && ( (right_str == "" ) || (regex_match(right_str, regex("\\s*\\{\\s*\\}\\s*"))) );
   
  return result;  
}

bool c_man_atom_prop_cli::get_labels(std::string left_str, 
                                     std::vector< std::string > &lables_pattern,
                                     c_man_atom_prop_item_cli::lbl_type &lt)
{
  regex left_rx("(^[wWrRpP]{0,1})\\((.+)\\)$");
  smatch rx_out;
 
  string str_lbl;
  if(regex_match(left_str, rx_out, left_rx))
  {
    str_lbl = rx_out[2];
    lt = c_man_atom_prop_item_cli::ltPlain;
    if( (rx_out[1] == "W") || (rx_out[1] == "w") )
      lt = c_man_atom_prop_item_cli::ltWC;
    
    if( (rx_out[1] == "R") || (rx_out[1] == "r") )
      lt = c_man_atom_prop_item_cli::ltRegex;
  }
  else
  {
    str_lbl = left_str;
    lt = c_man_atom_prop_item_cli::ltWC;
  }  
  
  lables_pattern.clear();
  lables_pattern = split_string_vector(str_lbl, " ");
  
  bool result = lables_pattern.size() > 0;
  
  if( lt == c_man_atom_prop_item_cli::ltRegex )
  {  
    for(int i = 0; i < lables_pattern.size(); i++)
    { 
      try
      { regex rx(lables_pattern[i]); }        
      catch (boost::regex_error& e)
      {  
        cerr << "Regex " << lables_pattern[i] << " is wrong." << endl;
        cerr << e.what() << endl;
        result = false; 
      }

      if( !result )
        break;
    }
  }

  cout << result << endl;  
  
  if( !result )
    lables_pattern.clear();
    
  return result;
}


bool c_man_atom_prop_cli::parse_input_item(std::string inp_str)
{
  bool result = true;

  int delim_pos = -1;
  int delim_count = 0;
  
  for(int i = 0; i < inp_str.length(); i++)
  {
    if( inp_str[i] == ':' )
    {
      bool double_colon = false;
      if( i + 1 < inp_str.length() )
      {  
        if( inp_str[i + 1] == ':' )
        {  
          inp_str.erase(i + 1, 1);
          double_colon = true;
        }  
      }  
      
      if(!double_colon )
      {  
        delim_pos = i;
        delim_count++;
      }  
    }  
  }  
 
  if(verbose_level >= 4)
  {  
    cout << inp_str << endl;
    cout << "delim_pos " << delim_pos << endl;
    cout << "delim_count " << delim_count << endl;
  }
  
  if((delim_count != 1) || (delim_pos == 0) || (delim_pos == inp_str.length() - 1) )
    return false;
  
  string left_str  = inp_str.substr(0, delim_pos);
  string right_str = inp_str.substr(delim_pos + 1);
  
  trim(left_str);
  trim(right_str);
  
  if( (left_str == "") || (right_str == "") )
    return false;
  
  vector<string> lbls;
  
  c_man_atom_prop_item_cli::lbl_type lblt;
  c_man_atom_prop_item_cli prop_item;
  
  if( !get_labels(left_str, lbls, lblt) )
    return false;
  
  if( !get_params(right_str, prop_item) )
    return false;
  
  for(size_t i = 0; i < lbls.size(); i++)
  {
    prop_item.label = lbls[i];
    prop_item.l_type = lblt;
    vc_raw.push_back(prop_item);
  }  
    
  return result;
}

bool c_man_atom_prop_cli::parse_input(const std::vector<std::string> &inp, std::string &param_error)
{
  bool result;
  
  vc_raw.clear();
  
  param_error = "";
  
  for(size_t i = 0; i < inp.size(); i++)
  {  
    result = parse_input_item(inp[i]);
    if(!result) 
    {  
      param_error = inp[i];
      break;
    }  
  }  

  return result;
}

void c_man_atom_prop_cli::convert_properties(const std::set<std::string> &labels)
{
  for(set<string>::const_iterator it_lbl =  labels.begin(); 
                                  it_lbl != labels.end(); it_lbl++)
  {
    for(vector< c_man_atom_prop_item_cli >::iterator it_p  = vc_raw.begin();
                                                     it_p != vc_raw.end(); it_p++)
    {
      bool the_label = false;
      the_label = the_label  || ( (it_p->l_type == c_man_atom_prop_item_cli::ltPlain) && 
                                  (*it_lbl == it_p->label) ); 
      the_label = the_label  || ( (it_p->l_type == c_man_atom_prop_item_cli::ltWC) && 
                                  (match_wildcard(it_p->label, *it_lbl)) ); 
      the_label = the_label  || ( (it_p->l_type == c_man_atom_prop_item_cli::ltRegex) && 
                                  ( regex_match(*it_lbl, regex(it_p->label)) ) );
      if( the_label )
        data_map[*it_lbl].assign(*it_p);
    }  
  }
}



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
