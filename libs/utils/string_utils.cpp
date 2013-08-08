/* 
 * File:   string_utils.cpp
 * Author: kokho
 * 
 * Created on May 2, 2012, 3:49 PM
 */

#include <vector>
#include <sstream>
#include <set>
#include <assert.h>
#include <map>

#include "string_utils.h"

std::string &trim(std::string &str)
{
  int i,j,start,end;

  //ltrim
  for (i = 0; (str[i] != 0 && str[i] <= 32); )
    i++;
  
  start=i;

  //rtrim
  for(i = 0, j = 0; str[i] != 0; i++)
    j = ((str[i]<=32)? j + 1 : 0);
  
  end = i - j;
  str = str.substr(start, end - start);
  
  return str;
}

std::string trim_delim(const std::string str, const std::string delimeter)
{
  std::string result = str;
  
  while (result.substr(0, delimeter.size()) == delimeter)
    result = result.substr(delimeter.size(), std::string::npos);
  
  while( (result.length() > 0) && (result.substr(result.length() - delimeter.size(), delimeter.size()) == delimeter) )
    result = result.substr(0, result.length() - delimeter.length());

  return result;
}


std::vector<std::string> split_string_vector(std::string str, const std::string delimeter)
{
  std::vector<std::string> result;
  
  while( str != "" )
  {
    str = trim_delim(str, delimeter);
    int pos = str.find(delimeter);
    result.push_back(str.substr(0, pos));
    if ( pos != std::string::npos )
      str = str.substr(pos + delimeter.size());
    else
      str = "";
  }
    
  return result;
}

std::set<std::string> split_string_set(std::string str, const std::string delimeter)
{
  std::set<std::string> result;
  result.clear();
  
  std::vector<std::string> sv = split_string_vector(str, delimeter);
  
  for(int i = 0; i < sv.size(); i++)
  {
    result.insert(sv[i]);
  }  
  
  return result;
}        

bool is_double(const std::string& s, double& r_double)
{
  std::istringstream i(s);
  
  if (i >> r_double)
  {
    return true;
  }
  r_double = 0.0;
  return false;
}

std::string wildcard_to_regex(const std::string &wildcard)
{
  using namespace std;

  std::string to_change = "+|()[]{}\\.^$";
  std::string result = "";
  
  for(int i = 0; i < wildcard.length(); i++)
  {
    char curr_char = wildcard[i];
    if( to_change.find(curr_char) != string::npos )
      result += "\\" + string(&curr_char);
    else if (curr_char == '*')
      result += ".*";
    else if (curr_char == '?')
      result += ".{1}";      
    else  
      result += curr_char;
  }  
            
  return result;          
}        

bool accept_pos(int pos, const std::string &cmp_str, const std::string &sub_str)
{
  bool result;
  
  result = (cmp_str.length() - pos) >= sub_str.length();
  if(!result)
    return result;
    
  result = true;
  for(int i = 0; i < sub_str.length(); i++)
  {
    result = (cmp_str[i + pos] == sub_str[i]) || (sub_str[i] == '?');
    if(!result)
      break;
  }  
  
  return result;  
}


bool match_wildcard(const std::string &wc_str, const std::string &cmp_str)
{
  std::vector<std::string> spv = split_string_vector(wc_str, "*");
  
  if(wc_str[0] == '*')
    spv.insert(spv.begin(), "");

  if(wc_str[wc_str.length() - 1] == '*')
    spv.insert(spv.end(), "");
  
  std::vector<int> pos;
  pos.resize(spv.size() + 1, 0);
  
  int curr_item = 0;
  do
  {
    int accepted_pos = -1;
    for(int i = pos[curr_item]; i < cmp_str.length(); i++)     
    {
      if(accept_pos(i, cmp_str, spv[curr_item]))
      {
        accepted_pos = i;
        break;
      }
    }
    
    if(accepted_pos >= 0)
    {  
      pos[curr_item] = accepted_pos + 1;
      pos[curr_item + 1] = accepted_pos + spv[curr_item].length();
      curr_item++;
    }  
    else
      curr_item--;
        
  } while ( (curr_item >= 0) && (curr_item < spv.size()) );
 
  return curr_item > 0;
}

std::string get_index_str(const int index, const int index_max, 
                          const std::string filler)
{
  assert(filler.length() == 1);
  
  std::stringstream ss;

  ss << index;
  
  std::string result = ss.str();
  
  ss.str("");
  ss << index_max;
  
  std::string ref_string = ss.str();
  
  while (result.length() < ref_string.length())
    result.insert(0, filler);
    
  return result;
}        
