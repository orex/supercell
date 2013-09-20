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

#include <iostream>

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
  
  result = sub_str == "";
  if(result)
    return result;
  
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

class wc_token
{
protected:
  const std::string &match_string;
  int pos_begin;
  bool reseted;
  virtual bool fit_next_intrln() = 0;  
public:
  wc_token(const std::string & match_string_v) : 
           match_string(match_string_v), reseted(false) {};
  virtual void reset(int pos_begin_v)
  { 
    pos_begin = pos_begin_v;
    reseted = true;
  };
  virtual int next_token_pos() = 0;
  bool fit_next()
  {
    if( reseted ) 
      return fit_next_intrln();
    else
      return false;
  }
};

class wc_star_token : public wc_token
{
protected:
  int size;    
  int max_size;

  virtual bool fit_next_intrln()
  {
    size++;
    return (size <= max_size + 2) || (pos_begin + size < match_string.size() + 1);
  }
 
public:  
  wc_star_token(const std::string & match_string_v, int max_size_v) :
                   wc_token(match_string_v), max_size(max_size_v) {};
  virtual void reset(int pos_begin_v)
  { 
    wc_token::reset(pos_begin_v);
    size = 0; 
  }
  virtual int next_token_pos()
  {
    //-1 because of size++ in fit next. 
    return pos_begin + size - 1;
  }
};

class wc_end_token : public wc_token
{
protected:
  virtual bool fit_next_intrln()
  {
    reseted = false;
    return pos_begin == match_string.size();
  }
    
public:
  wc_end_token(const std::string & match_string_v) :
                  wc_token(match_string_v){};

  virtual int next_token_pos()
  {
    return pos_begin + 1;
  }
  
};
/*
class wc_word_token : public wc_token
{
protected:
  std::string word;    
public:  
  wc_word_tocken(const std::string & match_string_v, const std::string &word_v) :
                   wc_token(match_string_v), word(word_v) {};
  virtual int next_token_pos()
  {
    return pos_begin + word.size();
  }
                   
  virtual bool fit()
  {
    _finish = true;
    bool result = pos_begin + word.size() < match_string.size();
    if(!result)
      return result;  
    
    return word == match_string.substr(pos_begin, word.size());
  }
};*/


class wc_char_token : public wc_token
{
protected:
  char char_match;  

  virtual bool fit_next_intrln()
  {
    reseted = false;
    bool result = pos_begin < match_string.size();
    if(!result)
      return result;  
    
    result = match_string[pos_begin] == char_match; 
    return result;
  }
  
public:
  wc_char_token(const std::string & match_string_v, char c) :
                  wc_token(match_string_v), char_match(c) {};
  
  virtual int next_token_pos()
  {
    return pos_begin + 1;
  }
};


class wc_any_char_token : public wc_token
{
protected:
  virtual bool fit_next_intrln()
  {
    reseted = false;
    return pos_begin < match_string.size();
  }
    
public:
  wc_any_char_token(const std::string & match_string_v) :
                  wc_token(match_string_v) {};

  virtual int next_token_pos()
  {
    return pos_begin + 1;
  }
};

bool match_wildcard(const std::string &wc_str, const std::string &cmp_str)
{
  std::vector<wc_token *> spv;
  
  for(int i = 0; i < wc_str.size(); i++)
  {    
    wc_token * wct;
    
    if(wc_str[i] == '*')
      wct = new wc_star_token(cmp_str, cmp_str.size());
    else if(wc_str[i] == '?')
      wct = new wc_any_char_token(cmp_str);
    else
      wct = new wc_char_token(cmp_str, wc_str[i]);
    
    spv.push_back(wct);
  }
  spv.push_back(new wc_end_token(cmp_str));
  
  int curr_token = 0;
  spv[0]->reset(0);
  while( (curr_token >= 0) && (curr_token < spv.size()) )
  {
    if( spv[curr_token]->fit_next() )
    {
      if( curr_token + 1 < spv.size() )
        spv[curr_token + 1]->reset(spv[curr_token]->next_token_pos());
      curr_token++;              
    }
    else
      curr_token--;
  }    
    
  for(int i = 0; i < spv.size(); i++)
    delete spv[i];
  
  return curr_token >= 0;
}

/*
bool match_wildcard(const std::string &wc_str, const std::string &cmp_str)
{
  std::vector<std::string> spv1 = split_string_vector(wc_str, "*");
  
  std::vector<std::string> spv;
  for(int i = 0; i < spv1.size(); i++)
  {
    spv.push_back(spv1[i]);
    spv.push_back("");
  }
  spv.pop_back();  
  
  if(wc_str[0] == '*')
    spv.insert(spv.begin(), "");

  if(wc_str[wc_str.length() - 1] == '*')
    spv.push_back("");
  
  std::vector<int> pos;
  pos.resize(spv.size(), 0);
  
  int curr_item = 0;
  do
  {
    int accepted_pos = -1;
    
    if(spv[curr_item] == "")
      accepted_pos = pos[curr_item];
    else
    {  
      for(int i = pos[curr_item]; i <= cmp_str.length(); i++)     
      {
        if(accept_pos(i, cmp_str, spv[curr_item]))
        {
          accepted_pos = i;
          break;
        }
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
*/