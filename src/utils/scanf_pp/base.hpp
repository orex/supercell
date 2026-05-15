/* 
 * File:   base.hpp
 * Author: kirill
 *
 * Created on September 8, 2013, 3:53 PM
 */

#ifndef SCANF_PP_BASE_HPP
#define	SCANF_PP_BASE_HPP

#include <string>
#include <regex>
#include <boost/lexical_cast.hpp>
#include <cassert>


namespace scanf_pp
{
  class regex_scanf
  {
  private:
    std::string match_str;
    std::regex rx;
    std::smatch what;
    int curr_pos;
    
    bool all_vars_read()
    {
      return (curr_pos < 0) || (curr_pos == what.size());
    }
    
  public:
    explicit regex_scanf(const std::string &rx_str) : rx(rx_str), curr_pos(-1) { };
    explicit regex_scanf(const char * rx_str) : rx(rx_str), curr_pos(-1) { };
    
    bool regex_match(const std::string &rx_str)
    { 
      assert(all_vars_read());
      match_str = rx_str;
      bool result = std::regex_match(match_str, what, rx);
      
      curr_pos   = result ? 1 : -1;
      
      if(!result)
        match_str.clear();
      
      return result; 
    };
    
    bool regex_match(const char * rx_str)
    { return regex_match(std::string(rx_str)); };
    
    template<typename T>
    friend regex_scanf& operator>>(regex_scanf& left, T &right);
    friend regex_scanf& operator>>(regex_scanf& left, const void * right);
    
    virtual ~regex_scanf()
    { 
      assert(all_vars_read());
    }
  };
  
  template<typename T>
  regex_scanf& operator>>(regex_scanf& left, T &right)
  {
    assert(left.what.size() > left.curr_pos);
    right = boost::lexical_cast<T>(left.what[left.curr_pos]);
    left.curr_pos++;
    return left;
  }
  
  regex_scanf& operator>>(regex_scanf& left, const void * right)
  {
    assert(right == NULL);
    left.curr_pos++;
    return left;
  }
  
}

#endif	/* SCANF_PP_BASE_HPP */

