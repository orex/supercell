#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE string_utils_test

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cstdlib>
#include <ctime>

#include "string_utils.h"

BOOST_AUTO_TEST_SUITE(StringUtilsTest)

BOOST_AUTO_TEST_CASE(Test_Wildcard)
{
  BOOST_CHECK(match_wildcard("", ""));
  BOOST_CHECK(!match_wildcard("", "1"));
    
  BOOST_CHECK(match_wildcard("*", ""));
  BOOST_CHECK(match_wildcard("*", "1"));
  
  BOOST_CHECK(match_wildcard("**", "1"));
  BOOST_CHECK(match_wildcard("*?", "1"));
  BOOST_CHECK(!match_wildcard("*?", ""));
  
  BOOST_CHECK(!match_wildcard("Ca*x", "Ca1x1"));
  BOOST_CHECK(match_wildcard("Ca*", "Ca"));
  BOOST_CHECK(!match_wildcard("Ca1*", "Ca"));  
  BOOST_CHECK(match_wildcard("Ca*", "Ca1"));    
  BOOST_CHECK(match_wildcard("Ca*x", "Ca1x")); 
  BOOST_CHECK(match_wildcard("*Ca?x*", "Ca1xq")); 
  
}

BOOST_AUTO_TEST_SUITE_END()
