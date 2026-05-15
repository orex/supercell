/* 
 * File:   string_utils.h
 * Author: kokho
 *
 * Created on May 2, 2012, 3:49 PM
 */

#ifndef STRING_UTILS_H
#define	STRING_UTILS_H

#include <string>
#include <vector>
#include <set>

std::string& trim(std::string &str);
std::string trim_delim(const std::string str, const std::string delimeter = " ");
std::vector<std::string> split_string_vector(std::string str, const std::string delimeter = "\t");
std::set<std::string> split_string_set(std::string str, const std::string delimeter = "\t");
bool is_double(const std::string& s, double& r_double);

std::string wildcard_to_regex(const std::string &wildcard);

bool match_wildcard(const std::string &wc_str, const std::string &cmp_str);

std::string get_index_str(const int index, const int index_max, 
                          const std::string filler = "0");

#endif	/* STRING_UTILS_H */

