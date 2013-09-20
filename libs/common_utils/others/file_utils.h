/* 
 * File:   file_utils.h
 * Author: kokho
 *
 * Created on December 3, 2010, 12:18 PM
 */

#ifndef FILE_UTILS_H
#define	FILE_UTILS_H

#include <string>
#include <fstream>

long slength(std::istream &ist);

bool getLineFromStream(std::string &str, std::istream &istr, const int max_line_size = 500);

bool getLineFromFile(std::string &str, std::ifstream &f_in, const int max_line_size = 500);

std::string  readAllFiletoString(const std::string filename);


/**
 * The function creates unique filename according to template.
 * @param dir temporary directory, if empty - 
 * 
 * @return name of the file with path
 */
std::string get_temporary_file_name(std::string dir = "", 
                                    std::string prefix = "",
                                    std::string suffix = "");
#endif	/* FILE_UTILS_H */

