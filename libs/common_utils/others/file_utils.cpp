/* 
 * File:   file_utils.cpp
 * Author: kokho
 * 
 * Created on December 3, 2010, 12:18 PM
 */

#include "file_utils.h"

#include <sstream>
#include <cstdlib>
#include <stdio.h>
#include <iostream>

#include <boost/filesystem.hpp>

bool getLineFromFile(std::string &str, std::ifstream &f_in, const int max_line_size)
{
  return getLineFromStream(str, f_in, max_line_size);
}


bool getLineFromStream(std::string &str, std::istream &istr, const int max_line_size)
{
  char sLine[max_line_size];
  istr.getline(sLine, max_line_size);
  str = sLine;
  return !istr.eof();
}

long slength(std::istream &ist)
{
  long result;
  long current_pos = ist.tellg();
	
  ist.seekg(0, std::ios::end);
  result = ist.tellg();
  ist.seekg(current_pos, std::ios::beg);
  return result;	
}

std::string  readAllFiletoString(const std::string filename)
{
  std::ifstream t(filename.c_str());
  std::stringstream buffer;
  buffer << t.rdbuf();
  return buffer.str();
}

std::string check_folder(const std::string dir)
{
  std::string result = "";

  boost::system::error_code ec;
  
  namespace bfs = boost::filesystem;
          
  if( bfs::is_directory(dir.c_str(), ec) )
  {  
    if(ec.value() == 0)
      result = dir;
  }
  
  if(result != "")
    if(result[result.size() - 1] != '/')
      result += "/";
  
  return result;  
}

std::string get_temporary_file_name(std::string dir, 
                                    std::string prefix,
                                    std::string suffix)
{
  static int call_count = 0;
  const std::string alphabet = "AaBbCcDdEeFfGgHhIiJjKkLlMmNnOoPpQqRrSsTtUuVvWwXxYyZz";
  
  using namespace boost;
  std::string result;
  
  srand(time(NULL) + call_count);

  char * env_tmpdir = getenv("TMPDIR");
  if(dir == "" && env_tmpdir != NULL)
    dir = check_folder(env_tmpdir);

  if(dir == "")
    dir = check_folder(P_tmpdir);
  
  if(dir == "")
    dir = check_folder("/tmp");
  
  if(dir == "")
    dir = "./";
  
  for (int i = 0; i < 10000; i++) 
  {
    std::string rnd_name = "";
    
    for(int j = 0; j < 6; j++)
      rnd_name += alphabet[rand() % alphabet.size()];
    
    std::string file_name = dir + prefix + rnd_name + suffix;

    int acs = access(file_name.c_str(), F_OK);
    if( acs < 0 ) 
    {  
      if( errno == ENOENT ) 
      {
        result = file_name;
        break;
      }  
      else if (errno != EEXIST) 
      {
        result = "";
        break;
      }
    }  
  }
  
  call_count++;
  
  return result;
}        
