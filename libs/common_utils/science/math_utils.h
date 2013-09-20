/* 
 * File:   math_utils.h
 * Author: kirill
 *
 * Created on June 6, 2013, 2:22 PM
 */

#ifndef MATH_UTILS_H
#define	MATH_UTILS_H

#include <vector>

class math_utils 
{
public:  
  static std::vector<double> brd_Lorenzian(const std::vector<double> &data, const double gamma);
  static std::vector<double> brd_Gaussian(const std::vector<double> &data, const double sigma, const double precision = 4.0);
  //static std::vector<double> move_data(const std::vector<double> &data, const double gamma);  
};

#endif	/* MATH_UTILS_H */

