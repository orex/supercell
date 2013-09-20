/* 
 * File:   nmr_utils.cpp
 * Author: kirill
 * 
 * Created on June 6, 2013, 2:07 PM
 */

#include "nmr_utils.h"

#include <cmath>

bool is_spin(double spin, bool &is_half)
{
  bool result;
  
  result = (spin > 0) && (spin <= 4.5) && (spin - floor(2 * spin)< 1.0E-6 ) ;
  
  if (result)
    is_half = spin-floor(spin) > 0.4 ;
  
  return result ;
}

nmr_utils::nmr_utils()
{
}

nmr_utils::nmr_utils(const nmr_utils& orig)
{
}

nmr_utils::~nmr_utils()
{
}

