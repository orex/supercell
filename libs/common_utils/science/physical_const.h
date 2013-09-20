/* 
 * File:   physical_const.h
 * Author: kirill
 *
 * Created on June 10, 2013, 1:55 PM
 */

#ifndef PHYSICAL_CONST_H
#define	PHYSICAL_CONST_H

#include <math.h>

namespace physical_const
{
  namespace SI
  {
    const double e_charge      = 1.6021765E-19;
    const double plank_nu      = 6.626069E-34;
    const double plank_w       = plank_nu / (2 * M_PI);
    const double e_mass        = 9.109382E-31;
    const double e0            = 8.854187817E-12;
  
    const double bohr_r        = 4 * M_PI * e0 * pow(plank_w, 2)/
                                 (e_mass * pow(e_charge, 2));
    const double hartree_eng   = e_mass * pow(e_charge, 4) / pow(4 * M_PI * e0 * plank_w, 2);
  }
};

#endif	/* PHYSICAL_CONST_H */

