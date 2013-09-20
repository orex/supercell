/* 
 * File:   nmr_utils.h
 * Author: kirill
 *
 * Created on June 6, 2013, 2:07 PM
 */

#ifndef NMR_UTILS_H
#define	NMR_UTILS_H

bool is_spin(double spin, bool &is_half);

class nmr_utils 
{
public:
  nmr_utils();
  nmr_utils(const nmr_utils& orig);
  virtual ~nmr_utils();
private:
};

#endif	/* NMR_UTILS_H */

