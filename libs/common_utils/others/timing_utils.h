/* 
 * File:   timing_utils.h
 * Author: kokho
 *
 * Created on March 13, 2011, 7:27 PM
 */

#ifndef TIMING_UTILS_H
#define	TIMING_UTILS_H

#include <time.h>

class timing_utils
{
protected:
  timeval tv;
public:
  timing_utils();
  timing_utils(const timing_utils& orig);
  double delta_t();
  virtual ~timing_utils();
};

#endif	/* TIMING_UTILS_H */

