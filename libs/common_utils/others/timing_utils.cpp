/* 
 * File:   timing_utils.cpp
 * Author: kokho
 * 
 * Created on March 13, 2011, 7:27 PM
 */

#include <stdlib.h>

#include "timing_utils.h"
#include <sys/time.h>

timing_utils::timing_utils()
{
  gettimeofday(&tv, NULL);
}

timing_utils::timing_utils(const timing_utils& orig)
{
}

double timing_utils::delta_t()
{
  timeval tv1, dtv;
  gettimeofday(&tv1, NULL);
  timersub(&tv1, &tv, &dtv);

  return dtv.tv_sec + double(dtv.tv_usec) * 1E-6;
}

timing_utils::~timing_utils() {
}

