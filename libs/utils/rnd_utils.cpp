/* 
 * File:   rnd_utils.cpp
 * Author: kokho
 * 
 * Created on February 7, 2011, 5:19 PM
 */

#include "rnd_utils.h"
#include <sys/time.h>

boost::mt19937 create_rnd_gen()
{
  boost::mt19937 result;
  struct timeval tv;
  gettimeofday(&tv, NULL);
  result.seed(static_cast<unsigned> (tv.tv_usec));

  return result;
}

double get_rnd_value_in_interval(boost::mt19937 &rng, double min, double max)
{
  double result;

  boost::uniform_real<double> uniform(min, max);
  boost::variate_generator< boost::mt19937&, boost::uniform_real<double> > uniform_sampler(rng, uniform);

  result = uniform_sampler();

  return result;
}

int get_rnd_int_value_in_interval(boost::mt19937 &rng, int min, int max)
{
  int result;

  boost::uniform_int<int> uniform(min, max);
  boost::variate_generator< boost::mt19937&, boost::uniform_int<int> > uniform_sampler(rng, uniform);

  result = uniform_sampler();

  return result;
}



