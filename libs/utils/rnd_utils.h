/* 
 * File:   rnd_utils.h
 * Author: kokho
 *
 * Created on February 7, 2011, 5:19 PM
 */

#ifndef RND_UTILS_H
#define	RND_UTILS_H

#include <boost/random.hpp>

boost::mt19937 create_rnd_gen();
double get_rnd_value_in_interval(boost::mt19937 &rng, double min, double max);
int get_rnd_int_value_in_interval(boost::mt19937 &rng, int min, int max);

#endif	/* RND_UTILS_H */

