/* 
 * File:   TStatisticCollect.cpp
 * Author: kokho
 * 
 * Created on October 8, 2010, 3:20 PM
 */

#include "TStatisticCollect.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>

using namespace std;

//----------------------------------- TCommonStatisticCollect -----------------

TCommonStatisticCollect::TCommonStatisticCollect()
{
}

void TCommonStatisticCollect::set_data_common(int num_points_v)
{
  points.resize(num_points_v);
  points_xiyi.resize(num_points_v);
  out_range_value = false;
  clear();
}

void TCommonStatisticCollect::add_point_common(int timestep, int value, double weight)
{
  bool out_cur_range = (value < 0) || (value >= num_points());

  if(! out_cur_range )
  {
     points[value] += weight;
     points_xiyi[value] += weight * timestep;
  }

  out_range_value = out_range_value || out_cur_range;
}

bool TCommonStatisticCollect::get_trend_at_point(int num_of_timepoints, int stat_point, double &trend_a, double &trend_b)
{
  bool result = num_of_timepoints > 1;

  if ( result )
  {
    double notm = double(num_of_timepoints);
    double sum_xi  = notm * (notm - 1) / 2.0;
    double sum_xi2 = notm * (notm - 1) * (2 * notm - 1) / 6.0;
    double sum_yi = points[stat_point] * notm;
    double sum_xiyi = points_xiyi[stat_point] * notm;
    //a * sum_xi2 + b * sum_xi = sum_xiyi
    //a * sum_xi  + b * notm   = sum_yi
    double det   = sum_xi2  * notm   - sum_xi * sum_xi;
    double det_a = sum_xiyi * notm   - sum_yi * sum_xi;
    double det_b = sum_xi2  * sum_yi - sum_xi * sum_xiyi;

    assert(det);

    trend_a = det_a/det;
    trend_b = det_b/det;
  }

  return result;
}

bool TCommonStatisticCollect::out_range()
{
  return out_range_value;
}

int TCommonStatisticCollect::num_points()
{
  return points.size();
}

void TCommonStatisticCollect::copy_data(std::vector<double> &pv)
{
  pv.resize(num_points());
  for(int i = 0; i < num_points(); i++)
  {
    pv[i] = points[i];
  }
}

void TCommonStatisticCollect::clear()
{
  for(int i = 0; i < num_points(); i++)
  {
    points[i] = 0.0;
    points_xiyi[i] = 0.0;
  }
}

void TCommonStatisticCollect::norm_by_value(double value)
{
  for(long i = 0; i < num_points(); i++)
  {
    points[i] /= value;
    points_xiyi[i] /= value;
  }
}
void TCommonStatisticCollect::norm_by_distribution(TCommonStatisticCollect &ccd)
{
  assert(points.size() == ccd.points.size());
  for(long i = 0; i < num_points(); i++)
  {
    points[i] /= ccd.points[i];
    points_xiyi[i] /= ccd.points[i];
  }
}

void TCommonStatisticCollect::norm_by_one()
{
  norm_by_value(norm());
}

double TCommonStatisticCollect::norm()
{
  double result = 0;
  for(long i = 0; i < num_points(); i++)
    result += points[i];

  return result;
}

void TCommonStatisticCollect::convert_to_mean_stat()
{
  double sum_points = 0;
  double sum_pointsxy = 0;
  
  double sumr_points = 0;
  double sumr_pointsxy = 0;
  
  for(int i = 0; i < num_points(); i++)
  {
    sum_points    += points[i];
    sum_pointsxy  += points_xiyi[i];
  
    sumr_points   += delta_value() * i * points[i];
    sumr_pointsxy += delta_value() * i * points_xiyi[i];

    points[i] = sumr_points / sum_points;
    points_xiyi[i] += sumr_pointsxy / sum_pointsxy;
  }
}

void TCommonStatisticCollect::convert_to_sum_stat()
{
  for(int i = 1; i < num_points(); i++)
  {
    points[i] += points[i - 1];
    points_xiyi[i] += points_xiyi[i - 1];
  }
}


//----------------------------------- TIntStepStatisticCollect ----------------

TIntStepStatisticCollect::TIntStepStatisticCollect()
{

}

TIntStepStatisticCollect::TIntStepStatisticCollect(int num_points_v)
{
  set_data(num_points_v);
}

double TIntStepStatisticCollect::delta_value()
{
  return 1.0;
}

void TIntStepStatisticCollect::set_data(int num_points_v)
{
  set_data_common(num_points_v);
}

void TIntStepStatisticCollect::add_point(int timestep, int value, double weight)
{
  add_point_common(timestep, value, weight);
}

int TIntStepStatisticCollect::save_to_file(std::string file_name, int num_of_timepoints)
{
  ofstream ofile(file_name.c_str(), fstream::out);

  if(! ofile)
    return -2;

  for(long i = 0; i < num_points(); i++)
  {
    double t_a, t_b;
    get_trend_at_point(num_of_timepoints, i, t_a, t_b);

    ofile << i << "\t"
          << points[i] << "\t"
          << t_b << "\t"
          << (t_b + t_a * double(num_of_timepoints - 1)) << endl;
  }

  return 0;
}


//----------------------------------- TRealStepStatisticCollect ----------------

TRealStepStatisticCollect::TRealStepStatisticCollect()
{

}

TRealStepStatisticCollect::TRealStepStatisticCollect(double min_value_v, double max_value_v, int num_points_v)
{
  set_data(min_value_v, max_value_v, num_points_v);
}

TRealStepStatisticCollect::TRealStepStatisticCollect(double min_value_v, double max_value_v, double step)
{
  set_data(min_value_v, max_value_v, (max_value_v - min_value_v) / step );
}

void TRealStepStatisticCollect::set_data(double min_value_v, double max_value_v, int num_points_v)
{
  set_data_common(num_points_v);
  min_value = min_value_v;
  max_value = max_value_v;
}

int TRealStepStatisticCollect::save_to_file(std::string file_name, int num_of_timepoints)
{
  ofstream ofile(file_name.c_str(), fstream::out);
  if(! ofile)
    return -2;

  for(long i = 0; i < num_points(); i++)
  {
    double t_a, t_b;
    get_trend_at_point(num_of_timepoints, i, t_a, t_b);
    ofile << min_value + double(i) * delta_value() << "\t"
          << points[i] << "\t"
          << t_b << "\t"
          << (t_b + t_a * double(num_of_timepoints - 1)) << endl;
  }

  return 0;
}

void TRealStepStatisticCollect::procedure_stat(TProcedureStatistic &pc)
{
  std::vector<double> xv;
  
  for(int i = 0; i < num_points(); i++)
    xv.push_back(min_value + double(i) * delta_value());
  
  pc.proc(xv, points);
}

void TRealStepStatisticCollect::add_point(int timestep, double value, double weight )
{
  long cur_point = floor(double(num_points()) * (value - min_value) / (max_value - min_value));

  add_point_common(timestep, cur_point, weight);
}

double TRealStepStatisticCollect::delta_value()
{
  return (max_value - min_value)/ double(num_points());
}



