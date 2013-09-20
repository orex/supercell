/* 
 * File:   TStatisticCollect.h
 * Author: kokho
 *
 * Created on October 8, 2010, 3:20 PM
 */

#ifndef TSTATISTICCOLLECT_H
#define	TSTATISTICCOLLECT_H

#include <string>
#include <vector>

class TProcedureStatistic
{
public:
  virtual void proc(const std::vector<double> &xv, std::vector<double> &yv) = 0;
};

class TCommonStatisticCollect
{
private:
    bool out_range_value;
protected:
    std::vector<double> points;
    std::vector<double> points_xiyi;

    int   num_points();
    void  set_data_common(int num_points);
    void  add_point_common(int timestep, int value, double weight);
    bool  get_trend_at_point(int num_of_timepoints, int stat_point, double &trend_a, double &trend_b);
    
    virtual double delta_value() = 0;    
public:
    TCommonStatisticCollect();

    void   convert_to_mean_stat();
    void   convert_to_sum_stat();
    void   norm_by_value(double value);
    void   norm_by_distribution(TCommonStatisticCollect &ccd);
    void   norm_by_one();
    double norm();
    void   clear();
    void copy_data(std::vector<double> &pv);
    bool out_range();
};


class TIntStepStatisticCollect : public TCommonStatisticCollect
{
private:

protected:
    virtual double delta_value();
public:

    TIntStepStatisticCollect();
    TIntStepStatisticCollect(int num_points_v);
    void set_data(int num_points_v);
    int save_to_file(std::string file_name, int num_of_timepoints);
    void add_point(int timestep, int value, double weight);
};

class TRealStepStatisticCollect : public TCommonStatisticCollect
{
private:

protected:
    double min_value;
    double max_value;

public:

    TRealStepStatisticCollect();
    TRealStepStatisticCollect(double min_value_v, double max_value_v, int num_points_v);
    TRealStepStatisticCollect(double min_value_v, double max_value_v, double step);    
    void set_data(double min_value_v, double max_value_v, int num_points_v);
    int save_to_file(std::string file_name, int num_of_timepoints);
    void add_point(int timestep, double value, double weight);
    void procedure_stat(TProcedureStatistic &pc);
    virtual double delta_value();
};


#endif	/* TSTATISTICCOLLECT_H */

