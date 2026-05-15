/* 
 * File:   comb_points.h
 * Author: kirill
 *
 * Created on December 16, 2013, 9:27 PM
 */

#ifndef COMB_POINTS_H
#define	COMB_POINTS_H

#include <vector>
#include <list>
#include <set>
#include <map>

class cmb_group
{
public:  
  std::set<int> indexes;
  double max_dist;
};

typedef std::vector<cmb_group> groups_vc;

typedef std::vector<std::vector<int> > index_conn;

typedef std::vector<double> vc_dist;
typedef std::vector<std::pair<double, int>> map_dist;

struct cmb_dist
{
  int index_cntr;
  vc_dist   dst_array;
  map_dist  dst_map;
};


class points_clusters
{
protected:
  double tol_list;
 protected:
  void get_possible_connections(index_conn &ic, double tol,
                                int min_cntr_points = 4) const;
  int  verify_connections(index_conn &ic) const;
  static void assign_groups(const index_conn &ic, groups_vc &grp);
  void create_groups_internal(groups_vc &vc, double tol_list_v, int min_cntr_points);  
public:
  virtual int get_points_size() const = 0;
  virtual double get_distance(int i, int j) const = 0;
  virtual bool is_connected(int i, int j) const
  {return get_distance(i, j) <= tol_list; };
public:
  virtual void assign_max_dist(groups_vc &gc);
  virtual double min_dist_between_groups(groups_vc &gc); 
};

#endif	/* COMB_POINTS_H */

