/* 
 * File:   comb_points.cpp
 * Author: kirill
 * 
 * Created on December 16, 2013, 9:27 PM
 */

#include "comb_points.h"

#include <map>
#include <cassert>
#include <cmath>
#include <algorithm>

using namespace std;

void points_clusters::get_possible_connections(index_conn &ic, double tol,
                                               int min_cntr_points) const
{
  ic.resize(get_points_size());
  if( get_points_size() < 10 || min_cntr_points == 0 ) {
    for(int i = 0; i < ic.size(); i++) {
      ic[i].clear();
      for(int j = 0; j < ic.size(); j++) {
        if( i != j )
          ic[i].emplace_back(j);
      }
    }
    return;
  }

  vector<cmb_dist> dist_cntr;

  min_cntr_points = std::min<int>(min_cntr_points, get_points_size());

  int ci = 0;
  for(int i = 0; i < min_cntr_points; i++) {
    dist_cntr.emplace_back();
    auto & cd = dist_cntr.back();
    cd.index_cntr = ci;
    cd.dst_array.resize(get_points_size());
    cd.dst_map.resize(get_points_size());
    pair<double, int> max_dist(0, 0);
    for(int j = 0; j < get_points_size(); j++) {
      double dst = get_distance(ci, j);
      cd.dst_array[j] = dst;
      cd.dst_map[j] = {dst, j};
      double td = 0;
      for(int k = 0; k < dist_cntr.size(); k++) {
        td += std::sqrt(get_distance(dist_cntr[k].index_cntr, j));
      }
      max_dist = std::max(max_dist, {td, j});
    }
    sort(cd.dst_map.begin(), cd.dst_map.end());
    ci = max_dist.second;
  }
  
  std::vector<int> ms;
  ms.resize(ic.size(), 0);
  std::vector<std::pair<map_dist::const_iterator, map_dist::const_iterator> > vm(dist_cntr.size());
  for(int i = 0; i < ic.size(); i++) {
    for(int j = 0; j < dist_cntr.size(); j++) {
      double dst = dist_cntr[j].dst_array[i];
      const auto & cd = dist_cntr[j].dst_map;
      auto itl = lower_bound(cd.cbegin(), cd.cend(), make_pair(dst - tol, 0));
      auto ith = upper_bound(itl, cd.cend(), make_pair(dst + tol, 0));
      vm[j] = {itl, ith};
    }
    sort(vm.begin(), vm.end(), [](const auto a, const auto b) {
      return std::distance(a.first, a.second) < std::distance(b.first, b.second);
    });

    int mxc = std::min<int>(3, vm.size());
    for(int j = 0; j < mxc; j++) {
      for(auto it = vm[j].first; it != vm[j].second; ++it)
        if( j == 0 )
          ms[it->second] = 1;
        else
          ms[it->second]++;
    }
    ic[i].clear();
    ic[i].reserve(std::distance(vm.front().first, vm.front().second));
    for(auto it = vm.front().first; it != vm.front().second; ++it) {
      if( i != it->second && ms[it->second] == mxc )
        ic[i].emplace_back(it->second);
    }
  }
}

int points_clusters::verify_connections(index_conn &ic) const {
  int result = 0;

  for (int i = 0; i < ic.size(); i++) {
    auto it = std::remove_if(ic[i].begin(), ic[i].end(),
                             [this, i](int a) -> bool {
                               return !is_connected(i, a);
                             });
    ic[i].erase(it, ic[i].end());
    result += ic[i].size();
  }
  return result;
}

void points_clusters::assign_groups(const index_conn &ic, groups_vc &grp)
{
  vector<int> group_num;
  map<int, int> group_st; 
  
  group_num.resize(ic.size());
  for(int i = 0; i < group_num.size(); i++)
    group_num[i] = i;
  
  //Set groups num by lowers connection  
  bool changed;
  do
  {  
    changed = false;
    for(int i = 0; i < ic.size(); i++)
    {
      // // Avoid connections down and self connections
      //set<int>::const_iterator b_it = upper_bound(ic[i].begin(), ic[i].end(), i);
      int min_grp_num = group_num[i];
      int max_grp_num = group_num[i];
      for(auto it = ic[i].cbegin(); it != ic[i].cend(); ++it)
      {  
        min_grp_num = min(min_grp_num, group_num[*it]);
        max_grp_num = max(max_grp_num, group_num[*it]);
      }
      assert(min_grp_num <= max_grp_num);
      if( min_grp_num <  max_grp_num ) 
      {  
        for(auto it = ic[i].cbegin(); it != ic[i].cend(); ++it)
          group_num[*it] = min_grp_num;
        changed = true;
      }  
    }
  }
  while( changed );
  
  //set groups num 0..1..grp_num - 1
  int grp_num = 0;
  for(int i = 0; i < group_num.size(); i++)
  {
    if(group_st.count(group_num[i]) == 0)
    {
      group_st[group_num[i]] = grp_num;
      grp_num++;
    }  
  }
  grp.resize(grp_num);
  
  for(int i = 0; i < group_num.size(); i++)
  {
    int index_group = group_st[group_num[i]];
    grp[index_group].indexes.insert(i);
  }

  for(int i = 0; i < grp.size(); i++)
  {
    int min_con = ic[*grp[i].indexes.begin()].size();
    int max_con = min_con;
    
    for(set<int>::const_iterator it  = grp[i].indexes.begin(); 
                                 it != grp[i].indexes.end(); ++it)
    {  
      min_con = min(min_con, int(ic[*it].size()) );
      max_con = max(max_con, int(ic[*it].size()) );
    }
    
    assert(min_con <= max_con);
    assert(max_con < grp[i].indexes.size());
    
  }
}

void points_clusters::create_groups_internal(groups_vc &vc, double tolerance_v, int min_cntr_points)
{
  tol_list = tolerance_v;
  
  index_conn ic;
  get_possible_connections(ic, tol_list, min_cntr_points);
  verify_connections(ic);
  assign_groups(ic, vc);
}

void points_clusters::assign_max_dist(groups_vc &gc)
{
  for(int i = 0; i < gc.size(); i++)
  {
    double md = 0;
    for(set<int>::const_iterator it  = gc[i].indexes.begin(); 
                                 it != gc[i].indexes.end(); ++it)
    {
      for(set<int>::const_iterator jt  = gc[i].indexes.begin(); 
                                   jt != gc[i].indexes.end(); ++jt)
        md = max(md, get_distance(*it, *jt));
    }  
    gc[i].max_dist = md;
  }
}

double points_clusters::min_dist_between_groups(groups_vc &gc)
{
  //TODO: Make it faster!!!
  double result = -1;
  for(int i = 0; i < gc.size(); i++)
  {
    for(set<int>::const_iterator it  = gc[i].indexes.begin(); 
                                 it != gc[i].indexes.end(); ++it)
    {
      for(int j = i + 1; j < gc.size(); j++)
      {
        for(set<int>::const_iterator jt  = gc[j].indexes.begin(); 
                                     jt != gc[j].indexes.end(); ++jt)
        {  
          if(result >= 0)
            result = min(result, get_distance(*it, *jt));
          else
            result = get_distance(*it, *jt);
        }  
      }
    }  
  }

  return result;  
}

