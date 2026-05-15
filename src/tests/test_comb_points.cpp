#define BOOST_TEST_MODULE comb_points_test

#include <boost/test/unit_test.hpp>

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>

#include "cryst_tools/comb_points.h"

using namespace std;

class points_ND : public points_clusters
{
protected:
  vector< vector<double> > vc;
public:
  int get_points_size() const override;
  double get_distance(int i, int j) const override;
  points_ND(int num_of_points, int dimension);
  void create_groups(groups_vc &out, double tol_list_v, int min_cntr_points)
  { create_groups_internal(out, tol_list_v, min_cntr_points); }
};

points_ND::points_ND(int num_of_points, int dimension)
{
  vc.resize(num_of_points);
  for(size_t i = 0; i < vc.size(); i++)
  {
    vc[i].resize(dimension);
    for(int j = 0; j < dimension; j++)
      vc[i][j] = double(rand()) / double(RAND_MAX);
  }
}

int points_ND::get_points_size() const
{
  return vc.size();
}

double points_ND::get_distance(int i, int j) const
{
  assert(vc[i].size() == vc[j].size());
  double result = 0;
  for(size_t k = 0; k < vc[i].size(); k++)
    result += (vc[i][k] - vc[j][k]) * (vc[i][k] - vc[j][k]);
  return std::sqrt(result);
}

static bool groups_equivalent(const groups_vc &g1, const groups_vc &g2)
{
  if (g1.size() != g2.size())
    return false;
  for(size_t i = 0; i < g1.size(); i++)
  {
    if (g1[i].indexes != g2[i].indexes)
      return false;
    if (g1[i].max_dist != g2[i].max_dist)
      return false;
  }
  return true;
}

static bool groups_intersect(const points_ND &pt, const cmb_group &g1,
                             const cmb_group &g2, const double tol)
{
  for(auto it = g1.indexes.begin(); it != g1.indexes.end(); ++it)
    for(auto jt = g2.indexes.begin(); jt != g2.indexes.end(); ++jt)
      if (pt.get_distance(*it, *jt) < tol) {
        cout << "intersect on " << *it << " - " << *jt;
        return true;
      }
  return false;
}

static bool group_connected(const points_ND &pt, cmb_group g, double tol)
{
  vector<int> cn_ind;
  cn_ind.push_back(*g.indexes.begin());
  g.indexes.erase(g.indexes.begin());

  bool changed;
  do {
    changed = false;
    for(int i = int(cn_ind.size()) - 1; i >= 0; i--) {
      for(auto it = g.indexes.begin(); it != g.indexes.end(); ++it) {
        if (pt.get_distance(*it, cn_ind[i]) < tol) {
          cn_ind.push_back(*it);
          g.indexes.erase(it);
          changed = true;
          break;
        }
      }
      if (changed) break;
    }
  } while (changed);

  return g.indexes.empty();
}

BOOST_AUTO_TEST_SUITE(CombPointsTest)

BOOST_AUTO_TEST_CASE(Test_random_points_comb_3D)
{
  const double tol = 0.1;
  time_t time_rand = time(NULL);
  srand(time_rand);
  cout << "Random initialization. Remember if test fails. " << endl;
  cout << "  " << time_rand << endl;

  const int points[] = {0, 1, 10, 100, 500, 1000, 1500, 2000, 3000, 5000};
  const int points_size = sizeof(points) / sizeof(points[0]);

  for(int i = 0; i < points_size; i++)
  {
    points_ND pt(points[i], 3);
    groups_vc gvc_ref;
    for(int j = 1; j < 7; j++) {
      groups_vc gvc;
      pt.create_groups(gvc, tol, j);
      pt.assign_max_dist(gvc);
      cout << "Num reference points: " << j << endl;
      cout << "Number of groups: " << gvc.size() << endl;
      if (j > 1)
        BOOST_CHECK(groups_equivalent(gvc, gvc_ref));
      gvc_ref = gvc;
    }

    BOOST_CHECK_LE(gvc_ref.size(), size_t(points[i]));

    multiset<int> ms;
    for(size_t j = 0; j < gvc_ref.size(); j++)
      ms.insert(gvc_ref[j].indexes.begin(), gvc_ref[j].indexes.end());
    BOOST_CHECK_EQUAL(ms.size(), size_t(points[i]));

    for(int j = 0; j < points[i]; j++)
      BOOST_CHECK(ms.count(j) == 1);

    for(size_t j = 0; j < gvc_ref.size(); j++)
      for(size_t k = j + 1; k < gvc_ref.size(); k++)
        BOOST_CHECK_EQUAL(groups_intersect(pt, gvc_ref[j], gvc_ref[k], tol), false);

    for(size_t j = 0; j < gvc_ref.size(); j++)
      BOOST_CHECK(group_connected(pt, gvc_ref[j], tol));
  }
}

BOOST_AUTO_TEST_SUITE_END()
