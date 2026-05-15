/* 
 * File:   cryst_tools.cpp
 * Author: kirill
 * 
 * Created on December 16, 2013, 4:00 PM
 */

#include <cassert>
#include <set>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "cryst_tools.h"
#include "comb_points.h"
#include "science/math_constants.h"

using namespace std;
using namespace Eigen;
using namespace cryst_tools;

int Matrix_compare(const Matrix3d &m1, const Matrix3d &m2)
{
  int result = 0;
  for(int i = 0; i < 9; i++)  
  {
    double d = m1(i / 3, i % 3) - m2(i / 3, i % 3);
    if( abs(d) > 1E-4)
    {
      result = int(1.1 * d / abs(d));
      break;
    }  
  }
  
  return result;
}

int Vector_compare(const Vector3d &v1, const Vector3d &v2)
{
  int result = 0;
  for(int i = 0; i < 3; i++)  
  {
    double d = v1(i) - v2(i);
    if( abs(d) > 1E-4)
    {
      result = int(1.1 * d / abs(d));
      break;
    }  
  }
  
  return result;
}

struct Matrix3d_comp_cl 
{
  bool operator() (const Matrix3d &m1, const Matrix3d &m2) const
  {  return Matrix_compare(m1, m2) < 0;}
};

struct Vector3d_comp_cl 
{
  bool operator() (const Vector3d &v1, const Vector3d &v2) const
  {  return Vector_compare(v1, v2) < 0;}
};


Eigen::Vector3d cryst_tools::norm_frac(const Eigen::Vector3d &frac)
{
  Vector3d result;

  for(int i = 0; i < 3; i++)
    result[i] = frac[i] - floor(frac[i]);
  
  return result;
}


std::vector<Eigen::Matrix3d> cryst_tools::get_cell_symmetries(const Eigen::Matrix3d &cell, const double tol)
{
  std::vector<Eigen::Matrix3d> result;
  set<Eigen::Matrix3d, Matrix3d_comp_cl> real_space_symm;

  Vector3d ve = SelfAdjointEigenSolver<Matrix3d>(cell.transpose() * cell).eigenvalues();
  pair<double, double> mm = minmax({ve.x(), ve.y(), ve.z()});
  const int range = sqrt(std::abs(mm.second / mm.first)) + 1;
  const int mult = 2 * range + 1;

  std::array<double, 3> cell_lengths{cell.col(0).norm(), cell.col(1).norm(), cell.col(2).norm()};
  std::array<vector<Vector3i>, 3> conv_vect;
  for(int k = 0; k < mult * mult * mult; k++)
  {
    Vector3i vf;
    int p = k;
    for(int i = 0; i < 3; i++)
    {
      vf[i] = p % mult - range;
      p /= mult;
    }
    double nm = (cell * vf.cast<double>()).norm();
    for(int i = 0; i < 3; i++)
    {
      if(std::abs(nm - cell_lengths[i]) < tol ) {
        conv_vect[i].push_back(vf);
      }
    }
  }

  Matrix3d r_cell = cell.inverse();

  Matrix3i B;
  for(int i = 0; i < conv_vect[0].size(); i++)
  {
    B.col(0) = conv_vect[0][i];
    for(int j = 0; j < conv_vect[1].size(); j++)
    {
      B.col(1) = conv_vect[1][j];
      for(int k = 0; k < conv_vect[2].size(); k++)
      {
        B.col(2) = conv_vect[2][k];
        if( std::abs(B.determinant()) != 1) continue;
        Matrix3d A = cell * B.cast<double>() * r_cell;

        if((A.transpose() * A - Matrix3d::Identity()).norm() < tol)
        {
          if(real_space_symm.count(A) == 0)
          {
            real_space_symm.insert(A);
            result.push_back(B.cast<double>());
          }
        }
      }
    }
  }

  return result;
}

class ps_shifts : public points_clusters 
{
protected:
  vector< Vector3d > vc;
  vector< int > base_index;
  Matrix3d cell;
  double tol_cart;
public:
  virtual int get_points_size() const;
  virtual double get_distance(int i, int j) const;
  virtual bool is_connected(int i, int j) const;
public:
  explicit ps_shifts(const Matrix3d &cell_v);
  void ps_add_shift(int index, const Vector3d &shift);
  Vector3d average_vector(const cmb_group &cbg);
  bool total_shift(const cmb_group &cbg, int poins_num) const;
  void create_groups(groups_vc &vc, double tol_cart_v);    
  virtual void assign_max_dist(groups_vc &gc);  
};

ps_shifts::ps_shifts(const Matrix3d &cell_v)
{
  base_index.resize(0);
  vc.resize(0);
  cell = cell_v;
}

int ps_shifts::get_points_size() const
{
  return vc.size();
}

double ps_shifts::get_distance(int i, int j) const
{
  return min_frac(vc[i] - vc[j]).norm();
}

bool ps_shifts::is_connected(int i, int j) const
{
  return (cell * min_frac(vc[i] - vc[j])).norm() <= tol_cart;
}

void ps_shifts::create_groups(groups_vc &vc, double tol_cart_v)
{
  tol_cart = tol_cart_v;
  
  Matrix3d r_cell = cell.inverse();
  
  Matrix3d mc = r_cell.transpose() * r_cell;

  double sum = 0;
  for(int i = 0; i < 9; i++)  
  {
    sum += abs(mc(i / 3, i % 3));
  }
  tol_list = sqrt(sum) * tol_cart;  
  create_groups_internal(vc, tol_list, 10);
}

void ps_shifts::ps_add_shift(int index, const Vector3d &shift)
{
  vc.push_back(min_frac(shift));
  base_index.push_back(index);
}

Vector3d ps_shifts::average_vector(const cmb_group &cbg)
{
  Vector3d result;
  
  result.setZero();
  
  for(set<int>::const_iterator it  = cbg.indexes.begin();
                               it != cbg.indexes.end(); ++it)
    result += min_frac(vc[*it] - vc[*cbg.indexes.begin()]);

  result = vc[*cbg.indexes.begin()] + result / double(cbg.indexes.size());
  
  return result;
}

bool ps_shifts::total_shift(const cmb_group &cbg, int poins_num) const
{
  bool result;
  
  result = cbg.indexes.size() == poins_num;
  
  #ifndef NDEBUG
  if(result)
  {
    vector<int> vi;
    for(set<int>::const_iterator it  = cbg.indexes.begin();
                                 it != cbg.indexes.end(); ++it)
      vi.push_back(base_index[*it]);
    
    sort(vi.begin(), vi.end());
    assert(vi.size() == poins_num);
    for(int i = 0; i < vi.size(); i++)
      assert(vi[i] == i);
  }  
  
  #endif  
  
  return result;
}

void ps_shifts::assign_max_dist(groups_vc &gc)
{
  for(int i = 0; i < gc.size(); i++)
  {
    double md = 0;
    for(set<int>::const_iterator it  = gc[i].indexes.begin(); 
                                 it != gc[i].indexes.end(); ++it)
    {
      for(set<int>::const_iterator jt  = gc[i].indexes.begin(); 
                                   jt != gc[i].indexes.end(); ++jt)
        md = max(md, (cell * min_frac(vc[*it] - vc[*jt])).norm());
    }  
    gc[i].max_dist = md;
  }  
}

std::vector<Eigen::Vector3d> cryst_tools::get_shifts(const Eigen::Matrix3d &sym_matrix_B,
                                                     const Eigen::Matrix3d &cell,
                                                     const std::vector<Eigen::Vector3d> &frac_coords,
                                                     double tolerance, bool all_symmetric)
{
  ps_shifts ps(cell);
  
  for(int i = 0; i < frac_coords.size(); i++)  
  {
    for(int j = 0; j < frac_coords.size(); j++)
    {
      Vector3d shift;
      shift = frac_coords[j] - sym_matrix_B * frac_coords[i];
      ps.ps_add_shift(i, shift);
    }  
  }  
  
  groups_vc gvc;
  
  ps.create_groups(gvc, tolerance);
  ps.assign_max_dist(gvc);

  vector<Vector3d> result;
  for(int i = 0; i < gvc.size(); i++)
  {  
    assert(gvc[i].max_dist < tolerance);
    if(!all_symmetric)
      result.push_back(ps.average_vector(gvc[i]));
    else if( ps.total_shift(gvc[i], frac_coords.size()) )
      result.push_back(ps.average_vector(gvc[i]));
  }
  
  result.shrink_to_fit();
    
  return result;
}

std::vector<Eigen::Vector3d> cryst_tools::shifts_intersect(const Eigen::Matrix3d &cell,  
                                                           const std::vector<Eigen::Vector3d> &frac_c1,
                                                           const std::vector<Eigen::Vector3d> &frac_c2,
                                                           double tolerance)
{
  std::vector<Eigen::Vector3d> result;
  result.clear();
  
  for(int i = 0; i < frac_c1.size(); i++)
  {
    for(int j = 0; j < frac_c2.size(); j++)
    {
      Vector3d dist = min_frac(frac_c2[j] - frac_c1[i]);
      if((cell * dist).norm() <= tolerance)
        result.push_back(frac_c1[i] + 0.5 * dist);
    }
  }  
  
  return result;
}

vector_Affine3d cryst_tools::get_all_symmetries(const Eigen::Matrix3d &cell,
                                                             const vc_sets &frac_coords,
                                                             const std::vector<bool> &all_symm, 
                                                             const double tol)
{
  assert(all_symm.size() == frac_coords.size());

  vector_Affine3d result;
  result.clear();
  
  vector<Matrix3d> symm = get_cell_symmetries(cell, tol);

  for(int i = 0; i < symm.size(); i++)
  {
    vector<Vector3d> shifts_total;
    for(int j = 0; j < frac_coords.size(); j++)
    {  
      vector<Vector3d> shifts = get_shifts(symm[i], cell, frac_coords[j], tol, all_symm[j]);

      if(j == 0)
        shifts_total = shifts;
      else
        shifts_total = shifts_intersect(cell, shifts_total, shifts, tol);
    }
    
    for(int j = 0; j < shifts_total.size(); j++)
    {
      Affine3d mo;
      
      mo.setIdentity();
      
      mo.linear() = symm[i];
      mo.translation() = shifts_total[j];
      
      result.push_back(mo);
    }
  }
  return result;  
}



vector_Affine3d cryst_tools::get_all_symmetries(const Eigen::Matrix3d &cell,
                                                             const vc_sets &frac_coords,
                                                             const double tol)
{
  vector<bool> bc;
  bc.resize(frac_coords.size(), true);
  
  return get_all_symmetries(cell, frac_coords, bc, tol);
}        

void cryst_tools::min_dist::set_cell(const Eigen::Matrix3d &cell_v, double tol)
{
  box_size.setZero();
  
  cell = cell_v;
  r_cell = cell.inverse();
  
  Matrix3d Q = cell.inverse().transpose() * cell.inverse();
  
  r2_dist_direct = 1E12;
  
  for(int i = -4; i < 5; i++)
  {
    for(int j = -4; j < 5; j++)
    {
      for(int k = -4; k < 5; k++)
      {
        Vector3i vi_ind(i, j, k);
        Vector3d v_ind(i, j, k);
        Vector3d Qv = Q * v_ind;
        double d = v_ind.transpose() * Qv - (abs(Qv[0]) + abs(Qv[1]) + abs(Qv[2]));
        
        if(vi_ind.squaredNorm() != 0)
          r2_dist_direct = min(r2_dist_direct, (cell * v_ind).squaredNorm());
        
        if( d < -tol )
        {  
          for(int l = 0; l < 3; l++)
            box_size[l] = max(box_size[l], abs(vi_ind[l]));
        }  
      }  
    }  
  }
  
  cell_assigned = true;  
}

Eigen::Vector3d cryst_tools::min_dist::operator ()(const Eigen::Vector3d &vd) const
{
  assert(cell_assigned);
  
  Eigen::Vector3d result = cell * min_frac(r_cell * vd);
  
  double sq_norm_min = vd.squaredNorm();
  
  if(sq_norm_min <= r2_dist_direct )
    return result;
  
  for(int i = -box_size[0]; i < box_size[0] + 1; i++)
  {
    for(int j = -box_size[1]; j < box_size[1] + 1; j++)
    {
      for(int k = -box_size[2]; k < box_size[2] + 1; k++)
      {
        if(abs(i) + abs(j) + abs(k) == 0 )
          continue;
        
        Vector3d vc = vd - cell * Vector3d(i, j, k);
        double d = vc.squaredNorm();
        if(d < sq_norm_min )
        {
          sq_norm_min = d;
          result = vc;
        }  
      }
    }
  }  
  
  return result;
}

std::vector<Eigen::Vector3d> cryst_tools::min_dist::get_img_dist(const Eigen::Vector3d &vd, 
                                                                 const Eigen::Vector3i &box) const
{
  assert(cell_assigned);

  std::vector<Eigen::Vector3d> result;
  
  for(int i = -box[0]; i <= box[0]; i++)
  {
    for(int j = -box[1]; j <= box[1]; j++)
    {
      for(int k = -box[2]; k <= box[2]; k++)
      {
        Vector3d vc = vd - cell * Vector3d(i, j, k);
        result.push_back(vc);
      }
    }
  }  
  
  return result;
}        

Eigen::Vector3d cryst_tools::min_dist::average_vector(const std::vector<Eigen::Vector3d> &av) const
{
  assert(cell_assigned);
  
  Vector3d result;
  
  result.setZero();
  
  for(int i = 0; i < av.size(); i++)
    result += (*this)(av[i] - av[0]);

  result = av[0] + result / double(av.size());
  result = (*this)(result);
  
  return result;
  
}

//Ewald sum
void cryst_tools::ewald_sum::set_cell(const Eigen::Matrix3d &cell_v)
{
  cell = cell_v;
  res_cell = 2 * cnst::pi * (cell.inverse()).transpose();
  volume = std::abs(cell.determinant());
}

void cryst_tools::ewald_sum::set_precision(double N, double precision_v)
{
  constexpr int max_periods = (2 * max_cell_check_range + 1) *
                              (2 * max_cell_check_range + 1) *
                              (2 * max_cell_check_range + 1);

  precision = precision_v;
  
  // pow(1/V, 2/3) == 1 / V^(2/3) == 1 / cbrt(V*V). Using cbrt avoids the libm
  // pow symbol (cbrt is at GLIBC_2.2.5, pow's default moved to 2.29 in 2019).
  eta = cnst::pi / std::cbrt(volume * volume);
  sqeta = sqrt(eta);

  for(r_max = 1 / sqeta;
       std::erfc(r_max * sqeta) / r_max > precision / max_periods;
       r_max += 0.1 / sqeta) {};

  double max_input_shift = 2 * (cell.col(0).norm() +
                     cell.col(1).norm() + cell.col(2).norm());

  real_data_shifts.clear();
  rec_data.clear();
  for(int i = -max_cell_check_range; i <= max_cell_check_range; i++) {
    for (int j = -max_cell_check_range; j <= max_cell_check_range; j++) {
      for (int k = -max_cell_check_range; k <= max_cell_check_range; k++) {
        Vector3d p(i, j, k);
        Vector3d shift = cell * p;
        if( shift.norm() - max_input_shift < r_max )
          real_data_shifts.emplace_back(shift);

        if( (i != 0) || (j != 0) || (k != 0) ) {
          Vector3d recl_K = res_cell * p;
          double K2 = recl_K.squaredNorm();
          double exp_recl_term = exp(-K2 / (4 * eta)) / K2;
          if (exp_recl_term > precision / max_periods)
            rec_data.emplace_back(recl_K, exp_recl_term);
        }
      }
    }
  }
  std::sort(
      real_data_shifts.begin(), real_data_shifts.end(),
      [](const auto &a, const auto &b) -> bool { return a.norm() < b.norm(); });
}

double cryst_tools::ewald_sum::get_energy(const Eigen::Vector3d &vd) const {
  double vd_norm = vd.norm();
  double real_term = 0;
  bool self_iteraction = false;
  for (const auto &s : real_data_shifts) {
    // For sorted vectors no needs to continue all length will be higher
    if (s.norm() - vd_norm > r_max)
      break;

    double r = (s + vd).norm();
    if (r < 1E-3) {
      self_iteraction = true;
    } else if (r < r_max) {
      real_term += std::erfc(sqeta * r) / r;
    }
  }

  // K-space summation
  double rec_term = 0;
  for (const auto &cv : rec_data) {
    double Kvd = cv.recl_K.dot(vd);
    rec_term += std::cos(Kvd) * cv.exp_recl_term;
  }

  double result = 0.5 * real_term + 0.5 * 4 * cnst::pi / volume * rec_term;
  if(self_iteraction)
    result -= sqrt(eta / cnst::pi);
  
  return result;
}

Eigen::MatrixXd cryst_tools::ewald_sum::potential_matrix(const std::vector<Eigen::Vector3d> &vd) const
{
  Eigen::MatrixXd result;
  
  result.resize(vd.size(), vd.size());
  result.setZero();

  for(int i = 0; i < vd.size(); i++)
  {
    for(int j = i + 1; j < vd.size(); j++)
    {
      double eng = get_energy(vd[i] - vd[j]);
      result(i, j) = eng;
      result(j, i) = eng;
    }
  }
  double eng_zero = get_energy(Vector3d::Zero());
  for(int i = 0; i < vd.size(); i++) {
    result(i, i) = eng_zero;
  }
  return result;
}        
