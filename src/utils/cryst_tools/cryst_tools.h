/* 
 * File:   cryst_tools.h
 * Author: kirill
 *
 * Created on December 16, 2013, 4:00 PM
 */

#ifndef CRYST_TOOLS_H
#define	CRYST_TOOLS_H

#include <vector>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>

typedef std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d> > vector_Affine3d;


namespace cryst_tools {
  typedef std::vector< std::vector<Eigen::Vector3d> > vc_sets;

  inline Eigen::Vector3d min_frac(const Eigen::Vector3d &frac) {
    Eigen::Vector3d result;

    for(int i = 0; i < 3; i++) {
      if( std::abs(frac[i]) < 0.5 )
        result[i] = frac[i];
      else
        result[i] = std::remainder(frac[i], 1.0);
    }

    return result;
  }
  Eigen::Vector3d norm_frac(const Eigen::Vector3d &frac);

  std::vector<Eigen::Matrix3d> get_cell_symmetries(const Eigen::Matrix3d &cell, const double tol);

  std::vector<Eigen::Vector3d> get_shifts(const Eigen::Matrix3d &sym_matrix_B,
          const Eigen::Matrix3d &cell,
          const std::vector<Eigen::Vector3d> &frac_coords,
          double tolerance, bool all_symmetric);

  std::vector<Eigen::Vector3d> shifts_intersect(const Eigen::Matrix3d &cell,
          const std::vector<Eigen::Vector3d> &frac_c1,
          const std::vector<Eigen::Vector3d> &frac_c2,
          double tolerance);

  vector_Affine3d get_all_symmetries(const Eigen::Matrix3d &cell,
          const vc_sets &frac_coords,
          const double tol);

  vector_Affine3d get_all_symmetries(const Eigen::Matrix3d &cell,
          const vc_sets &frac_coords,
          const std::vector<bool> &all_symm,
          const double tol);

};

//Minimal distance utilites
namespace cryst_tools {

  class min_dist {
  protected:
    double r2_dist_direct;
    Eigen::Vector3i box_size;
    Eigen::Matrix3d cell;
    Eigen::Matrix3d r_cell;
    bool cell_assigned;
  public:
    min_dist() : cell_assigned(false) {};
    void set_cell(const Eigen::Matrix3d &cell_v, double tol = 1E-3);
    Eigen::Vector3d operator ()(const Eigen::Vector3d &vd) const;
    std::vector<Eigen::Vector3d> get_img_dist(const Eigen::Vector3d &vd, 
                                              const Eigen::Vector3i &box) const;
    Eigen::Vector3d average_vector(const std::vector<Eigen::Vector3d> &av) const;
  };
  
};

//Ewald summation utilites
namespace cryst_tools {

  class ewald_sum {
   private:
    static constexpr int max_cell_check_range = 20;
    struct recep_data_t {
      Eigen::Vector3d recl_K;
      double exp_recl_term;
      recep_data_t(Eigen::Vector3d recl_K_v, double exp_recl_term_v)
          : recl_K(recl_K_v), exp_recl_term(exp_recl_term_v) {};
    };

  private:
    Eigen::Matrix3d cell;
    Eigen::Matrix3d res_cell;
    double volume;
    
    double precision;
    double r_max;
    double eta;
    double sqeta;

    std::vector<Eigen::Vector3d> real_data_shifts;
    std::vector<recep_data_t> rec_data;

  public:
    ewald_sum()
    { precision = -1; };
    void set_cell(const Eigen::Matrix3d &cell_v);
    void set_precision(double N, double precision_v = 1E-8);
    double get_energy(const Eigen::Vector3d &vd) const;
    Eigen::MatrixXd potential_matrix(const std::vector<Eigen::Vector3d> &vd) const;
  };
  
};



#endif	/* CRYST_TOOLS_H */

