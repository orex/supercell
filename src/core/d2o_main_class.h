/* 
 * File:   d2o_main_class.h
 * Author: kirill
 *
 * Created on July 19, 2013, 2:25 PM
 */

#ifndef D2O_MAIN_CLASS_H
#define	D2O_MAIN_CLASS_H

#include <string>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <random>

#include "common_types.h"
#include "permut_process_t.h"
#include "cif_io.h"

#include <boost/utility.hpp>
#include <boost/random.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/circular_buffer.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#ifdef LIBARCHIVE_ENABLED
#include <archive.h>
#endif    


struct c_occup_item
{
  std::string label;
  double occup_target;
  int num_of_atoms_sc;
  double charge;
  int el_num;

  c_occup_item() = delete;
  c_occup_item(const atom_t &a, double charge_v) : label(a.label), occup_target(a.occupancy),
                                                   num_of_atoms_sc(0), charge(charge_v), el_num(a.el_num) {};

  bool operator < ( const c_occup_item &r ) const
  { return ( label < r.label ); }
};

struct t_comb_descr {
  std::vector<int> prm_indexes;
  std::vector<int> group_indexes;
};

class c_occup_group
{
private:
  bool _fixed;  
  bool _fixed_fast;
  
public:
  std::vector<c_occup_item> items;
  std::vector<Eigen::Vector3d> positions;  //Cartesian
  std::vector<std::vector<int>> symms_sets;
  double max_dis_within_group;
  c_occup_group(): _fixed(false), _fixed_fast(false), max_dis_within_group(0)
  { ; };
  double get_total_occup_input() const;
  int    get_total_num_occup_sites() const;
  int64_t get_number_of_combinations() const;
  void set_fixation(bool fix = true)
  { _fixed = fix; };
  inline bool is_fixed() const
  { return _fixed || (get_number_of_combinations() == 1); }
  inline bool is_empty() const {
    int numa = 0;
    for(const auto &itm : items) {
      numa += itm.num_of_atoms_sc;
    }
    return numa == 0;
  }

  inline bool is_fixed_fast() const
  { return _fixed_fast; }
  inline void set_fixed_fast()
  { _fixed_fast = is_fixed(); }
  inline int number_of_sites() const
  { return positions.size(); };
  inline void add_item(const atom_t &oa, double charge) {
      items.emplace_back(oa, charge);
  }
};

struct site_charges
{
  double input_charge;
  double curr_charge;
  double occup;
  double cif_mult;
};

class struct_processor {
public:
  static const std::string coulomb_energy_suffix;
  struct_processor(const std::string &prefix_, int64_t tot_comb_) : prefix(prefix_),
  index_length(boost::lexical_cast<std::string>(tot_comb_).length()) {};
  std::string file_name(const struct_info &si, const std::string &sampl_type = "") const;
  std::string get_q_file_name(const std::string &suffix = "") const
  {
    return prefix + coulomb_energy_suffix + ( suffix.empty() ? std::string("") : std::string("_") ) + suffix + ".txt";
  }

  static std::string get_energy_str( const struct_info &si, int prec = 3) {
    std::stringstream result;
    result << std::fixed << std::setprecision(prec) << si.energy << " eV";
    return result.str();
  }

  static std::string get_energy_line(const std::string &file_name, const struct_info &si, int prec = 3) {
    return file_name + "\t" + get_energy_str(si, prec);
  }

  static std::string get_struct_title(const struct_info &si) {
    return std::string("Supercell generated structure. E_col = ") + get_energy_str(si, 10);
  }

private:
  std::string prefix;
  int index_length;
};

class c_struct_sel_containers : public c_struct_sel
{
 private:
   rnd_indexer_t rnd_indexer;
public:
  std::vector<struct_info> rnd_container;
  
  std::vector<struct_info> low_container;
  std::vector<struct_info> high_container;

  std::vector<struct_info> first_container;
  boost::circular_buffer<struct_info> last_container;

  std::vector<struct_info> weight_container;
  
public:
  c_struct_sel_containers(std::random_device::result_type seed) : rnd_indexer(seed) {};
  //symm_op 1 for non-merging
  bool sampling_active() const
  { return !save_all();  };
  void set_containers_prop(int64_t total_comb, int symm_op);
  void add_structure(const struct_info &si);
  void prepare_to_store();
};


struct q_energy_reduced {
 private:
   static constexpr double to_eV = 14.39964547842567205854;
 public:
  Eigen::MatrixXd qf;
  Eigen::VectorXd vf;
  double cf;
  std::vector<int> d_asn;
  std::vector<std::vector<double>> charges;
  inline double calc_energy(const t_vec_comb &vc) const {
    assert(vc.size() == d_asn.size() && qf.cols() == vc.size());
    Eigen::VectorXd v = Eigen::VectorXd::Zero(vc.size());
    for(int i = 0; i < vc.size(); i++)
      v[i] = charges[d_asn[i]][vc[i]];
    return to_eV * ( cf + 2 * v.dot(vf) + v.transpose() * qf * v);
  };
};

class d2o_main_class
{
private:
  struct rangi
  {
    int group_index, atom_index;
    int curr_value, min_value, max_value;
  };

private:
  int verbose_level;
  bool charge_balancing;

  cryst_structure_t orig_cst, supercell_cst;

  c_man_atom_prop * manual_properties;
  c_struct_sel_containers ss_p;

  std::vector<c_occup_group> occup_groups;

  double r_tolerance;
  double min_dist_between_groups;

  std::map<std::string, site_charges> scs;

  bool calc_q_energy;
  bool create_q_file;
  Eigen::MatrixXd q_energy;
  std::ofstream f_q_calc;

  #ifdef LIBARCHIVE_ENABLED
  archive *tar_container;
  #endif

  bool create_tar_container(const std::string &tar_fname);
  bool add_file_to_tar(const std::string &fname, const std::stringstream &strm);
  bool tar_enabled();
  bool close_tar_container();

  bool fix_groups();

  bool read_cryst_structure(std::string file_name);
  static std::pair<int, int> correct_rms_range(const int total_sites,
                                               const double occup,
                                               const double x2);

  int64_t total_combinations();
  double ss_charge_by_occup_groups();

  std::vector< rangi > get_rangi_array(const double x2);
  bool get_atoms_population();
  bool create_super_cell(int a, int b, int c);
  bool create_occup_groups();
  bool show_groups_information();
  bool process_charges(charge_balance cb);
  bool calculate_q_matrix();
  bool create_symmetry_operations_groups();
  static std::vector<int> index_symmetries(const Eigen::Matrix3d &cell,
                                    const Eigen::Affine3d &af,
                                    const std::vector<Eigen::Vector3d> &pos);
  template<class ConstIterator>
  void store_cont_cif(const ConstIterator& begin,
                      const ConstIterator& end,
                      const struct_processor &sp,
                      const t_comb_descr &cd,
                      const std::string &sampl_type = "")
  {
    for(ConstIterator it = begin; it != end; it++)
      write_struct(sp, *it, cd, sampl_type);
  }

  template<class ConstIterator>
  void store_cont_eng(const ConstIterator& begin,
                      const ConstIterator& end,
                      const struct_processor &sp,
                      const std::string &sampl_type = "")
  {
    if(begin == end)
      return;

    std::stringstream fq;
    std::string qfname = sp.get_q_file_name(sampl_type);
    for(ConstIterator it = begin; it != end; it++)
    {
      std::string fname_str = sp.file_name(*it, sampl_type);
      fq << struct_processor::get_energy_line(fname_str, *it) + "\n";
    }
    add_file_to_tar(qfname, fq);
  }

  bool write_struct(const struct_processor &sp, const struct_info &si, const t_comb_descr & cd,
                    const std::string &sampl_type = "");

  bool store_sampling(const std::string &output_base_name, const t_comb_descr &cd, int64_t tot_comb);
  bool write_files(const std::string &output_base_name, bool dry_run, bool merge_confs);
  std::tuple<t_symm_set, t_vec_comb, t_comb_descr> create_init_perm_structs() const;
  q_energy_reduced reduce_q_matrix(const t_comb_descr &cd) const;
  bool add_confs_to_mol(cif_output &co, const t_comb_descr &cd, const t_vec_comb &ppc) const;
  static std::string get_formula(const cryst_structure_t &cs);
  std::string get_formula_by_groups();
  bool set_labels_to_manual();
public:
  static inline double charge_tol() { return 1E-1; };
  static inline double occup_tol() { return 2E-3; };
  static inline double symm_tol() { return 1E-2; };
public:
  d2o_main_class() = delete;
  d2o_main_class(std::random_device::result_type seed) : ss_p(seed), charge_balancing(false) {};
  void set_verbosity(int vb)
  { verbose_level = vb; };
  bool process(std::string input_file_name, bool dry_run,
               const std::vector<int> &supercell_shape,
               charge_balance cb, double tolerance_v,
               bool merge_confs, bool calc_q_energy_v, bool create_q_file_v,
               c_man_atom_prop &manual_properties,
               const c_struct_sel &ss_p,
               std::string output_base_name,
               std::string output_tar_name);
  bool check_properties_consistency() const;
};

#endif	/* D2O_MAIN_CLASS_H */

