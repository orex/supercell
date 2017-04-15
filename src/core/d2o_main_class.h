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

#include <openbabel/mol.h>

#include "common_types.h"

#include <boost/utility.hpp>
#include <boost/random.hpp>
#include <boost/format.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "obabel/obabel_utils.h"

#ifdef LIBARCHIVE_ENABLED
#include <archive.h>
#endif    


struct c_occup_item
{
  std::string label;
  double occup_target;
  int num_of_atoms_sc;
  double charge;
  OpenBabel::OBAtom * obp;
  
  c_occup_item(OpenBabel::OBAtom *ob, double charge_v);

  c_occup_item(const c_occup_item &orig)
  { 
    label = orig.label;
    occup_target = orig.occup_target;
    num_of_atoms_sc = orig.num_of_atoms_sc;
    charge = orig.charge;
    obp = new OpenBabel::OBAtom();
    obp->Duplicate(orig.obp);
  }

  ~c_occup_item()
  { delete obp; }
  
  
  bool operator < ( const c_occup_item &r ) const
  { return ( label < r.label ); }
};

typedef std::vector<int> symm_set;

class c_occup_group
{
protected:
  bool _fixed;  
  bool _fixed_fast;
  
public:
  std::vector<c_occup_item> items;
  std::vector<OpenBabel::vector3> positions;  //Cartesian
  std::vector<symm_set> symms_sets;
  c_occup_group(): _fixed(false), _fixed_fast(false)
  { ; };
  double get_total_occup_input() const;
  int    get_total_num_occup_sites() const;
  int64_t get_number_of_combinations() const;
  void set_fixation(bool fix = true)
  { _fixed = fix; };
  bool is_fixed() const 
  { return _fixed || (get_number_of_combinations() == 1); }
  bool is_fixed_fast() const 
  { return _fixed_fast; }
  void set_fixed_fast()
  { _fixed_fast = is_fixed(); }
  double max_dis_within_group;
  int number_of_sites() const
  { return positions.size(); };
  void add_item(OpenBabel::OBAtom * oba, double charge);
};

struct site_charges
{
  double input_charge;
  double curr_charge;
  double occup;
  double cif_mult;
};

typedef std::vector< std::vector<int> > t_vec_comb;

class struct_info
{
public:
  t_vec_comb cmb;
  double energy;
  int index;
  //weight should be equal 0 for non-merge run.
  int weight;
  
  struct_info(): energy(0), index(0), weight(0) {};

  std::string file_name(const std::string &prefix, int tot_comb, 
                        const std::string &sampl_type = "") const;
  
  bool operator < (const struct_info &second ) const
  {
    return (energy < second.energy) || ( (energy == second.energy) && (index < second.index) );
  }
};

class c_struct_sel_containers : public c_struct_sel
{
protected:  
  boost::mt19937 rnd;
  int min_comb;
  int symm_op;
  double probability;
public:
  std::list<struct_info> rnd_container;
  
  std::set<struct_info> low_container;
  std::set<struct_info> high_container;

  std::vector<struct_info> first_container;
  std::list<struct_info> last_container;
  
public:  
  //symm_op 1 for non-merging
  bool sampling_active() const
  { return !save_all();  };
  void set_containers_prop(int total_comb, int symm_op_v);  
  void add_structure(const struct_info &si, int str_left);
  void prepare_to_store();
};

class d2o_main_class 
{
protected:
  struct rangi
  {
    int group_index, atom_index;
    int curr_value, min_value, max_value;
  };
  
protected:
  int verbose_level;
  bool charge_balancing;
  
  OpenBabel::OBMol mol_initial;
  OpenBabel::OBMol mol_supercell;
  
  c_man_atom_prop * manual_properties;
  c_struct_sel_containers ss_p;

  std::vector<c_occup_group> occup_groups;

  double r_tolerance;
  double min_dist_between_groups;

  std::map<std::string, site_charges> scs;
  
  bool calc_q_energy;
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
  
  bool read_molecule(std::string file_name);
  void correct_rms_range(const int total_sites, 
                         const double occup, 
                         const double x2,
                         int &min_value,
                         int &max_value);

  int64_t total_combinations();
  double ss_charge_by_occup_groups();
  
  std::vector< rangi > get_rangi_array(const double x2);
  bool get_atoms_population();
  bool create_super_cell(int a, int b, int c);
  void FillUnitCell_rmdup(OpenBabel::OBMol * mol);
  bool create_occup_groups();
  bool show_groups_information();
  bool process_charges(charge_balance cb);
  bool calculate_q_matrix();
  double calculate_q_energy(const t_vec_comb &mc);
  bool create_symmetry_operations_groups();
  std::vector<int> index_symmetries( OpenBabel::OBUnitCell * uc, 
                                     const Eigen::Affine3d &af, 
                                     const std::vector<OpenBabel::vector3> &pos);
  template<class ConstIterator>
  void store_cont_cif(const ConstIterator& begin, 
                       const ConstIterator& end,
                       const std::string &prefix, int tot_comb, 
                       const std::string &sampl_type = "")
  {
    for(ConstIterator it = begin; it != end; it++)
      write_struct(*it, prefix, tot_comb, sampl_type);
  }

  template<class ConstIterator>
  void store_cont_eng(const ConstIterator& begin, 
                      const ConstIterator& end,
                      const std::string &prefix, int tot_comb, 
                      const std::string &sampl_type = "")
  {
    if(begin == end)
      return;
    
    std::stringstream fq;
    std::string qfname = get_q_file_name(prefix, sampl_type);
    for(ConstIterator it = begin; it != end; it++)
    {  
      std::string fname_str = it->file_name(prefix, tot_comb, sampl_type);
      fq << boost::format("%1%\t%2$.3f eV\n") % fname_str % it->energy;
    }  
    add_file_to_tar(qfname, fq);
  }
  
  bool write_struct(const struct_info &si,
                     const std::string &prefix, int tot_comb, 
                     const std::string &sampl_type = "");
  
  bool store_samling(std::string output_base_name, int tot_comb);
  bool write_files(std::string output_base_name, bool dry_run, bool merge_confs);
  bool check_comb_unique(const t_vec_comb &mc, int &merged_comb);
  bool create_comb(const symm_set &sc, const std::vector<int> &cmb_in, std::vector<int> &cmb_out);
  bool init_atom_change_mol(OpenBabel::OBMol *cmol, const struct_info &si);
  bool add_confs_to_mol(OpenBabel::OBMol *cmol, const t_vec_comb &ppc);
  std::string get_formula(OpenBabel::OBMol &mol);
  std::string get_formula_by_groups();
  bool set_labels_to_manual();
  std::string get_q_file_name(const std::string &output_base_name, const std::string &suffix);
public:
  static const double charge_tol = 1E-1;
  static const double occup_tol = 2E-3;
  static const double symm_tol = 1E-2;
public:
  d2o_main_class();
  void set_verbosity(int vb)
  { verbose_level = vb; };
  bool process(std::string input_file_name, bool dry_run,
               const std::vector<int> &scs,
               charge_balance cb, double tolerance_v, 
               bool merge_confs, bool calc_q_energy_v,
               c_man_atom_prop &manual_properties,
               const c_struct_sel &ss_p,
               std::string output_base_name,
               std::string output_tar_name);
  virtual ~d2o_main_class();
};

#endif	/* D2O_MAIN_CLASS_H */

