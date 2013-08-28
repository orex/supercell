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

#include <openbabel/mol.h>

#include "common_types.h"

#include "boost/utility.hpp"

struct c_occup_item
{
  std::string label;
  double occup_target;
  int num_of_atoms_sc;
  double charge;
  OpenBabel::OBAtom * obp;
  
  c_occup_item()
  { obp = new OpenBabel::OBAtom(); }
  
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


class c_occup_group
{
protected:
  bool _fixed;  
public:
  std::vector<c_occup_item> items;
  double get_total_occup_input() const;
  int    get_total_num_occup_sites() const;
  int64_t get_number_of_combinations() const;
  void set_fixation(bool fix = true)
  { _fixed = fix; };
  bool is_fixed() const 
  { return _fixed || (get_number_of_combinations() == 1); }
  double max_dis_within_group;
  int number_of_sites;
};

struct site_charges
{
  double input_charge;
  double curr_charge;
  double occup;
  double cif_mult;
};        

class lbl_order
{
public:
  std::string lbl1;
  std::string lbl2;

  double tolerance;
  int64_t mult2;
 
  lbl_order(std::string l1, std::string l2, double tol, int64_t m2) : 
            tolerance(tol), mult2(m2)
  {
    lbl1 = min(l1, l2);
    lbl2 = max(l1, l2);
  };

  bool operator < ( const lbl_order &r ) const
  { return ( mult2 < r.mult2 ); }
};

class map_comp_item : private boost::noncopyable
{
protected:
  static int compare_distances(const std::vector<float> &f1, const std::vector<float> &f2, float tol);  
protected:
  OpenBabel::OBUnitCell * uc;
  OpenBabel::OBUnitCell * unitcell();
  
  std::vector<lbl_order> * labels_order;
  std::map<std::string, std::vector<float> > * dsts;
  //output is sorted
  std::vector<float> get_lengths_by_labels(const lbl_order &lbl);
  int mult;
public:  
  map_comp_item(std::vector<lbl_order> &lbl);
  OpenBabel::OBMol *mol;
  
  void inc_mult()
  { mult++; };
  
  int get_mult()
  { return mult; };
  
  static int comp( map_comp_item &r1, map_comp_item &r2 );
  ~map_comp_item() { delete dsts; delete mol; };
};

class d2o_main_class 
{
public: //Static
  static OpenBabel::vector3 get_minimal_distance(OpenBabel::vector3 dist, 
                                                 OpenBabel::OBUnitCell * unitcell);
  
  static OpenBabel::vector3 center_mass(const std::vector<OpenBabel::vector3> &atoms_pos,
                                        OpenBabel::OBUnitCell * unitcell,
                                        const double tol);
  
protected:
  struct rangi
  {
    int group_index, atom_index;
    int curr_value, min_value, max_value;
  };
  
protected:
  int verbose_level;
  
  OpenBabel::OBMol mol_initial;
  OpenBabel::OBMol mol_supercell;  
  
  c_man_atom_prop * manual_properties;

  std::vector<c_occup_group> occup_groups;

  double r_tolerance;
  double min_dist_between_groups;

  std::map<std::string, site_charges> scs;

  bool fix_groups();
  
  bool read_molecule(std::string file_name);
  void correct_rms_range(const int total_sites, 
                         const double occup, 
                         const double x2,
                         int &min_value,
                         int &max_value);

  int64_t total_combinations();
  
  bool add_to_list(std::list<map_comp_item *> &mpis, map_comp_item * mpi);  
  std::vector< rangi > get_rangi_array(const double x2);
  std::vector<lbl_order> set_lbl_order();
  bool get_atoms_population();
  bool create_super_cell(int a, int b, int c);
  bool create_occup_groups();
  bool show_groups_information();
  bool process_charges(charge_balance cb, bool verbose = false);
  bool write_files(std::string output_base_name, double n_store, bool dry_run, bool merge_confs);
  bool init_atom_change_mol(OpenBabel::OBMol *cmol);
  bool add_confs_to_mol(OpenBabel::OBMol *cmol, const std::map<int, std::vector<int> > &ppc);
  std::string get_formula(OpenBabel::OBMol &mol);
  std::string get_formula_by_groups();
  bool set_labels_to_manual();
public:
  static const double charge_tol = 1E-1;
  static const double occup_tol = 2E-3;
public:
  d2o_main_class();
  void set_verbosity(int vb)
  { verbose_level = vb; };
  bool process(std::string input_file_name, bool dry_run,
               const std::vector<int> scs,
               charge_balance cb, double tolerance_v, 
               bool merge_confs,
               c_man_atom_prop &manual_properties,
               double n_store,                       
               std::string output_base_name);
  virtual ~d2o_main_class();
};

#endif	/* D2O_MAIN_CLASS_H */

