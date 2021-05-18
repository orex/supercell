/* 
 * File:   d2o_main_class.cpp
 * Author: kirill
 * 
 * Created on July 19, 2013, 2:25 PM
 */

#include "d2o_main_class.h"

#include <iostream>
#include <iomanip>
#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>
#include <boost/algorithm/string.hpp>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <gemmi/elem.hpp>
#define PT_GETSYMBOL(atomic_num) gemmi::Element(atomic_num).name()

#include "science/combinatorics.h"

#include "others/rnd_utils.h"
#include "others/string_utils.h"
#include "containers/array_common.hpp"
#include "containers/hash_unique.h"
#include "permut_process_t.h"


#include "cryst_tools/comb_points.h"
#include "cryst_tools/cryst_tools.h"

#include <Eigen/Dense>
#include <tbb/pipeline.h>
#include <tbb/concurrent_queue.h>

#ifdef LIBARCHIVE_ENABLED
#include <archive_entry.h>
#endif


using namespace std;

#if defined(LIBARCHIVE_ENABLED) && (!defined(LIBARCHIVE_PATCH_DISABLE))
#include "archive_extent_patch.h"
#endif

bool d2o_main_class::create_tar_container(const std::string &fname)
{
  bool result = true;
  #ifdef LIBARCHIVE_ENABLED
  if( ! fname.empty() )
  {  
    tar_container = archive_write_new();
    result = archive_write_set_format_filter_by_ext(tar_container, fname.c_str()) == ARCHIVE_OK;
    if( result )
      result = archive_write_open_filename(tar_container, fname.c_str()) == ARCHIVE_OK;
  }
  else
    tar_container = nullptr;
  #endif
  return result;
}

bool d2o_main_class::add_file_to_tar(const std::string &fname, const std::stringstream &strm)
{
  #ifdef LIBARCHIVE_ENABLED    
  if( tar_container != nullptr )
  {  
    archive_entry *entry = archive_entry_new();
    archive_entry_set_pathname(entry, fname.c_str());
    archive_entry_set_size(entry, strm.str().size());
    archive_entry_set_filetype(entry, AE_IFREG);
    archive_entry_set_perm(entry, 0644);
    archive_write_header(tar_container, entry);
    archive_write_data(tar_container, strm.str().c_str(), strm.str().size());
    archive_entry_free(entry);
  }
  else
  {
  #endif  
    std::ofstream out_file( fname.c_str() );
    out_file << strm.rdbuf();
    out_file.close();
  #ifdef LIBARCHIVE_ENABLED    
  }  
  #endif
  return true;
}

bool d2o_main_class::tar_enabled()
{
  #ifdef LIBARCHIVE_ENABLED
  return tar_container != nullptr;
  #else
  return false;
  #endif 
}

bool d2o_main_class::close_tar_container()
{
  #ifdef LIBARCHIVE_ENABLED
  if( tar_container != nullptr)
  {  
    bool b1, b2;
    b1 = archive_write_close(tar_container) == ARCHIVE_OK;
    b2 = archive_write_free(tar_container) == ARCHIVE_OK;
    return b1 && b2;
  }  
  else    
    return true;
  #else
  return true;
  #endif
}


std::string struct_processor::file_name(const struct_info &si, const std::string &sampl_type) const
{
  stringstream result;

  result << prefix + "_i" + sampl_type << std::setfill('0') << std::setw(index_length) << std::internal << si.index;

  if( si.weight > 0 )
    result << "_w" << si.weight;
  
  result << ".cif";
  
  return result.str();
}

void c_struct_sel_containers::add_structure(const struct_info &si)
{
  //Random
  if(str_random_count() > 0)
  {  
    if( rnd() <= probability )
      rnd_container.push_back(si);
  }
  
  //First  
  if( first_container.size() < str_first_count() )
    first_container.push_back(si);
  
  //Last
  if( str_last_count() > 0 ) {
    last_container.push_back(si);
  }  
  
  //Low
  if(str_low_count() > 0) {
    if (low_container.size() < str_low_count()) {
      low_container.push_back(si);
      std::push_heap(low_container.begin(), low_container.end());
    } else {
      if (si < low_container.front()) {
        low_container.push_back(si);
        std::push_heap(low_container.begin(), low_container.end());
        std::pop_heap(low_container.begin(), low_container.end());
        low_container.pop_back();
      }
    }  
  }  
  
  //High
  if(str_high_count() > 0) {
    auto cmp = [](const auto &a, const auto &b) -> bool {
      return b < a;
    };

    if (high_container.size() < str_high_count()) {
      high_container.push_back(si);
      std::push_heap(high_container.begin(), high_container.end(), cmp);
    } else {
      if ( high_container.front() < si ) {
        high_container.push_back(si);
        std::push_heap(high_container.begin(), high_container.end(), cmp);
        std::pop_heap(high_container.begin(), high_container.end(), cmp);
        high_container.pop_back();
      }
    }  
  }  
   
  //Weight
  if( si.weight <= str_weight_limit() )
    weight_container.push_back(si);
}

void c_struct_sel_containers::prepare_to_store()
{
  assert(str_first_count() >= first_container.size());
  assert(str_last_count()  >= last_container.size());
  assert(str_high_count()  >= high_container.size());
  assert(str_low_count()   >= low_container.size());  
  
  random_thin_to(rnd_container, str_random_count(), rnd);
  std::sort(low_container.begin(), low_container.end());
  std::sort(high_container.begin(), high_container.end());

}


void c_struct_sel_containers::set_containers_prop(int64_t total_comb, int symm_op_v)
{
  symm_op = symm_op_v;
  min_comb = max<int64_t>(total_comb / symm_op, 1);
  
  double real_prob = 4 * double(str_random_count()) / double(min_comb);
  if( real_prob >= 1.0 ) {
    probability = rnd_engine_t::max();
  } else {
    probability = static_cast<std::uint64_t>( real_prob * rnd_engine_t::max());
  }
  
  rnd_container.clear();
  first_container.clear();
  last_container.clear();
  high_container.clear();
  low_container.clear();
  weight_container.clear();
  
  rnd_container.reserve(5 * str_random_count());
  first_container.reserve(str_first_count());
  last_container.set_capacity(str_last_count());
  high_container.reserve(str_high_count() + 1);
  low_container.reserve(str_low_count() + 1);

}

double c_occup_group::get_total_occup_input() const
{
  double result = 0.0;
  
  for(std::vector<c_occup_item>::const_iterator it = items.begin(); it != items.end(); ++it)
    result += (*it).occup_target;
  
  return result;
}

int c_occup_group::get_total_num_occup_sites() const
{
  int result = 0;
  
  for(std::vector<c_occup_item>::const_iterator it = items.begin(); it != items.end(); ++it)
    result += (*it).num_of_atoms_sc;
  
  return result;
}

int64_t c_occup_group::get_number_of_combinations() const
{
  int64_t result = 1;

  if( !_fixed )
  {
    std::map<int,int> nm;

    int sum = 0;  
    for(int i = 0; i < items.size(); i++)
    {  
      nm[i] = items[i].num_of_atoms_sc;
      sum  += items[i].num_of_atoms_sc;
    }  

    if(sum <= number_of_sites())
    {  
      nm[-1] = number_of_sites() - sum;
    
      bool good = num_combinations(nm, result);
      if(!good)
        result = -1;
    }
    else
      result = 0;
  }
  
  return result;
}

std::string d2o_main_class::get_formula_by_groups()
{
  map<string, int> formula_map;
  
  for(int i = 0; i < occup_groups.size(); i++)
  {
    for(int j = 0; j < occup_groups[i].items.size(); j++)
    {
      const c_occup_item &ci = occup_groups[i].items[j]; 
      string atom_symbol = PT_GETSYMBOL(ci.el_num);
      formula_map[atom_symbol] += ci.num_of_atoms_sc;
    }  
  }
  
  string result;
  
  for(map<string, int>::iterator it = formula_map.begin();
                                 it != formula_map.end(); ++it)
  {
    boost::format fmter("%1%%2% ");
    fmter % it->first % it->second;
    result += fmter.str();
  }  
  
  trim(result);
  
  return result;      
}

std::string d2o_main_class::get_formula(const cryst_structure_t &cs)
{
  map<string, double> formula_map;
  
  for(const auto &a : cs.atoms) {
    string atom_symbol = PT_GETSYMBOL(a.el_num);
    formula_map[atom_symbol] += a.occupancy;
  }
  
  string result;
  
  for(map<string, double>::iterator it = formula_map.begin();
                                    it != formula_map.end(); ++it)
  {
    boost::format fmter("%1%%2% ");
    fmter % it->first % it->second;
    result += fmter.str();
  }  
  
  trim(result);
  
  return result;      
}

bool d2o_main_class::add_confs_to_mol(cif_output &co, const t_comb_descr &cd, const t_vec_comb &ppc) const {
  for (int i = 0; i < occup_groups.size(); i++) {
    const c_occup_group &curr_group = occup_groups[i];
    int gpos = std::distance(cd.group_indexes.cbegin(), std::find(cd.group_indexes.cbegin(), cd.group_indexes.cend(), i));
    bool fixed_group = gpos == cd.group_indexes.size();
    assert( fixed_group || cd.prm_indexes[gpos + 1] - cd.prm_indexes[gpos] == curr_group.positions.size() );
    int cinx = cd.prm_indexes[gpos];
    for (int j = 0; j < curr_group.positions.size(); j++) {
      for (int k = 0; k < curr_group.items.size(); k++) {
        if (fixed_group ||  ppc[j + cinx] == k + 1 ) {
          double occup_value = fixed_group ? double(curr_group.items[k].num_of_atoms_sc) /
                double(curr_group.number_of_sites()) : 1.0;
	
          co.add_atom(curr_group.items[k].el_num, curr_group.items[k].label,
                      curr_group.positions[j], occup_value);
            }  
	  }
        }
      }

  return true;
}

bool d2o_main_class::calculate_q_matrix()
{
  using namespace Eigen;
  using namespace cryst_tools;
  
  vector<Vector3d> all_pos;
  for(auto & occup_group : occup_groups)
  {
    for(int j = 0; j < occup_group.positions.size(); j++)
      all_pos.push_back(occup_group.positions[j]);
  }
  
  if(verbose_level >= 2)
    cout << "Start Coulomb matrix (" << all_pos.size() << "x" << all_pos.size() << ") calculation." << endl;

  ewald_sum es;
  
  es.set_cell(supercell_cst.unit_cell.cell());
  es.set_precision(all_pos.size(), 1E-7);
  
  q_energy = es.potential_matrix(all_pos);
 
  if(verbose_level >= 2)
    cout << "Coulomb matrix calculation finished." << endl;
  
  return true;
}

bool d2o_main_class::create_symmetry_operations_groups()
{
  using namespace cryst_tools;
  using namespace Eigen;

  Eigen::Matrix3d m_cell = supercell_cst.unit_cell.cell().inverse();

  vc_sets vc;
  vc.resize(occup_groups.size());
  
  vector<bool> bc;
  bc.resize(occup_groups.size());
  for(int i = 0; i < occup_groups.size(); i++)
  {
    vc[i].resize(occup_groups[i].positions.size());
    bc[i] = occup_groups[i].is_fixed();
    for(int j = 0; j < occup_groups[i].positions.size(); j++)
    {  
      vc[i][j] = m_cell * occup_groups[i].positions[j];
    }
  }
  
  vector_Affine3d syms;
  
  syms = get_all_symmetries(supercell_cst.unit_cell.cell(), vc, bc, symm_tol());
  
  if(verbose_level >= 1)
    cout << syms.size() << " symmetry operation found for supercell." << endl;
  
  bool good_set = true;
  
  for(int i = 0; i < occup_groups.size(); i++)    
  {
    occup_groups[i].symms_sets.resize(syms.size());
    for(int j = 0; j < syms.size(); j++)
    {  
      occup_groups[i].symms_sets[j] = index_symmetries(supercell_cst.unit_cell.cell(), syms[j],
                                      occup_groups[i].positions);
      if( occup_groups[i].is_fixed() )
      {  
        good_set = find(occup_groups[i].symms_sets[j].begin(), 
                        occup_groups[i].symms_sets[j].end(), -1) == 
                        occup_groups[i].symms_sets[j].end();
      }          
      if(!good_set)
        break;
    }
    if(!good_set)
      break;
  }
  
  if(verbose_level >= 2)
    cout << "Symmetries operation assigned to groups." << endl;
  
  return good_set;
}

std::vector<int> d2o_main_class::index_symmetries(const Eigen::Matrix3d &cell,
                                                  const Eigen::Affine3d &af, 
                                                  const vector<Eigen::Vector3d> &pos)
{
  using namespace Eigen;
  
  std::vector<int> result;
  Matrix3d r_cell = cell.inverse();
  Affine3d cart_tr, cl_d, cl_i;
  
  cl_d.setIdentity();
  cl_d.linear() = cell;
  
  cl_i.setIdentity();
  cl_i.linear() = cell.inverse();
  
  cart_tr = cl_d * af * cl_i;
  
  std::vector<Vector3d> pose;  
  for(int i = 0; i < pos.size(); i++)
    pose.push_back(pos[i]);
  
  result.resize(pos.size(), -1);
  
  for(int i = 0; i < pose.size(); i++)
  {
    int index_sym = -1;
    int index_count = 0;
    for(int j = 0; j < pose.size(); j++)
    {
      Vector3d ve = pose[j] - cart_tr * pose[i];
      ve = cell * cryst_tools::min_frac(r_cell * ve);
      if(ve.norm() < symm_tol())
      {
        index_sym = j;
        index_count++;
      }  
    }
    assert(index_count <= 1);
    result[i] = index_sym;
  }  
  
  return result;
}

std::tuple<t_symm_set, t_vec_comb, t_comb_descr> d2o_main_class::create_init_perm_structs() const {
  t_comb_descr psm;
  t_vec_comb init_cmb;
  psm.prm_indexes.emplace_back(0);
  for(int i = 0; i < occup_groups.size(); i++) {
    if( !occup_groups[i].is_fixed() ) {
      map<base_prm_t, int> mvc;
  
      for(int j = 0; j < occup_groups[i].items.size(); j++)
        mvc[j + 1] = occup_groups[i].items[j].num_of_atoms_sc;
  
      mvc[0] = occup_groups[i].number_of_sites() -
          occup_groups[i].get_total_num_occup_sites();

      vector<base_prm_t> vp = create_start_combination(mvc);
      psm.group_indexes.emplace_back(i);
      psm.prm_indexes.emplace_back(psm.prm_indexes.back() + vp.size());
      for(int j = 0; j < vp.size(); j++)
        init_cmb.emplace_back(vp[j]);
    }
  }

  t_symm_set sm(occup_groups[0].symms_sets.size(), init_cmb.size());
  if( !occup_groups[0].symms_sets.empty() ) {
    int smpos = 0;
    for (int i = 0; i < occup_groups.size(); i++) {
      if (std::find(psm.group_indexes.cbegin(), psm.group_indexes.cend(), i) != psm.group_indexes.cend()) {
        for (int j = 0; j < occup_groups[i].symms_sets.size(); j++) {
          for (int k = 0; k < occup_groups[i].symms_sets[j].size(); k++) {
            sm.at(j, k + smpos) = occup_groups[i].symms_sets[j][k] + smpos;
          }
        }
        assert(occup_groups[i].symms_sets[0].size() == occup_groups[i].positions.size());
        smpos += occup_groups[i].symms_sets[0].size();
      }
    }
    assert(smpos == init_cmb.size());
  }
  return std::tuple<t_symm_set, t_vec_comb, t_comb_descr>(sm, init_cmb, psm);
}

q_energy_reduced d2o_main_class::reduce_q_matrix(const t_comb_descr &cd) const {
  using namespace Eigen;
  q_energy_reduced result;

  VectorXd q_v;
  q_v.resize(q_energy.cols());
  q_v.setZero();
  vector<int> dw;
  
  int q_v_pos = 0;
  int grp_id = 0;
  for(int i = 0; i < occup_groups.size(); i++) {
    const c_occup_group &curr_group = occup_groups[i];
    bool fixed_group = std::find(cd.group_indexes.cbegin(), cd.group_indexes.cend(), i) == cd.group_indexes.cend();
    if( !fixed_group ) {
      result.charges.emplace_back(curr_group.items.size() + 1, 0);
      for (int k = 0; k < curr_group.items.size(); k++)
        result.charges.back()[k + 1] = curr_group.items[k].charge;
    }
    for(int j = 0; j < curr_group.positions.size(); j++)
    {  
      double charge = 0;
      if( fixed_group ) {
        for (int k = 0; k < curr_group.items.size(); k++) {
           double occup_value = double(curr_group.items[k].num_of_atoms_sc) /
                  double(curr_group.number_of_sites());
          charge += occup_value * curr_group.items[k].charge;
      }
    } 
      q_v[q_v_pos] = charge;
      if( !fixed_group ) {
        dw.emplace_back(q_v_pos);
        result.d_asn.emplace_back(grp_id);
  }
      q_v_pos++;
    }
    if( !fixed_group )
      grp_id++;
  }
  assert(q_v_pos == q_energy.cols());
  assert(dw.size() == cd.prm_indexes.back());
  assert(grp_id == cd.group_indexes.size());
  
  result.cf = q_v.transpose() * q_energy * q_v;

  VectorXd m1 = q_energy * q_v;
  result.vf.resize(dw.size());
  for(int i = 0; i < dw.size(); i++)
    result.vf[i] = m1[dw[i]];

  result.qf.resize(dw.size(), dw.size());
  for(int i = 0; i < dw.size(); i++) {
    for(int j = 0; j < dw.size(); j++) {
      result.qf(i, j) = q_energy(dw[i], dw[j]);
    }
  }
  
  //0.5 not to count twice pairs
  //11.4 - to eV

  return result;
}

bool d2o_main_class::write_struct(const struct_processor &sp, const struct_info &si, const t_comb_descr & cd,
                                  const std::string &sampl_type)
{
  bool result = true;
  
  stringstream ss;
  cif_output co(ss, supercell_cst,
                std::string("E_e = ") + struct_processor::get_energy_str(si, 10));
  
  add_confs_to_mol(co, cd, si.cmb);

  string f_name = sp.file_name(si, sampl_type);
  
      add_file_to_tar(f_name, ss);
  
  if( !result )
    cerr << "An error occurred during storing of \"" << f_name << "\" file." << endl;

  return result;
}

bool d2o_main_class::store_sampling(const string &output_base_name, const t_comb_descr &cd, int64_t tot_comb)
{
  ss_p.prepare_to_store();

  struct_processor str_proc(output_base_name, tot_comb);
  
  store_cont_cif(ss_p.first_container.begin(),
                 ss_p.first_container.end(),
                 str_proc, cd, "f");

  store_cont_cif(ss_p.last_container.begin(),
                 ss_p.last_container.end(),
                 str_proc, cd, "a");

  store_cont_cif(ss_p.low_container.begin(),
                 ss_p.low_container.end(),
                 str_proc, cd, "l");

  store_cont_cif(ss_p.high_container.begin(),
                 ss_p.high_container.end(),
                 str_proc, cd, "h");

  store_cont_cif(ss_p.rnd_container.begin(),
                 ss_p.rnd_container.end(),
                 str_proc, cd, "r");

  store_cont_cif(ss_p.weight_container.begin(),
                 ss_p.weight_container.end(),
                 str_proc, cd, "w");
  
  if( calc_q_energy )
  {
    store_cont_eng(ss_p.first_container.begin(),
                   ss_p.first_container.end(),
                   str_proc, "f");

    store_cont_eng(ss_p.last_container.begin(),
                   ss_p.last_container.end(),
                   str_proc, "a");

    store_cont_eng(ss_p.low_container.begin(),
                   ss_p.low_container.end(),
                   str_proc, "l");

    store_cont_eng(ss_p.high_container.begin(),
                   ss_p.high_container.end(),
                   str_proc, "h");

    store_cont_eng(ss_p.rnd_container.begin(),
                   ss_p.rnd_container.end(),
                   str_proc, "r");

    store_cont_eng(ss_p.weight_container.begin(),
                   ss_p.weight_container.end(),
                   str_proc, "w");

  }  
    
  return true;
}

class generate_comb_t {
 private:
  tbb::concurrent_queue<permut_process_t *> prc_queue;
  const vector<int> &permi;
  t_vec_comb init_cmb, cur_combs;
  int64_t _combination_left;
  int64_t _total_index;
  int64_t _index;
  int packet_size;
  bool finish;
 public:
  generate_comb_t(int nthreads, int64_t __total_combinations, int64_t __packet_size,
                  const t_symm_set &sms, const vector<int> &__permi,
                  const t_vec_comb &init) : permi(__permi),
                                            init_cmb(init),
                                            cur_combs(init),
                                            _combination_left(__total_combinations),
                                            _total_index(0),
                                            _index(0),
                                            packet_size(__packet_size),
                                            finish(false) {
    for (int i = 0; i < nthreads; i++) {
      auto pp = new permut_process_t(sms, permi, packet_size);
      prc_queue.emplace(pp);
    }
  };
  
  inline int64_t combination_left() const {
    return _combination_left;
  }
  inline int64_t total_index() const {
    return  _total_index;
  }
  inline int64_t index() const {
    return _index;
  }
  
  permut_process_t *generator(tbb::flow_control &fc) {
    if (finish) {
      fc.stop();
      return nullptr;
    }
    t_vec_comb next_combs = cur_combs;
    t_vec_comb last_comb;
    _total_index += next_k_complex_permutation(next_combs, last_comb, packet_size, permi);
    permut_process_t *result;
    bool popped = prc_queue.try_pop(result);
    assert(popped);
    result->set_proc_range(cur_combs, last_comb);
    finish = next_combs == init_cmb;
    swap(cur_combs, next_combs);
    return result;
  }

  permut_process_t *indexer(permut_process_t *pp) {
    for(int i = 0; i < pp->ps_size; i++) {
      auto &x = pp->ps[i];
      x.index = _index++;
      _combination_left -= x.weight;
    }
    return pp;
  }
    
  void finalizer(permut_process_t *pp) {
    prc_queue.push(pp);
  }
    
  ~generate_comb_t() {
    while (!prc_queue.empty()) {
      permut_process_t *pp;
      bool b = prc_queue.try_pop(pp);
      assert(b);
      delete pp;
    }
  };
};

bool d2o_main_class::write_files(const string &output_base_name, bool dry_run, bool merge_confs)
{
  if(dry_run && (!merge_confs) )
    return true;
  
  t_vec_comb init_cmb;
  t_comb_descr psm;
  t_symm_set sms;
  std::tie(sms, init_cmb, psm) = create_init_perm_structs();
  q_energy_reduced qrd;
  if( calc_q_energy ) {
    qrd = reduce_q_matrix(psm);
  }  
  
  int64_t tot_comb = total_combinations();
  struct_processor str_proc(output_base_name, tot_comb);

  if(!dry_run && !tar_enabled() )
  {
    string del_command = "rm -f " + output_base_name + "*.cif";
    int rc = system(del_command.c_str());
    if( (verbose_level >= 2) && (rc == 0) )
      cout << "Output files was deleted successfully" << endl;
  }
  
  int syms_num = max<int>(occup_groups[0].symms_sets.size(), 1);
  ss_p.set_containers_prop(tot_comb, syms_num);

  const int packet_size = 10000;
  const int num_threads = 64;

  generate_comb_t gc(num_threads, tot_comb, packet_size, sms, psm.prm_indexes, init_cmb);

  auto gen_f = tbb::make_filter<void, permut_process_t *>(tbb::filter::serial_in_order,
                                                          [&gc](tbb::flow_control &fc) -> permut_process_t * {
                                                            return gc.generator(fc);
                                                          });

  auto parr_proc_f =
      tbb::make_filter<permut_process_t *, permut_process_t *>(
          tbb::filter::parallel,
          [dry_run, &qrd, merge_confs, this](permut_process_t *p) -> permut_process_t * {
    if( merge_confs )
              p->process_merge();
    else
              p->process_no_merge();
            if (!dry_run && calc_q_energy) {
              for (int i = 0; i < p->ps_size; i++)
                p->ps[i].energy = qrd.calc_energy(p->ps[i].cmb);
    }
            return p;
          });
    
  auto ctime = std::chrono::system_clock::now();
      
  auto serial_after_proc_f = tbb::make_filter<permut_process_t *, void>(
      tbb::filter::serial_in_order,
      [this, &gc, dry_run, &str_proc, &psm, &ctime, tot_comb](permut_process_t *p) -> void {
        gc.indexer(p);
        if( !dry_run ) {
          for (int i = 0; i < p->ps_size; i++) {
            const auto &x = p->ps[i];
            if ( create_q_file ) {
              string fname_str = str_proc.file_name(x);
              f_q_calc << struct_processor::get_energy_line(fname_str, x) + "\n";
        }
            if (ss_p.sampling_active()) {
              ss_p.add_structure(x);
            } else {
              write_struct(str_proc, x, psm);
      }
    }
    }  
        if( verbose_level >= 2 ) {
          auto mt = std::chrono::system_clock::now();
          double dt = std::chrono::duration<double>(mt - ctime).count();
          if( dt > 1.0 ) {
            ctime = mt;
            constexpr double mx = 5.0;
            double tm = double(gc.total_index()) + mx * double(tot_comb - gc.combination_left());
            cout << "Finished " << round(1000.0 * tm / ((mx + 1) * tot_comb)) / 10.0 << "%. "
                 << "Stored " << gc.index() << " configurations. Left "
                 << gc.combination_left() << "          \r"<< std::flush;
    }  
      }
        gc.finalizer(p);
      });

  auto pipeline_start_time = std::chrono::system_clock::now();
  tbb::parallel_pipeline(num_threads, gen_f & parr_proc_f & serial_after_proc_f);
  auto pipeline_finish_time = std::chrono::system_clock::now();
  
  if( verbose_level >= 2)
    cout <<  endl;
  
  if( verbose_level >= 1 ) {
    double dt = std::chrono::duration<double>(pipeline_finish_time - pipeline_start_time).count();
    double dt_h = std::floor(dt / 3600.0);
    dt -= dt_h * 3600;
    double dt_m = std::floor( dt / 60.0);
    dt -= dt_m * 60;
    std::stringstream ss;
    ss << "Total enumeration time: " << dt_h << ":"
       << std::setfill('0') << std::setw(2) << dt_m << ":" << dt;
    std::cout << ss.str() << std::endl;
  }


  if( gc.total_index() != tot_comb )
  {  
    cerr << "ERROR: Number of combinations (" << tot_comb
         << ")  is not equal of total index (" << gc.total_index() << ")." << std::endl;
    return false;
  }  
  
  if(merge_confs && (verbose_level >= 1) )
    cout << "Combinations after merge: " << gc.index() << endl;
  
  if(gc.combination_left() != 0)
    cerr << "ERROR: Combination left " << gc.combination_left() << " != 0 " << endl;
  
  if( ss_p.sampling_active() && (!dry_run) )
    store_sampling(output_base_name, psm, tot_comb);
  
  return gc.combination_left() == 0;
}

std::pair<int, int> d2o_main_class::correct_rms_range(const int total_sites,
                                       const double occup, 
                                                      const double x2)
{
  int min_value = total_sites;
  int max_value = 0;
  double x2_min = 1E18;
  int x2_min_index = -1;

  for(int i = 0; i <= total_sites; i++)
  {
    double x2_curr = pow( double(i)/double(total_sites) - occup, 2);
    
    if(x2_min > x2_curr)    
    {
      x2_min = x2_curr;
      x2_min_index = i;
    }  
    
    if(x2_curr < x2)
    {  
      min_value = min(min_value, i);  
      max_value = max(max_value, i);
    }          
  }
  
  min_value = min(min_value, x2_min_index);  
  max_value = max(max_value, x2_min_index);
  
  assert(min_value <= max_value);
  assert((min_value >= 0) && (max_value <= total_sites));
  return {min_value, max_value};
}

std::vector< d2o_main_class::rangi > d2o_main_class::get_rangi_array(const double x2)
{
  vector< rangi > rc;

  for(int i = 0; i < occup_groups.size(); i++)
  {
    c_occup_group &curr_group = occup_groups[i];
    
    if( (curr_group.items.size() == 1) && 
        (!(*manual_properties)[curr_group.items[0].label].population.is_initialized()) &&
        (abs(1 - curr_group.get_total_occup_input()) < 1E-4) )
    {
      rangi rd;
      rd.group_index = i;
      rd.atom_index = 0;
      rd.min_value = curr_group.number_of_sites();
      rd.max_value = rd.min_value;
      rd.curr_value = rd.min_value;
      rc.push_back(rd);
    }
    else
    {  
      for(int j = 0; j < occup_groups[i].items.size(); j++)
      {
        rangi rd;
        rd.group_index = i;
        rd.atom_index = j;
        if (!(*manual_properties)[curr_group.items[j].label].population.is_initialized()) {
          std::tie(rd.min_value, rd.max_value) =
          correct_rms_range(curr_group.number_of_sites(), 
                                occup_groups[i].items[j].occup_target, x2);
          rd.curr_value = rd.min_value;
        }
        else
        {
          int value = (*manual_properties)[curr_group.items[j].label].population.get();
          rd.min_value  = value;
          rd.max_value  = value;
          rd.curr_value = value;
        }
        rc.push_back(rd);        
      }
    }  
  }  
  
  return rc;
}

int64_t d2o_main_class::total_combinations()
{
  bool overflow = false;
  std::vector<int64_t> tm;
  for(vector<c_occup_group>::const_iterator it = occup_groups.begin(); 
                                            it != occup_groups.end(); ++it)
  {
    int64_t tc = it->get_number_of_combinations();
    overflow = overflow || (tc < 0);
    tm.push_back(tc);
  }
  
  if(overflow)
    return -1;
  else
  {
    int64_t result = 1;
    bool good = safe_multiplication(tm, result);
    return good ? result : -1;
  }
}

double d2o_main_class::ss_charge_by_occup_groups()
{
  double result = 0;
  for(vector<c_occup_group>::const_iterator it = occup_groups.begin(); 
                                            it != occup_groups.end(); ++it)
  {
    for(int j = 0; j < it->items.size(); j++)
      result += it->items[j].charge * it->items[j].num_of_atoms_sc;
  }
  
  return result;  
}

bool d2o_main_class::get_atoms_population()
{
  const double x2_min_inf = 1E18;
  vector< rangi > rc_min;
  double x2_min = x2_min_inf;
  double x2_range = 1E-5;

  while(x2_min > x2_range)
  {  
    vector< rangi > rc = get_rangi_array(x2_range);

    int64_t num_of_cmb_to_check = 1;
    
    for(int i = 0; i < rc.size(); i++)
      num_of_cmb_to_check *= rc[i].max_value - rc[i].min_value + 1;
    
    if( (num_of_cmb_to_check > 2E6 ) && (verbose_level >= 1) )
    {
      cout << "WARN: The number of combinations to check is high " 
           << num_of_cmb_to_check << endl;
      cout << "WARN: Structure is too big and/or problems with charges values." 
           << num_of_cmb_to_check << endl;
    }  
    
    if(verbose_level >= 4)
    {  
      cout << "RMS to calculate: " << x2_range << " on " << x2_min << endl;
      cout << "Combinations to check: " << num_of_cmb_to_check << endl;
    }  
    
    bool last_conf = false;
    while(!last_conf)
    {

      //check charge balance and occupancy
      vector<int> ocp_t;
      ocp_t.resize(occup_groups.size(), 0);
      double charge = 0.0;
      for(int i = 0; i < rc.size(); i++)
      {  
        const c_occup_item &cp = occup_groups[rc[i].group_index].items[rc[i].atom_index];
        charge += cp.charge * rc[i].curr_value;
        ocp_t[rc[i].group_index] += rc[i].curr_value;
        if(verbose_level >= 5)
          cout << "rc " << cp.label << ": " << rc[i].min_value  << " < " << 
                                               rc[i].curr_value << " < " << 
                                               rc[i].max_value  << endl; 
      }

      //check occupancy 
      bool overoccup = false;
      for(int i = 0; i < occup_groups.size(); i++)
      {
        overoccup = overoccup || (occup_groups[i].number_of_sites() < ocp_t[i]);
        if( overoccup ) 
        {
          if(verbose_level >= 5)
            cout << "Over occup" << endl;
          break;
        }  
      }
      
      bool underoccup = false;

      for(int i = 0; i < occup_groups.size(); i++)
      {
        bool man_occup_group = false;
        
        for(int j = 0; j < occup_groups[i].items.size(); j++)
        {
           man_occup_group = (*manual_properties)[occup_groups[i].items[j].label].population.is_initialized();
           if( man_occup_group )
             break;
        }
        if( man_occup_group )
          underoccup = false;
        else
        {  
          underoccup = (occup_groups[i].number_of_sites() > ocp_t[i]) &&
                       (abs(1.0 - occup_groups[i].get_total_occup_input()) < occup_tol());
        }          
        
        if((verbose_level >= 5) && underoccup)
          cout << "Under occup" << endl;

        if( underoccup ) break;
      }
      
      if(verbose_level >= 5)
      {
        for(int i = 0; i < occup_groups.size(); i++)
          cout << "Curr occup is " << ocp_t[i] << ". Total is " << occup_groups[i].number_of_sites() << endl;
        cout << "charge: " << charge << endl;
      }  
      
      if( (!charge_balancing || (abs(charge) < charge_tol())) && (!overoccup) && (!underoccup))
      {  
        double rms_curr = 0;
        //calculate RMS
        for(int i = 0; i < rc.size(); i++)
        {  
          double group_sites = occup_groups[rc[i].group_index].number_of_sites();
          const c_occup_item & cp = occup_groups[rc[i].group_index].items[rc[i].atom_index];
          double rms_item = pow( double(rc[i].curr_value)/group_sites - cp.occup_target ,2);

          rms_curr += rms_item;
        }

        if(rms_curr < x2_min)
        {
          rc_min = rc;
          x2_min = rms_curr;
        }
      }
      //Next configuration increment
      last_conf = true;
      for(int i = rc.size() - 1; i >= 0; i--)
      {
        rc[i].curr_value++;
        if(rc[i].curr_value > rc[i].max_value)
          rc[i].curr_value = rc[i].min_value;
        else
        {
          last_conf = false;
          break;
        }  
      }  
    }
    x2_range *= 2;
  }
  
  for(int i = 0; i < rc_min.size(); i++)
    occup_groups[rc_min[i].group_index].items[rc_min[i].atom_index].num_of_atoms_sc = rc_min[i].curr_value;
  
  return rc_min.size() > 0;
}

bool d2o_main_class::process_charges(charge_balance cb)
{
  scs.clear();
  for(const auto &a : orig_cst.atoms) {
    if( scs.count(a.label) > 0 )
  {
      if(!std::isnan(scs[a.label].input_charge) || !std::isnan(a.charge))
        assert(scs[a.label].input_charge  == a.charge);

      scs[a.label].occup += a.occupancy;
      scs[a.label].cif_mult++;
    } else {
      scs[a.label].input_charge  = a.charge;
      scs[a.label].occup = a.occupancy;
      scs[a.label].cif_mult = 1;
    }  
  }
  
  for(std::map<std::string, site_charges>::iterator it = scs.begin(); it != scs.end(); ++it)
  {
    (*it).second.curr_charge = 0;
    if(! std::isnan((*it).second.input_charge) )
      (*it).second.curr_charge = (*it).second.input_charge;

    if( (*manual_properties)[(*it).first].charge.is_initialized())
      (*it).second.curr_charge = (*manual_properties)[(*it).first].charge.get();
  }

  double total_input_charge = 0;
  double total_used_charge = 0;
    
  for(std::map<std::string, site_charges>::iterator it = scs.begin(); it != scs.end(); ++it)
  {
    total_input_charge += (*it).second.input_charge * (*it).second.occup;
    total_used_charge  += (*it).second.curr_charge * (*it).second.occup;
  }

  if( verbose_level >= 1 )
    cout << "Current charge balance option is \"" << cb_names::get_name(cb) << "\"" << endl;
  
  if(( verbose_level >= 0) && (abs(total_used_charge) > charge_tol()) )
    cout << "WARN: Total charge of the system is not zero" << endl;
  
  switch(cb)
  {
    case charge_balance::cb_no:
      charge_balancing = false;
    break;
    case charge_balance::cb_yes:
      charge_balancing = true;
    break;
    case charge_balance::cb_try:
      charge_balancing = abs(total_used_charge) < charge_tol();
    break;
    default:
      assert(false);
    break;  
  }
  
  if(verbose_level >= 1)
  {  
    cout << "Total charge oxidation state (cif):  " << total_input_charge << endl;
    cout << "Total charge cell:   " << total_used_charge << endl; 
    cout << "Charge balancing:   " << (charge_balancing ? "yes" : "no") << endl;

    cout << "----------------------------------------------------------------" << endl;
    cout << "| Atom Label\t| \tcharge  \t| mult\t| occup x mult" << endl;
    cout << "| \t\t| Ox. state\t| Used\t| (cif)\t|\t\t " << endl;
    cout << "----------------------------------------------------------------" << endl;

    for(std::map<std::string, site_charges>::iterator it = scs.begin(); it != scs.end(); ++it)
    {

      cout << "|  " << (*it).first                << "\t\t|  "
                    << (*it).second.input_charge  << "\t\t|  "
                    << (*it).second.curr_charge   << "\t|  " 
                    << (*it).second.cif_mult      << "\t|  "              
                    << (*it).second.occup  << endl;
    }
    cout << "----------------------------------------------------------------" << endl;
    cout << endl;
  }
  
  return true;
}

bool d2o_main_class::fix_groups()
{
  bool result = true;
  
  for(vector< c_occup_group >::iterator itg  = occup_groups.begin();
                                        itg != occup_groups.end(); ++itg)
  {
    assert(!itg->items.empty());
    bool fixed_status = (*manual_properties)[itg->items[0].label].fixed.get_value_or(false);
    bool wrong_status = false;
    for(vector< c_occup_item >::iterator iti  = itg->items.begin();
                                         iti != itg->items.end(); ++iti)
    {
      wrong_status = fixed_status != (*manual_properties)[iti->label].fixed.get_value_or(false);
      if(wrong_status)
        break;
    }    
    
    if(wrong_status)
    {
      cerr << "Fixation error at group " << distance(occup_groups.begin(), itg) << endl;
      result = false;
      break;
    }
    else
      itg->set_fixation(fixed_status);
  }  
  
  return result;
}

class ob_comb_atoms : public points_clusters 
{
protected:
  const cryst_structure_t & cs_str;
  cryst_tools::min_dist min_dist_obm;
public:
  virtual int get_points_size() const override;
  virtual double get_distance(int i, int j) const override;
public:
  ob_comb_atoms(const cryst_structure_t & cryst_str, double tolerance) : cs_str(cryst_str) {
    min_dist_obm.set_cell(cs_str.unit_cell.cell());
    tol_list = tolerance;
  }
  Eigen::Vector3d average_vector(const cmb_group &cbg);
  void create_groups(groups_vc &vc);
};

int ob_comb_atoms::get_points_size() const
{
  return cs_str.atoms.size();
}

double ob_comb_atoms::get_distance(int i, int j) const {
  Eigen::Vector3d dist = cs_str.unit_cell.cell() * (cs_str.atoms[i].fract_pos
      - cs_str.atoms[j].fract_pos);
  dist = min_dist_obm(dist);
  return dist.norm();
}

void ob_comb_atoms::create_groups(groups_vc &vc)
{
  create_groups_internal(vc, tol_list, 2);
  assign_max_dist(vc);
}

Eigen::Vector3d ob_comb_atoms::average_vector(const cmb_group &cbg)
{
  vector<Eigen::Vector3d> vc;
    
  for(set<int>::const_iterator it  = cbg.indexes.begin();
                               it != cbg.indexes.end(); ++it)
    vc.push_back(cs_str.unit_cell.cell() * cs_str.atoms[*it].fract_pos);

  return min_dist_obm.average_vector(vc);
}


bool d2o_main_class::create_occup_groups()
{
  map<string, double> tot_occ;
  for(int i = 0; i < supercell_cst.atoms.size(); i++) {
    const auto &a = supercell_cst.atoms[i];
    if (tot_occ.count(a.label) == 0)
      tot_occ[a.label] = a.occupancy;
    else
      tot_occ[a.label] += a.occupancy;
        }
       
  groups_vc gvc;
  
  ob_comb_atoms obc(supercell_cst, r_tolerance);
  obc.create_groups(gvc);
  min_dist_between_groups = obc.min_dist_between_groups(gvc);
  
  assert(min_dist_between_groups > r_tolerance);
  
  //Check that all groups within tolerance
  bool unique_groups = true;
  for(int i = 0; i < gvc.size(); i++)  
  {
    if(gvc[i].max_dist > r_tolerance)
    {
      unique_groups = false;
      cerr << "ERROR: Group has no unique connection. " << endl;
    }  
  }
  
  if(!unique_groups)
    return false;
  
  typedef map< set<string>, c_occup_group> cc;
  cc coc;
  
  for(int i = 0; i < gvc.size(); i++) {
    Eigen::Vector3d avg_dist = obc.average_vector(gvc[i]);
    set<string> sc;
    c_occup_group ocg_temp;
    ocg_temp.max_dis_within_group = 0;
    for(set<int>::const_iterator it  = gvc[i].indexes.begin();
                                 it != gvc[i].indexes.end(); ++it) {
      const auto &oa = supercell_cst.atoms[*it];
      if( sc.insert(oa.label).second ) //element inserted
        ocg_temp.add_item(oa, scs[oa.label].curr_charge);
    }
    
    if( coc.count(sc) == 0 )
      coc[sc] = ocg_temp;

    coc[sc].positions.push_back(avg_dist);
    coc[sc].max_dis_within_group = max(coc[sc].max_dis_within_group, gvc[i].max_dist);
  }
  
  //Check that groups are single
  vector<string> vc;
  for(cc::const_iterator it  = coc.begin();
                         it != coc.end(); ++it)
    std::copy(it->first.begin(), it->first.end(), std::back_inserter(vc));
  
  array_common::delete_singles(vc);
  
  if(vc.size() > 0) {
    for(int i = 0; i < vc.size(); i++)
      cerr << "ERROR: Label " << vc[i] << " belong to 2 or more groups." << endl;
    return false;
  }

  occup_groups.clear();
  for(cc::const_iterator it  = coc.begin();
                         it != coc.end(); ++it)
    occup_groups.push_back(it->second);

  //Correct occup_target 
  for(int i = 0; i < occup_groups.size(); i++) {
    for(int j = 0; j < occup_groups[i].items.size(); j++)
      occup_groups[i].items[j].occup_target = 
          tot_occ[occup_groups[i].items[j].label]/double(occup_groups[i].number_of_sites());
  }

  return true;
}
 
bool d2o_main_class::check_properties_consistency() const {
  bool same_properties = true;
  for(int i = 0; i < supercell_cst.atoms.size(); i++)
  {
    const auto &atom_i = supercell_cst.atoms[i];
    for(int j = i + 1; j < supercell_cst.atoms.size(); j++)
    {
      const auto &atom_j = supercell_cst.atoms[j];
      if( atom_i.label == atom_j.label ) {
        if( atom_i.el_num != atom_j.el_num ) {
          same_properties = false;
          cerr << "ERROR: Label " << atom_i.label
               << " has 2 type of atoms "
               << PT_GETSYMBOL(atom_i.el_num) << " != "
               << PT_GETSYMBOL(atom_i.el_num) << endl;
        }
        if( std::abs(atom_i.occupancy - atom_j.occupancy) > occup_tol() ) {
          same_properties = false;
          cerr << "ERROR: Label " << atom_i.label
               << " has 2 different occupations "
               << atom_i.occupancy << " != "
               << atom_j.occupancy << endl;
        }
      }
    }
  }
  return same_properties;
}

bool d2o_main_class::show_groups_information()
{
  cout << "Chemical formula of the supercell: " << get_formula_by_groups() << endl;
  cout << "Total charge of supercell: " << ss_charge_by_occup_groups() << endl;
  cout << endl ;
  cout << "----------------------------------------------------" << endl ;
  cout << " Identification of groups of crystallographic sites " << endl ;
  cout << "----------------------------------------------------" << endl ;
  
  for(int i = 0; i < occup_groups.size(); i++)
  {
    cout << endl;
    cout << " Group " << i + 1 << " (" 
         << occup_groups[i].number_of_sites() 
         << " atomic positions in supercell):" << endl;
    for(int j = 0; j < occup_groups[i].items.size(); j++)
    {
      cout << "  * Site #" << j + 1 << ": "
           << occup_groups[i].items[j].label << " (occ. " 
           << occup_groups[i].items[j].occup_target << ")";
      
      if( occup_groups[i].is_fixed() )
        cout << " -> FIXED with occupancy " << 
                boost::format("%.3f") %
                ( double(occup_groups[i].items[j].num_of_atoms_sc) /
                  double(occup_groups[i].number_of_sites() ) ) << ".";
      else
        cout << " -> distributed over " << occup_groups[i].items[j].num_of_atoms_sc << 
		" positions out of " << occup_groups[i].number_of_sites() << " (actual occ.: "
		<< boost::format("%.3f") %
		( double(occup_groups[i].items[j].num_of_atoms_sc) /
                  double(occup_groups[i].number_of_sites() ) ) << ").";
      
      cout << endl;
    }

    if( occup_groups[i].max_dis_within_group > 1E-3)
    {
      cout << "  Crystallographic sites with different positions found for this group." << endl;
      cout << "  Maximum distance within the group: " << occup_groups[i].max_dis_within_group << " A." << endl;
    }
    // else
    //   cout << "  All atoms occupied the same site." << endl;

    int64_t g_occ = occup_groups[i].get_number_of_combinations();
    
    if( g_occ > 1)
      cout << "  Number of combinations for the group is " << g_occ  << endl;
    else if( g_occ < 0 )
      cout << "  WARNING: Number of combinations for the group is too high to work with." << endl;
    //  else
    //    cout << "  The atom position within the group are set unambiguously" << endl;
    
    if( (occup_groups[i].get_total_num_occup_sites() < occup_groups[i].number_of_sites() ) &&
        (abs(1.0 - occup_groups[i].get_total_occup_input()) < occup_tol()) )
    {        
      cout << "  WARNING: Vacancy introduced in a crystallographic site, which was originally fully-occupied." << endl;
    }
    
    /* if(occup_groups[i].is_fixed())
       cout << "  FIXED: No combinations will be applied." << endl;
    */

  }

  cout << endl;
  cout << "Minimal distance between atoms of two distinct groups: " << min_dist_between_groups << " A." << endl;

  int64_t t_comb = total_combinations();
  string t_comb_approx = "";
  string t_comb_str = t_comb > 0 ? boost::lexical_cast<std::string>(t_comb) : "+INF";
  if(t_comb > 1E5)
  {  
    boost::format fmt("(~%1$2.1e)");
    fmt % double(t_comb);
    t_comb_approx = fmt.str();
  }  
  
  cout << endl ;
  cout << "-------------------------------------------------" << endl ;
  cout << "The total number of combinations is " << t_comb_str
       << t_comb_approx << endl;
  cout << "-------------------------------------------------" << endl ;
  
  return true;
}

bool d2o_main_class::create_super_cell(int a, int b, int c)
{
  supercell_cst.unit_cell.set(orig_cst.unit_cell.cell() * Eigen::DiagonalMatrix<double, 3>(a, b, c));
  supercell_cst.atoms.clear();
  supercell_cst.atoms.reserve(orig_cst.atoms.size() * a * b *c);
  for (int i = 0; i < a; i++) {
    for (int j = 0; j < b; j++) {
      for (int k = 0; k < c; k++) {
        for(const auto &oa : orig_cst.atoms) {
          supercell_cst.atoms.emplace_back(oa);
          supercell_cst.atoms.back().fract_pos = Eigen::Vector3d((oa.fract_pos.x() + i) / a,
                                                                 (oa.fract_pos.y() + j) / b,
                                                                 (oa.fract_pos.z() + k) / c);
        }
      }
    }
  }
  
  if(verbose_level >= 1)
  { 
    cout << "Initial system:" << endl;
    cout << "  Chemical Formula: " << get_formula(orig_cst) << endl;
    cout << endl;
    
    cout << "Supercell system " << boost::format("(%1%x%2%x%3%)") %a %b %c << ":" << endl;
    Eigen::Vector3d sz = supercell_cst.unit_cell.lengths();
    cout << boost::format("  Size a=%1%, b=%2%, c=%3%") % sz.x() % sz.y() % sz.z() << endl;
    cout << endl;
  }  
  
  return true;
}


bool d2o_main_class::read_cryst_structure(std::string file_name)
{
  bool result = true;
  
  string msg;
  result = read_cif_file(file_name, orig_cst, msg);
  if( verbose_level >= 1 ) {
    boost::algorithm::replace_all(msg, "\n", "\n  ");
    std::cout << "CIF file info: " << std::endl;
    std::cout << "  " << msg << std::endl;
  }
  
  return result;  
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const boost::optional<T> &cav)
{
  if(cav.is_initialized())
    os << *cav;
  else
    os << "N/A";

  return os;
}

bool d2o_main_class::set_labels_to_manual()
{
  set<string> set_lbl;
  
  for(const auto &a : orig_cst.atoms) {
    set_lbl.emplace(a.label);
  }
  
  (*manual_properties).convert_properties(set_lbl);
  
  if( (verbose_level >= 1) && (manual_properties->data().size() > 0))
  {
    cout << "Manual properties" << endl;
    cout << "Label\t|fixed\t|charge\t|popul\t|" << endl;
    for(c_man_atom_prop::data_type::const_iterator it  = manual_properties->data().begin();
                                                   it != manual_properties->data().end(); ++it)
    {
      cout << it->first << "\t|"
           << it->second.fixed      << "\t|"
           << it->second.charge     << "\t|"
           << it->second.population << "\t|" << endl;
    }
    cout << endl;
  }  
   
  return true;
}

bool d2o_main_class::process(std::string input_file_name, bool dry_run,
                             const std::vector<int> &supercell_shape,
                             charge_balance cb, double tolerance_v,
                             bool merge_confs, bool calc_q_energy_v, bool create_q_file_v,
                             c_man_atom_prop &manual_properties_v,
                             const c_struct_sel &ss_p_v,        
                             std::string output_base_name,
                             std::string output_tar_name)
{
  assert(supercell_shape.size() == 3);
  
  r_tolerance = max(tolerance_v, 1.0E-6);
  manual_properties = &manual_properties_v;
  ss_p.assign_base(ss_p_v);
  calc_q_energy = calc_q_energy_v;
  create_q_file = create_q_file_v;
          
  if(!read_cryst_structure(input_file_name))
  {
    cerr << "Input file cannot be opened." << endl;
    return false;
  }
  
  if(!set_labels_to_manual())
  {
    cerr << "Manual labels cannot be set." << endl;
    return false;
  }
  
  
  if(!create_super_cell(supercell_shape[0], supercell_shape[1], supercell_shape[2]))
  {
    cerr << "Supercell cannot be created." << endl;
    return false;
  }
  
  if( !process_charges(cb) )
  {
    cerr << "Charge processing fails." << endl;
    return false;
  }
  
  //Check, that all labels has the same properties
  if( !check_properties_consistency() ) {
    cerr << "Sites properties are not consistent." << endl;
    return false;
  }
  
  if( !create_occup_groups() )
  {
    cerr << "Create occupation groups error." << endl;
    return false;
  }
  
  if( !fix_groups() )
  {
    cerr << "Error while fixing groups." << endl;
    return false;
  }  
  
  if( !get_atoms_population() )
  {
    cerr << "Create atoms population error:" << endl <<
            "  Change supercell size"        << endl <<
            "  Switch off charge balancing"  << endl <<
            "  Check manual population settings"  << endl;
    return false;
  }
  
  if(verbose_level >= 1)
   show_groups_information();
  
  int64_t tc = total_combinations();
  if(tc > static_cast<int64_t>(1E16) || tc < 0 )
  {
    cerr << "ERROR: Number of total combinations is too high to work with." << endl;
    return false;
  }

  if( tc == 0) {
    cerr << "ERROR: Number of total combinations is 0. Probably wrong pupulation values." << endl;
    return false;
  }
  
  if( merge_confs )
  {  
    if(!create_symmetry_operations_groups())
    {
      cerr << "ERROR: Symmetry operation creation failed." << endl;
      return false;
    }
  }

  if( create_q_file && !calc_q_energy)
  {
    cerr << "ERROR: Electrostatic energy file cannot be created without electrostatic energy calculation. " << endl;
    return false;
  }

  if( calc_q_energy )  
  {
    if(!charge_balancing)
    {
      cerr << "ERROR: Electrostatic energy cannot be calculated without charge balancing. " << endl;
      return false;
    }    
    if(!calculate_q_matrix())
    {
      cerr << "ERROR: Coulomb energy is not calculated." << endl;
      return false;
    }
    
    if( !dry_run && create_q_file )
    {
      struct_processor sp(output_base_name, tc);
      string fq_name = sp.get_q_file_name("");
      f_q_calc.open(fq_name.c_str(), fstream::out);
      if(!f_q_calc.is_open())
      {
        cerr << "ERROR: File \"" << fq_name << "\" is not open." << endl;
        return false;
      }
    }
  }
  
  if( (!calc_q_energy) && ( (ss_p.str_high_count() > 0) || (ss_p.str_low_count() > 0) ) )
  {
    cerr << "Energy sampling (h or l) impossible without Coulomb calculation enabled." << endl;
    return false;
  }
  
  if(!dry_run)
    create_tar_container(output_tar_name);
   
  if(!write_files(output_base_name, dry_run, merge_confs))
  {
    cerr << "Write files error." << endl;
    return false;
  }
  
  if(!dry_run)
    close_tar_container();
  
  return true;
}
