/* 
 * File:   d2o_main_class.cpp
 * Author: kirill
 * 
 * Created on July 19, 2013, 2:25 PM
 */

#include "d2o_main_class.h"

#include <openbabel/obconversion.h>
#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <algorithm>
#include <cmath>

#include <openbabel/chargemodel.h>

#include "science/combinatorics.h"

#include "others/rnd_utils.h"
#include "others/string_utils.h"
#include "containers/array_common.hpp"

#include "obabel/eigen2babel.h"
#include "cryst_tools/comb_points.h"
#include "cryst_tools/cryst_tools.h"

#include <Eigen/Dense>


using namespace OpenBabel;
using namespace std;

d2o_main_class::d2o_main_class()
{
}

c_occup_item::c_occup_item(OpenBabel::OBAtom *ob, double charge_v)
{ 
  obp = new OpenBabel::OBAtom(); 
  obp->Duplicate(ob);
  
  label = obp->GetData("original_label")->GetValue();
  occup_target = dynamic_cast<OBPairFloatingPoint *> (obp->GetData("_atom_site_occupancy"))->GetGenericValueDef(1.0);
  
  charge = charge_v;
}

double c_occup_group::get_total_occup_input() const
{
  double result = 0.0;
  
  for(std::vector<c_occup_item>::const_iterator it = items.begin(); it != items.end(); it++)
    result += (*it).occup_target;
  
  return result;
}

int c_occup_group::get_total_num_occup_sites() const
{
  int result = 0;
  
  for(std::vector<c_occup_item>::const_iterator it = items.begin(); it != items.end(); it++)
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
      result = num_combinations(nm);
    }
    else
      result = 0;
  }
  
  return result;
}

void c_occup_group::add_item(OpenBabel::OBAtom * oba, double charge)
{
  c_occup_item item(oba, charge);
  items.push_back(item);
}

bool d2o_main_class::init_atom_change_mol(OpenBabel::OBMol *cmol)
{
  cmol->Clear();
  
  OBUnitCell *uc = new OBUnitCell(*dynamic_cast<OBUnitCell *>(mol_supercell.GetData(OBGenericDataType::UnitCell)));
  assert(!cmol->HasData(OBGenericDataType::UnitCell));
  cmol->SetData(uc);
  
  return true;
}

std::string d2o_main_class::get_formula_by_groups()
{
  map<string, int> formula_map;
  
  for(int i = 0; i < occup_groups.size(); i++)
  {
    for(int j = 0; j < occup_groups[i].items.size(); j++)
    {
      const c_occup_item &ci = occup_groups[i].items[j]; 
      string atom_symbol = etab.GetSymbol(ci.obp->GetAtomicNum());
      formula_map[atom_symbol] += ci.num_of_atoms_sc;
    }  
  }
  
  string result = "";
  
  for(map<string, int>::iterator it = formula_map.begin();
                                 it != formula_map.end(); it++)
  {
    boost::format fmter("%1%%2% ");
    fmter % it->first % it->second;
    result += fmter.str();
  }  
  
  trim(result);
  
  return result;      
}

std::string d2o_main_class::get_formula(OpenBabel::OBMol &mol)
{
  map<string, double> formula_map;
  
  for(OBAtomIterator it = mol.BeginAtoms(); it != mol.EndAtoms(); ++it)
  {
    string atom_symbol = etab.GetSymbol((*it)->GetAtomicNum());
    
    double curr_occup = dynamic_cast<OBPairFloatingPoint *>
                        ((*it)->GetData("_atom_site_occupancy"))->GetGenericValueDef(1.0);

    formula_map[atom_symbol] += curr_occup;
  }
  
  string result = "";
  
  for(map<string, double>::iterator it = formula_map.begin();
                                    it != formula_map.end(); it++)
  {
    boost::format fmter("%1%%2% ");
    fmter % it->first % it->second;
    result += fmter.str();
  }  
  
  trim(result);
  
  return result;      
}

bool d2o_main_class::add_confs_to_mol(OpenBabel::OBMol *cmol, const t_vec_comb &ppc)
{
  cmol->BeginModify();
  for(int i = 0; i < occup_groups.size(); i++)
  {
    c_occup_group &curr_group = occup_groups[i];
    bool fixed_group = curr_group.is_fixed();    
    for(int j = 0; j < curr_group.positions.size(); j++)
    {
      for(int k = 0; k < curr_group.items.size(); k++)
      {
        bool the_item = false;
	if(ppc[i].size() > 0)
	  the_item = (k == ppc[i][j]);
	
	if( fixed_group || the_item )
        {  
          OBAtom atm;
          atm.Duplicate( curr_group.items[k].obp );
          
          double occup_value;
          
          if( fixed_group )
            occup_value = double(curr_group.items[k].num_of_atoms_sc) / 
                          double(curr_group.number_of_sites());
          else
            occup_value = 1.0;

          if( atm.HasData("_atom_site_occupancy") )
            dynamic_cast<OBPairFloatingPoint *>(atm.GetData("_atom_site_occupancy"))->SetValue(occup_value);
          else
          {
            OBPairFloatingPoint * obo = new OBPairFloatingPoint();
            obo->SetAttribute("_atom_site_occupancy");
            obo->SetValue(occup_value);
            atm.SetData(obo);
          }  
  
          atm.SetVector(curr_group.positions[j]);
          cmol->AddAtom(atm, true);
        }
      }
    }  
  }
  
  cmol->EndModify();

  return true;
}

bool d2o_main_class::calculate_q_matrix()
{
  using namespace Eigen;
  using namespace cryst_tools;
  OBUnitCell * uc = static_cast<OBUnitCell *>(mol_supercell.GetData(OBGenericDataType::UnitCell));
  Matrix3d cell = b2e_matrix<double>(uc->GetCellMatrix().transpose());
  
  vector<Vector3d> all_pos;
  for(int i = 0; i < occup_groups.size(); i++)
  {
    for(int j = 0; j < occup_groups[i].positions.size(); j++)
      all_pos.push_back(b2e_vector<double>(occup_groups[i].positions[j]));
  }
  
  if(verbose_level >= 2)
    cout << "Start Coulomb matrix (" << all_pos.size() << "x" << all_pos.size() << ") calculation." << endl;

  ewald_sum es;
  
  es.set_cell(cell);
  es.set_precision(all_pos.size(), 1E-7);
  
  q_energy = es.potential_matrix(all_pos);
 
  if(verbose_level >= 2)
    cout << "Coulomb matrix calculation finished." << endl;
  
  return true;
}

/*
   for(int i = 0; i < all_pos.size(); i++)
  {
    for(int j = 0; j < all_pos.size(); j++)
    {
      vector<Vector3d> vq = 
      md.get_img_dist(all_pos[i] - all_pos[j], Vector3i(45, 45, 15));
      for(int k = 0; k < vq.size(); k++)
      {  
        double d = vq[k].norm();
        if(d > symm_tol)
        {  
          q_energy(i, j) += 1.0/d;
          //if(i != j)
            //q_energy(j, i) += 1.0/d;
        }  
      }
    }
  }

 */

double d2o_main_class::calculate_q_energy(const t_vec_comb &mc)
{
  using namespace Eigen;
  
  VectorXd q_v;
  q_v.resize(q_energy.cols());
  q_v.setZero();
  
  int q_v_pos = 0;
  
  for(int i = 0; i < occup_groups.size(); i++)
  {
    c_occup_group &curr_group = occup_groups[i];
    bool fixed_group = curr_group.is_fixed();
    for(int j = 0; j < curr_group.positions.size(); j++)
    {
      double charge = 0;
      for(int k = 0; k < curr_group.items.size(); k++)
      {
        double occup_value = 0;

        bool the_item = false;
	if(mc[i].size() > 0)
	  the_item = (k == mc[i][j]);
        
        if( fixed_group || the_item )
        {  
          if( fixed_group )
            occup_value = double(curr_group.items[k].num_of_atoms_sc) / 
                          double(curr_group.number_of_sites());
          else
            occup_value = 1.0;
        }
        charge += occup_value * curr_group.items[k].charge;
      }
      q_v[q_v_pos] = charge;
      q_v_pos++;
    }  
  }

  assert(q_v_pos == q_energy.cols());
  
  //0.5 not to count twice pairs
  //11.4 - to eV
  return 14.4 * q_v.transpose() * q_energy * q_v;
}

bool d2o_main_class::create_symmetry_operations_groups()
{
  using namespace cryst_tools;
  using namespace Eigen;

  OBUnitCell * uc = static_cast<OBUnitCell *>(mol_supercell.GetData(OBGenericDataType::UnitCell));

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
      matrix3x3 m_cell = (uc->GetCellMatrix().transpose()).inverse();
      vector3 v3 = m_cell * occup_groups[i].positions[j];
      vc[i][j] = b2e_vector<double>(v3);
    }
  }
  
  vector<Affine3d> syms;
  
  Matrix3d cell;
  cell = b2e_matrix<double>(uc->GetCellMatrix().transpose());
  
  syms = get_all_symmetries(cell, vc, bc, symm_tol);
  
  if(verbose_level >= 1)
    cout << syms.size() << " symmetry operation found for supercell." << endl;
  
  bool good_set = true;
  
  for(int i = 0; i < occup_groups.size(); i++)    
  {
    occup_groups[i].symms_sets.resize(syms.size());
    for(int j = 0; j < syms.size(); j++)
    {  
      occup_groups[i].symms_sets[j] = index_symmetries(uc, syms[j], 
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

std::vector<int> d2o_main_class::index_symmetries(OpenBabel::OBUnitCell * uc, 
                                                  const Eigen::Affine3d &af, 
                                                  const std::vector<OpenBabel::vector3> &pos)
{
  using namespace Eigen;
  
  std::vector<int> result;
  Matrix3d cell = b2e_matrix<double>(uc->GetCellMatrix().transpose());
  Matrix3d r_cell = cell.inverse();
  Affine3d cart_tr, cl_d, cl_i;
  
  cl_d.setIdentity();
  cl_d.linear() = cell;
  
  cl_i.setIdentity();
  cl_i.linear() = cell.inverse();
  
  cart_tr = cl_d * af * cl_i;
  
  std::vector<Vector3d> pose;  
  for(int i = 0; i < pos.size(); i++)
    pose.push_back(b2e_vector<double>(pos[i]));
  
  result.resize(pos.size(), -1);
  
  for(int i = 0; i < pose.size(); i++)
  {
    int index_sym = -1;
    int index_count = 0;
    for(int j = 0; j < pose.size(); j++)
    {
      Vector3d ve = pose[j] - cart_tr * pose[i];
      ve = cell * cryst_tools::min_frac(r_cell * ve);
      if(ve.norm() < symm_tol)
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

bool d2o_main_class::create_comb(const symm_set &sc, std::vector<int> &cmb)
{
  assert(sc.size() == cmb.size());
  std::vector<int> cb;
  
  cb.resize(cmb.size());
  bool result = true;
  for(int i = 0; i < cb.size(); i++)
  {  
    if(sc[i] >= 0)
      cb[sc[i]] = cmb[i];
    else
      result = result && (cmb[i] < 0);
  }  
    
  cmb = cb;
  
  return result;
}

bool d2o_main_class::check_comb_unique(const t_vec_comb &mc, int &merged_comb)
{
  assert(occup_groups.size() > 0);
  assert(mc.size() == occup_groups.size());
  int syms_num = occup_groups[0].symms_sets.size();
  assert(syms_num > 0);
  
  set< t_vec_comb > st;

  bool not_first = false;
  for(int i = 0; i < syms_num; i++)
  {
    t_vec_comb nc;
    nc.resize(mc.size());
    bool good_cb = true;
    for(int j = 0; j < mc.size(); j++)
    {
      assert(syms_num == occup_groups[j].symms_sets.size());
      nc[j] = mc[j];
      
      if(!occup_groups[j].is_fixed())
        good_cb = good_cb && create_comb(occup_groups[j].symms_sets[i], nc[j]);
    }
    if(good_cb)
    {  
      st.insert(nc);
      if(nc < mc)
      {
        not_first = true;
	break;
      }
    }  
  }
  merged_comb = st.size();
  
  //assert(*st.begin() <= mc);
  
  return !not_first;
}

bool d2o_main_class::write_files(std::string output_base_name, bool dry_run, bool merge_confs)
{
  if(dry_run && (!merge_confs) )
    return true;
  
  t_vec_comb cur_combs;
  
  cur_combs.resize(occup_groups.size());
  for(int i = 0; i < occup_groups.size(); i++)
  {  
    cur_combs[i].clear();
    if(!occup_groups[i].is_fixed())
    {  
      map<int, int> mvc;

      for(int j = 0; j < occup_groups[i].items.size(); j++)
        mvc[j] = occup_groups[i].items[j].num_of_atoms_sc;
    
      mvc[-1] = occup_groups[i].number_of_sites() - 
                        occup_groups[i].get_total_num_occup_sites();
    
      cur_combs[i] = create_start_combination(mvc);
    }
  }  
  
  double tot_comb = total_combinations();
  
  if(!dry_run)
  {
    string del_command = "rm -f " + output_base_name + "*.cif";
    int rc = system(del_command.c_str());
    if( (verbose_level >= 2) && (rc == 0) )
      cout << "Output files was deleted successfully" << endl;
  }
  
  int combination_left = total_combinations();
  int total_index = 0;
  int index = 0;
  bool done;
  
  OBMol cmol;
  do
  {
    bool u_comb;
    int merged_conf_num;

    if( merge_confs )
      u_comb = check_comb_unique(cur_combs, merged_conf_num);
    else
    {  
      u_comb = true;
      merged_conf_num = 1;
    }

    if( u_comb )
    {  
      combination_left -= merged_conf_num;
      if(!dry_run)
      {
        init_atom_change_mol(&cmol);
        add_confs_to_mol(&cmol, cur_combs);

        OBConversion obc;

        string index_str = get_index_str(index, tot_comb - 1);
        string fname_str = output_base_name + "_ind" + index_str;
        if( merge_confs ) 
          fname_str += "w_" + boost::lexical_cast<string>(merged_conf_num);

        obc.SetOutFormat("cif");
        obc.WriteFile(&cmol, fname_str + ".cif");
        if(calc_q_energy)
          f_q_calc << boost::format("%1%\t%2$.3f eV\n") %
          fname_str % calculate_q_energy(cur_combs);
      }
      index++;
    }

    if( (total_index % 2000 == 0) && (total_index != 0) && (verbose_level >= 2))
    {  
      cout << "Finished " << round(double(total_index) / tot_comb * 100) << "%. " 
           << "Stored " << index << " configurations. Left " << combination_left << "          \r"<< endl;
    }  

    //Next combination
    done = true;    
    for(t_vec_comb::reverse_iterator rit  = cur_combs.rbegin(); 
                                     rit != cur_combs.rend(); ++rit)
    {
      if( std::next_permutation(rit->begin(), rit->end()) )
      {
        done = false;
        break;
      }
    }
    total_index++;
  }while(!done);
  
  if( verbose_level >= 2)
    cout <<  endl;
  

  if( total_index != total_combinations() )
  {  
    cerr << "ERROR: Number of combinations is not equal of total index." << endl;
    return false;
  }  
  
  if(merge_confs && (verbose_level >= 1) )
    cout << "Combinations after merge: " << index << endl;
  
  if(combination_left != 0)
    cerr << "ERROR: Combination left " << combination_left << " != 0 " << endl;
  
  return combination_left == 0;
}

void d2o_main_class::correct_rms_range(const int total_sites, 
                                       const double occup, 
                                       const double x2,
                                       int &min_value,
                                       int &max_value)
{
  min_value = total_sites;
  max_value = 0;
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
}

std::vector< d2o_main_class::rangi > d2o_main_class::get_rangi_array(const double x2)
{
  vector< rangi > rc;

  for(int i = 0; i < occup_groups.size(); i++)
  {
    c_occup_group &curr_group = occup_groups[i];
    
    if( (curr_group.items.size() == 1) && 
        (!(*manual_properties)[curr_group.items[0].label].population.assigned()) &&
        (abs(1 - curr_group.get_total_occup_input()) < 1E-4) )
    {
      rangi rd;
      rd.group_index = i;
      rd.atom_index = 0;
      rd.min_value = curr_group.number_of_sites();
      rd.max_value = curr_group.number_of_sites();
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
        if(!(*manual_properties)[curr_group.items[j].label].population.assigned())
        {  
          correct_rms_range(curr_group.number_of_sites(), 
                            occup_groups[i].items[j].occup_target, x2,
                            rd.min_value, rd.max_value);
          rd.curr_value = rd.min_value;
        }
        else
        {
          int value = (*manual_properties)[curr_group.items[j].label].population.value();
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
  int64_t result = 1;
  for(vector<c_occup_group>::const_iterator it = occup_groups.begin(); 
                                            it != occup_groups.end(); it++)
  {
    result *= it->get_number_of_combinations();
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
        c_occup_item cp = occup_groups[rc[i].group_index].items[rc[i].atom_index];
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
           man_occup_group = (*manual_properties)[occup_groups[i].items[j].label].population.assigned();
           if( man_occup_group )
             break;
        }
        if( man_occup_group )
          underoccup = false;
        else
        {  
          underoccup = (occup_groups[i].number_of_sites() > ocp_t[i]) &&
                       (abs(1.0 - occup_groups[i].get_total_occup_input()) < occup_tol);
        }          
        
        if((verbose_level >= 5) && underoccup)
          cout << "Under occup" << endl;

        if( underoccup ) break;
      }
      
      if(verbose_level >= 5)
      {
        for(int i = 0; i < occup_groups.size(); i++)
          cout << "Curr occup is " << ocp_t[i] << ". Total is " << occup_groups[i].number_of_sites() << endl;
      }  
      
      if(verbose_level >= 5)
        cout << "charge: " << charge << endl;

      
      if( (abs(charge) < charge_tol) && (!overoccup) && (!underoccup))
      {  
        double rms_curr = 0;
        //calculate RMS
        for(int i = 0; i < rc.size(); i++)
        {  
          double group_sites = occup_groups[rc[i].group_index].number_of_sites();
          c_occup_item cp = occup_groups[rc[i].group_index].items[rc[i].atom_index];
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
  
  return true;
}

bool d2o_main_class::process_charges(charge_balance cb, bool verbose)
{
  FOR_ATOMS_OF_MOL(a, mol_initial)
  {
    string label = a->GetData("_atom_site_label")->GetValue();
    double curr_input_charge = dynamic_cast<OBPairFloatingPoint *>(a->GetData("input_charge"))->GetGenericValueDef(NAN);
    double curr_occup = dynamic_cast<OBPairFloatingPoint *>(a->GetData("_atom_site_occupancy"))->GetGenericValueDef(1.0);
    if( scs.count(label) > 0 )
    {  
      if(!isnan(scs[label].input_charge) || !isnan(curr_input_charge))
        assert(scs[label].input_charge  == curr_input_charge);

      scs[label].occup += curr_occup;
      scs[label].cif_mult++;
    }
    else
    {
      scs[label].input_charge  = curr_input_charge;
      scs[label].occup = curr_occup; 
      scs[label].cif_mult = 1;
    }  
  }
  
  for(std::map<std::string, site_charges>::iterator it = scs.begin(); it != scs.end(); it++)
  {
    switch(cb)
    {        
      case cb_no:
        (*it).second.curr_charge = 0;
      break;  

      case cb_input:
      case cb_try:  
        if(! isnan((*it).second.input_charge) )
          (*it).second.curr_charge = (*it).second.input_charge;
      break;  
      
      default:
        assert(false);
      break;  
    }
    
    if( cb != cb_no )
    {  
      if( (*manual_properties)[(*it).first].charge.assigned()) 
        (*it).second.curr_charge = (*manual_properties)[(*it).first].charge.value();
    }   
  }

  double total_input_charge = 0;
  double total_used_charge = 0;
    
  for(std::map<std::string, site_charges>::iterator it = scs.begin(); it != scs.end(); it++)
  {
    total_input_charge += (*it).second.input_charge * (*it).second.occup;
    total_used_charge  += (*it).second.curr_charge * (*it).second.occup;
  }

  if( verbose_level >= 1 )
    cout << "Current charge balance option is \"" << cb_names::get_name(cb) << "\"" << endl;
  
  if((verbose >= 0) && (abs(total_used_charge) > charge_tol) )
    cout << "WARN: Total charge of the system is not zero" << endl;
  
  if( (cb == cb_try ) && (abs(total_used_charge) > charge_tol) )
  {
    for(std::map<std::string, site_charges>::iterator it = scs.begin(); it != scs.end(); it++)
      (*it).second.curr_charge = 0;
    if(verbose_level >= 1)
      cout << "Charge balancing is switched off." << endl;
    total_used_charge = 0;
  }

  
  if(verbose_level >= 1)
  {  
    cout << "Total charge oxydation state (cif):  " << total_input_charge << endl;
    cout << "Total charge used:   " << total_used_charge << endl << endl;

    cout << "| Atom Label\t| \tcharge  \t| mult\t| occup x mult" << endl;
    cout << "| \t\t| Ox. state\t| Used\t| (cif)\t|\t\t " << endl;
    for(std::map<std::string, site_charges>::iterator it = scs.begin(); it != scs.end(); it++)
    {

      cout << "|  " << (*it).first                << "\t\t|  "
                    << (*it).second.input_charge  << "\t\t|  "
                    << (*it).second.curr_charge   << "\t|  " 
                    << (*it).second.cif_mult      << "\t|  "              
                    << (*it).second.occup  << endl;
    }
    cout << endl;
  }
  
  return true;
}

bool d2o_main_class::fix_groups()
{
  bool result = true;
  
  for(vector< c_occup_group >::iterator itg  = occup_groups.begin();
                                        itg != occup_groups.end(); itg++)
  {
    assert(itg->items.size() > 0);
    bool fixed_status = (*manual_properties)[itg->items[0].label].fixed.value_def(false);
    bool wrong_status = false;
    for(vector< c_occup_item >::iterator iti  = itg->items.begin();
                                         iti != itg->items.end(); iti++)
    {
      wrong_status = fixed_status != (*manual_properties)[iti->label].fixed.value_def(false);
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
  OBMol * obm;
  ob_min_dist min_dist_obm;  
public:
  virtual int get_points_size() const;
  virtual double get_distance(int i, int j) const;
public:
  ob_comb_atoms(OBMol * mol_v, double tolerance);
  vector3 average_vector(const cmb_group &cbg);
  bool total_shift(const cmb_group &cbg, int poins_num);
  void create_groups(groups_vc &vc);
};

ob_comb_atoms::ob_comb_atoms(OBMol * mol_v, double tolerance)
{
  obm = mol_v;
  OBUnitCell * uc = (OBUnitCell *) obm->GetData(OBGenericDataType::UnitCell);
  tol_list = tolerance;
  min_dist_obm.set_cell(uc->GetCellMatrix().transpose());
}

int ob_comb_atoms::get_points_size() const
{
  return obm->NumAtoms();
}

double ob_comb_atoms::get_distance(int i, int j) const
{
  vector3 dist;
  
  dist = obm->GetAtom(i + 1)->GetVector() - 
         obm->GetAtom(j + 1)->GetVector();
  
  dist = min_dist_obm(dist);
  
  return dist.length();
}

void ob_comb_atoms::create_groups(groups_vc &vc)
{
  create_groups_internal(vc, tol_list, 2);
  assign_max_dist(vc);
}

vector3 ob_comb_atoms::average_vector(const cmb_group &cbg)
{
  vector<vector3> vc;
    
  for(set<int>::const_iterator it  = cbg.indexes.begin();
                               it != cbg.indexes.end(); ++it)
    vc.push_back( obm->GetAtom(*it + 1)->GetVector() );

  return min_dist_obm.average_vector(vc);
}


bool d2o_main_class::create_occup_groups()
{
  //Check, that all labels has the same properties
  bool same_properties = true;
  for(int i = 0; i < mol_supercell.NumAtoms(); i++)
  {
    OBAtom * atom_i = mol_supercell.GetAtom(i + 1);
    string label_i = atom_i->GetData("original_label")->GetValue();
    double occup_i = dynamic_cast<OBPairFloatingPoint *> (atom_i->GetData("_atom_site_occupancy"))->GetGenericValueDef(1.0);
    for(int j = i + 1; j < mol_supercell.NumAtoms(); j++)
    {
      OBAtom * atom_j = mol_supercell.GetAtom(j + 1);
      string label_j = atom_j->GetData("original_label")->GetValue();
      if(label_i == label_j)
      {
        if( atom_i->GetAtomicNum() != atom_j->GetAtomicNum() )
        {
          same_properties = false;
          cerr << "ERROR: Label " << label_i
               << " has 2 type of atoms " 
               << etab.GetSymbol(atom_i->GetAtomicNum()) << " != "   
               << etab.GetSymbol(atom_j->GetAtomicNum()) << endl;
        }
        
        double occup_j = dynamic_cast<OBPairFloatingPoint *> (atom_j->GetData("_atom_site_occupancy"))->GetGenericValueDef(1.0);
        if( abs(occup_i - occup_j) > occup_tol )
        {
          same_properties = false;
          cerr << "ERROR: Label " << label_i
               << " has 2 different occupations "
               << occup_i << " != "   
               << occup_j << endl;
        }
      }
    }
  }
  
  if(!same_properties)
    return false;
  
  groups_vc gvc;
  
  ob_comb_atoms obc(&mol_supercell, r_tolerance);
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
  
  for(int i = 0; i < gvc.size(); i++)  
  {
    vector3 avg_dist = obc.average_vector(gvc[i]);
    set<string> sc;
    c_occup_group ocg_temp;
    ocg_temp.max_dis_within_group = 0;
    for(set<int>::const_iterator it  = gvc[i].indexes.begin();
                                 it != gvc[i].indexes.end(); ++it)
    {
      OBAtom * a = mol_supercell.GetAtom(*it + 1);
      string label = a->GetData("original_label")->GetValue();
      if( sc.insert(label).second ) //element inserted
        ocg_temp.add_item(a, scs[label].curr_charge);
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
  
  if(vc.size() > 0)
  {
    for(int i = 0; i < vc.size(); i++)
      cerr << "ERROR: Label " << vc[i] << " belong to 2 or more groups." << endl;
    return false;
  }

  occup_groups.clear();
  for(cc::const_iterator it  = coc.begin();
                         it != coc.end(); ++it)
    occup_groups.push_back(it->second);
  
  return true;
}
 
/*
bool d2o_main_class::create_occup_groups()
{
  OBUnitCell *unitcell = (OBUnitCell *)mol_supercell.GetData(OBGenericDataType::UnitCell);
  
  set<OBAtom *> to_delete;
  
  //delete elements of the same label, which close to each other.  
  for(int i = 0; i < mol_supercell.NumAtoms(); i++)
  {
    OBAtom * atom_i = mol_supercell.GetAtom(i + 1);
    for(int j = i + 1; j < mol_supercell.NumAtoms(); j++)
    {
      OBAtom * atom_j = mol_supercell.GetAtom(j + 1);
      string label_i = atom_i->GetData("original_label")->GetValue();
      vector3 dist = atom_i->GetVector() - atom_j->GetVector();
      dist = get_minimal_distance(dist, unitcell);
      if(dist.length() < r_tolerance)
      {
        string label_j = atom_j->GetData("original_label")->GetValue();
        if((label_i == label_j))
          to_delete.insert(atom_j);
      }
    }
  }
  
  for(set<OBAtom *>::iterator it = to_delete.begin(); it != to_delete.end(); it++)
    mol_supercell.DeleteAtom(*it, false);
  
  
  vector< vector<c_occup_item> > element_data;
  vector<int> first_element;
  
  element_data.clear();
  element_data.resize(mol_supercell.NumAtoms());
  
  first_element.clear();
  first_element.resize(mol_supercell.NumAtoms(), 10000000);
  
  //Set closest array
  for(int i = 0; i < mol_supercell.NumAtoms(); i++)
  {
    OBAtom * atom_i = mol_supercell.GetAtom(i + 1);
    for(int j = 0; j < mol_supercell.NumAtoms(); j++)
    {
      OBAtom * atom_j = mol_supercell.GetAtom(j + 1);      
      vector3 dist = atom_i->GetVector() - atom_j->GetVector();
      dist = get_minimal_distance(dist, unitcell);
      if(dist.length() < r_tolerance)
      {  
        c_occup_item tp;
        tp.label = atom_j->GetData("original_label")->GetValue();
        tp.occup_target = dynamic_cast<OBPairFloatingPoint *> (atom_j->GetData("_atom_site_occupancy"))->GetGenericValue();
        tp.charge = scs[tp.label].curr_charge;
        tp.obp->Duplicate(atom_j);
        element_data[i].push_back(tp);
      }  
    }  
  }

  //check uniform occupation. 
  for(int i = 0; i < element_data.size(); i++)
  {
    assert(element_data[i].size() > 0);
    for(int j = 0; j < element_data.size(); j++)
    {
      bool find_all = true;
      bool find_one = false;
      
      for(int k = 0; k < element_data[i].size(); k++)
      {
        int find_count = 0;        
        for(int l = 0; l < element_data[j].size(); l++)
        {
          if(element_data[i][k].label == element_data[j][l].label)
          {  
            assert(element_data[i][k].occup_target == element_data[j][l].occup_target);
            find_count++;
          }  
        }
        
        if(find_count > 1)
        {
          cerr << "Multiple definition of label " + element_data[i][k].label << endl;
          
          cout << j << endl;          
          for(int l = 0; l < element_data[j].size(); l++)
            cout << element_data[j][l].label << endl;
          
          return false;
        }
        
        find_all = find_all && (find_count == 1);
        find_one = find_one || (find_count == 1);
      }
      
      assert((find_all != true) || (find_one != false));
      
      if(find_all != find_one)
      {
        cerr << "Sites are wrong. Check input file " << endl;
        return false;
      }
      
      if(find_all && find_one)
      {
        int index;
        
        index = min(i, j);
        index = min(index, first_element[i]);
        index = min(index, first_element[j]);
        
        first_element[i] = index;
        first_element[j] = index;
      }
    }
  }
  
  occup_groups.clear();
  vector<int> elem_group_num;
  
  elem_group_num.resize(element_data.size(), -1);
  
  //Create set of indexes and delete close sites;
  for(int i = 0; i < element_data.size(); i++)
  {  
    //sort(element_data[i].begin(), element_data[i].end());
    assert(first_element[i] <= i);
    if(first_element[i] == i)
    {  
      c_occup_group cop;
      cop.number_of_sites = 0;
      cop.max_dis_within_group = 0.0;
      cop.items = element_data[i];
      occup_groups.push_back(cop);
      elem_group_num[i] = occup_groups.size() - 1;
    }
    else
      elem_group_num[i] = elem_group_num[first_element[i]];
  }
  
  for(int i = 0; i < elem_group_num.size(); i++)
    assert((elem_group_num[i] >= 0) && (elem_group_num[i] < occup_groups.size()));

  //Assign group number to atoms 
  for(int i = 0; i < mol_supercell.NumAtoms(); i++)
  {
    OBPairInteger *obo = new OBPairInteger;

    obo->SetAttribute("group_number");
    obo->SetValue(elem_group_num[i]);
          
    mol_supercell.GetAtom(i + 1)->SetData(obo);
  }

  min_dist_between_groups = 1000.0;  
  //Delete close atoms and set maximum distance within group
  for(int i = 0; i < mol_supercell.NumAtoms(); i++)
  {
    OBAtom * atom_i = mol_supercell.GetAtom(i + 1);
    vector<vector3> pos_near;
    pos_near.push_back(atom_i->GetVector());
    for(int j = i + 1; j < mol_supercell.NumAtoms(); j++)
    {
      OBAtom * atom_j = mol_supercell.GetAtom(j + 1);      
      vector3 dist = atom_i->GetVector() - atom_j->GetVector();
      dist = get_minimal_distance(dist, unitcell);
      if(dist.length() < r_tolerance)
      {  
        int group_index = dynamic_cast<OBPairInteger *> (atom_i->GetData("group_number"))->GetGenericValue();
        occup_groups[group_index].max_dis_within_group = max(occup_groups[group_index].max_dis_within_group, dist.length());
        pos_near.push_back(atom_j->GetVector());
        mol_supercell.DeleteAtom(atom_j);
      }
      else
        min_dist_between_groups = min(min_dist_between_groups, dist.length());
    }  
    atom_i->SetVector(center_mass(pos_near, unitcell, r_tolerance));
  }
  
  //calculate statistic for atoms of the same type
  for(int i = 0; i < mol_supercell.NumAtoms(); i++)
  {
    OBAtom * atom_i = mol_supercell.GetAtom(i + 1);
    
    vector<vector3> pos_near;
    pos_near.push_back(atom_i->GetVector());
    
    for(set<OBAtom *>::iterator it = to_delete.begin(); it != to_delete.end(); it++)
    {  
      OBAtom * atom_j = *it;      
      vector3 dist = atom_i->GetVector() - atom_j->GetVector();
      dist = get_minimal_distance(dist, unitcell);
      if(dist.length() < r_tolerance)
      {  
        int group_index = dynamic_cast<OBPairInteger *> (atom_i->GetData("group_number"))->GetGenericValue();
        occup_groups[group_index].max_dis_within_group = max(occup_groups[group_index].max_dis_within_group, dist.length());
        pos_near.push_back(atom_j->GetVector());
      }
      else
        min_dist_between_groups = min(min_dist_between_groups, dist.length());
    }  
    atom_i->SetVector(center_mass(pos_near, unitcell, r_tolerance));
    
  }
  for(set<OBAtom *>::iterator it = to_delete.begin(); it != to_delete.end(); it++)
    delete *it;  

  //Calculate available positions for each group
  for(OBAtomIterator it = mol_supercell.BeginAtoms(); it != mol_supercell.EndAtoms(); it++)
  {
    int group_index = dynamic_cast<OBPairInteger *> ((*it)->GetData("group_number"))->GetGenericValue();
    assert( (group_index >= 0) && (group_index < occup_groups.size()) );
    occup_groups[group_index].number_of_sites++;
  }
  
  return true;
}
*/

bool d2o_main_class::show_groups_information()
{
  cout << "Chemical formula of the supercell: " << get_formula_by_groups() << endl;
  cout << "Minimal distance between groups: " << min_dist_between_groups << endl;
  
  for(int i = 0; i < occup_groups.size(); i++)
  {
    cout << " Group " << i << " has " 
         << occup_groups[i].number_of_sites() 
         << " sites:" << endl;
    for(int j = 0; j < occup_groups[i].items.size(); j++)
    {
      cout << "   Atom #" << j + 1 << " - "
           << occup_groups[i].items[j].label << "(occup: " 
           << occup_groups[i].items[j].occup_target << ")";
      
      if( occup_groups[i].is_fixed() )
        cout << " stored with occupancy " << 
                boost::format("%.3f") %
                ( double(occup_groups[i].items[j].num_of_atoms_sc) /
                  double(occup_groups[i].number_of_sites() ) );
      else
        cout << " has now " << occup_groups[i].items[j].num_of_atoms_sc << " sites.";
      
      cout << endl;
    }

    if( occup_groups[i].max_dis_within_group > 1E-3)
      cout << "  Maximum distance within the group is " << occup_groups[i].max_dis_within_group << endl;
    else
      cout << "  All atoms occupied the same site." << endl;
    
    if(occup_groups[i].get_number_of_combinations() != 1)
      cout << "  Number of combinations for the group is " << occup_groups[i].get_number_of_combinations() << endl;
    else
      cout << "  The atom position within the group are set unambiguously" << endl;
    
    if( (occup_groups[i].get_total_num_occup_sites() < occup_groups[i].number_of_sites() ) &&
        (abs(1.0 - occup_groups[i].get_total_occup_input()) < occup_tol) )
    {        
      cout << "  WARN: Vacancy introduced to fully occupied state." << endl;
    }
    
    if(occup_groups[i].is_fixed())    
      cout << "  FIXED: No combinations will be applied." << endl;
    
    cout << endl;
  }
  
  cout << "Total combinations is " << total_combinations() << endl;
  
  return true;
}

bool d2o_main_class::create_super_cell(int a, int b, int c)
{
  OBUnitCell *orig_unitcell = (OBUnitCell *)mol_initial.GetData(OBGenericDataType::UnitCell);
  
  vector<vector3>  cellVectors = orig_unitcell->GetCellVectors();
  
  OBUnitCell *super_unitcell = new OBUnitCell();
  super_unitcell->SetData(cellVectors[0] * a, cellVectors[1] * b, cellVectors[2] * c);
  mol_supercell.SetData(super_unitcell);
  
  orig_unitcell->FillUnitCell(&mol_initial);
  
  super_unitcell->SetSpaceGroup(1);
  
  for (int i = 0; i < a; i++) 
  {
    for (int j = 0; j < b; j++)  
    {
      for (int k = 0; k < c; k++)  
      {
        vector3 disp(  cellVectors[0].x() * i
                     + cellVectors[1].x() * j
                     + cellVectors[2].x() * k,
                       cellVectors[0].y() * i
                     + cellVectors[1].y() * j
                     + cellVectors[2].y() * k,
                       cellVectors[0].z() * i
                     + cellVectors[1].z() * j
                     + cellVectors[2].z() * k );
          
        for(OBAtomIterator it = mol_initial.BeginAtoms(); it != mol_initial.EndAtoms(); ++it)
        {
          OBAtom *new_atom = mol_supercell.NewAtom();
          new_atom->Duplicate(*it);
          
          //cout << dynamic_cast<OBPairFloatingPoint *> ((*it)->GetData("_atom_site_occupancy"))->GetGenericValue() << endl;
          
          string label = (*it)->GetData("_atom_site_label")->GetValue();
          
          OBPairData *obo = new OBPairData;

          obo->SetAttribute("original_label");
          obo->SetValue(label);
          
          new_atom->SetData(obo);
          
          OBPairData *obd = (OBPairData *) new_atom->GetData("_atom_site_label");
          
          obd->SetAttribute("_atom_site_label");
          obd->SetValue(label);/* + "_" +
                            boost::lexical_cast<string>(i) + "x" +
                            boost::lexical_cast<string>(j) + "x" +
                            boost::lexical_cast<string>(k));*/
          
          //cout << new_atom->GetData("_atom_site_label")->GetValue() << endl;
          
          vector3 ps = (*it)->GetVector() + disp;
          new_atom->SetVector(ps);
        }
      }
    }
  }
  
  if(verbose_level >= 1)
  { 
    cout << "Initial system:" << endl;
    cout << "  Chemical Formula: " << get_formula(mol_initial) << endl;
    cout << endl;
    
    cout << "Supercell system " << boost::format("(%1%x%2%x%3%)") %a %b %c << ":" << endl;
    cout << boost::format("  Size a=%1%, b=%2%, c=%3%") %super_unitcell->GetA() 
                                                        %super_unitcell->GetB() 
                                                        %super_unitcell->GetC() << endl;    
    cout << endl;
  }  
  
  return true;
}


bool d2o_main_class::read_molecule(std::string file_name)
{
  bool result = true;
  
  OBConversion obc;
  
  result = result && obc.SetInFormat(obc.FormatFromExt(file_name.c_str()));
  result = result && obc.ReadFile(&mol_initial, file_name);
  
  return result;  
}


bool d2o_main_class::set_labels_to_manual()
{
  set<string> set_lbl;
  
  for(OBAtomIterator it = mol_initial.BeginAtoms(); it != mol_initial.EndAtoms(); ++it)
  {
    string label = (*it)->GetData("_atom_site_label")->GetValue();
    set_lbl.insert(label);
  }
  
  (*manual_properties).convert_properties(set_lbl);
  
  if( (verbose_level >= 1) && (manual_properties->data().size() > 0))
  {
    cout << "Manual properties" << endl;
    cout << "Label\t|fixed\t|charge\t|popul\t|" << endl;
    for(c_man_atom_prop::data_type::const_iterator it  = manual_properties->data().begin();
                                                   it != manual_properties->data().end(); it++)
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
                             const std::vector<int> scs,
                             charge_balance cb, double tolerance_v,
                             bool merge_confs, bool calc_q_energy_v,
                             c_man_atom_prop &manual_properties_v,
                             std::string output_base_name)
{
  assert(scs.size() == 3);
  
  r_tolerance = max(tolerance_v, 1.0E-6);
  manual_properties = &manual_properties_v;  
  calc_q_energy = calc_q_energy_v;
          
  if(!read_molecule(input_file_name))
  {
    cerr << "Input file cannot be opened." << endl;
    return false;
  }
  
  if(!set_labels_to_manual())
  {
    cerr << "Manual labels cannot be set." << endl;
    return false;
  }
  
  
  if(!create_super_cell(scs[0], scs[1], scs[2]))
  {
    cerr << "Supercell cannot be created." << endl;
    return false;
  }
  
  if( !process_charges(cb) )
  {
    cerr << "Charge processing fails." << endl;
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
    cerr << "Create population error." << endl;
    return false;
  }
  
  if(verbose_level >= 1)
   show_groups_information();
  
  if(total_combinations() > 7E8)
  {
    cerr << "ERROR: Number of total combinations is too high to work with." << endl;
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

  if( calc_q_energy )  
  {
    if(!calculate_q_matrix())
    {
      cerr << "ERROR: Coulomb energy is not calculated." << endl;
      return false;
    }
    string f_name = output_base_name + "_coulomb_energy.txt";
    if(!dry_run)
    {  
      f_q_calc.open(f_name.c_str(), fstream::out);
      if(!f_q_calc.is_open())
      {
        cerr << "ERROR: File \"" << f_name << "\" is not open." << endl;
        return false;
      }
    }
  }  
  
  if(!write_files(output_base_name, dry_run, merge_confs))
  {
    cerr << "Write files error." << endl;
    return false;
  }
  
  return true;
}


d2o_main_class::~d2o_main_class()
{
}
