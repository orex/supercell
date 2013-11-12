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

#include "obabel/obabel_utils.h"


using namespace OpenBabel;
using namespace std;

map_comp_item::map_comp_item(std::vector<lbl_order> &lbl)
{ 
  labels_order = &lbl;
  mol = new OBMol;
  dsts = new std::map<std::string, std::vector<float> > ;
  uc = NULL;
  mult = 0;
}

OpenBabel::OBUnitCell * map_comp_item::unitcell()
{
  assert(mol->HasData(OBGenericDataType::UnitCell));
  if( uc == NULL)
    uc = dynamic_cast<OBUnitCell *>(mol->GetData(OBGenericDataType::UnitCell));
  
  return uc;
}


std::vector<float> map_comp_item::get_lengths_by_labels(const lbl_order &lbl)
{
  std::vector<float> result;
  
  string map_key = lbl.lbl1 + "  -  " + lbl.lbl2;
  
  assert(dsts != NULL);
  assert(labels_order != NULL);
  
  if(dsts->count(map_key) > 0)
  {  
    result = (*dsts)[map_key];
    return result;
  }
  
  bool homo_dist = lbl.lbl1 == lbl.lbl2;
  
  for(int i = 0; i < mol->NumAtoms(); i++)
  {
    OBAtom * atom_i = mol->GetAtom(i + 1);
    string label_i = atom_i->GetData("original_label")->GetValue();
    
    if(label_i != lbl.lbl1)
      continue;
    
    for(int j = homo_dist ? i + 1: 0; j < mol->NumAtoms(); j++)
    {
      OBAtom * atom_j = mol->GetAtom(j + 1);
      string label_j = atom_j->GetData("original_label")->GetValue();
    
      if(label_j != lbl.lbl2)
        continue;
      
      vector3 dist = atom_i->GetVector() - atom_j->GetVector();
      dist = get_minimal_distance(dist, unitcell());
      
      result.push_back(dist.length());
    }  
  }

  std::sort(result.begin(), result.end());
  (*dsts)[map_key] = result;
  
  return result;
}

int map_comp_item::compare_distances(const std::vector<float> &f1, const std::vector<float> &f2, float tol)
{
  assert(f1.size() == f2.size());
  
  int result = 0;
  
  for(int i = 0; i < f1.size(); i++)
  {
    result = f1[i] < f2[i] ? -1 : 1;

    if( abs(f1[i] - f2[i]) < tol)
      result = 0;
    
    if( result != 0 ) break;
  }

  return result;  
}

int map_comp_item::comp( map_comp_item &r1, map_comp_item &r2 )
{
  int cmp;
  assert(r1.labels_order->size() == r2.labels_order->size());
  for(int i = 0; i < r1.labels_order->size(); i++)
  {
    std::vector<float> f1, f2;
    f1 = r1.get_lengths_by_labels(r1.labels_order->at(i));
    f2 = r2.get_lengths_by_labels(r2.labels_order->at(i));
    
    cmp = compare_distances(f1, f2, r1.labels_order->at(i).tolerance);
    if(cmp != 0) break;
  }
  
  return cmp;
}

d2o_main_class::d2o_main_class()
{
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
  int result = 1;

  if( !_fixed )
  {
    std::map<int,int> nm;

    int sum = 0;  
    for(int i = 0; i < items.size(); i++)
    {  
      nm[i] = items[i].num_of_atoms_sc;
      sum  += items[i].num_of_atoms_sc;
    }  

    if(sum <= number_of_sites)
    {  
      nm[-1] = number_of_sites - sum;    
      result = num_combinations(nm);
    }
    else
      result = 0;
  }
  
  return result;
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

bool d2o_main_class::add_confs_to_mol(OpenBabel::OBMol *cmol, const std::map<int, std::vector<int> > &ppc)
{
  map<int, int> curr_pos;
  
  FOR_ATOMS_OF_MOL(a, mol_supercell)
  {
    int group_index = dynamic_cast<OBPairInteger *>(a->GetData("group_number"))->GetGenericValue();
    
    assert( ( group_index >= 0 ) && ( group_index < occup_groups.size() ) );
    
    c_occup_group &curr_group = occup_groups[group_index];
    
    for(int i  = 0; i < curr_group.items.size(); i++)
    {
      double occup_value;
      
      if( curr_group.is_fixed() )
      {
        assert( ppc.count(group_index) == 0 );
        occup_value = double(curr_group.items[i].num_of_atoms_sc) /
                      double(curr_group.number_of_sites);
      }
      else
      {
        assert( ppc.count(group_index) > 0);        
        assert( ppc.at(group_index).size() > 0);
        assert( curr_pos[group_index] < ppc.at(group_index).size() );
        if( ppc.at(group_index)[curr_pos[group_index]] != i )
          continue;
        occup_value = 1.0;
      }  
            
      //Adding atom
      OBAtom atm;
      atm.Duplicate( curr_group.items[i].obp );

      if( atm.HasData("_atom_site_occupancy") )
        dynamic_cast<OBPairFloatingPoint *>(atm.GetData("_atom_site_occupancy"))->SetValue(occup_value);
      else
      {
        OBPairFloatingPoint * obo = new OBPairFloatingPoint();
        obo->SetAttribute("_atom_site_occupancy");
        obo->SetValue(occup_value);
        atm.SetData(obo);
      }  
  
      atm.SetVector(a->GetVector());
      cmol->AddAtom(atm, true);
    }  
    
    curr_pos[group_index]++;
  }
 
  return true;
}

bool d2o_main_class::add_to_list(std::list<map_comp_item *> &mpis, map_comp_item * mpi)
{
  typedef std::list<map_comp_item *>::iterator c_it;
  
  c_it min_it = mpis.begin();
  c_it max_it = mpis.end();
  c_it curr_pos = mpis.end();
  
  bool fnd = false;
  while(min_it != max_it)
  {
    curr_pos = min_it;
    advance(curr_pos, (std::distance(min_it, max_it) / 2));
 
    int cmp = map_comp_item::comp(*(*curr_pos), *mpi);
    
    if(cmp == 0)
    {
      fnd = true;
      break;
    }
    else if(cmp > 0)
    {
      max_it = curr_pos;
    }
    else
    {  
      curr_pos++;
      min_it = curr_pos;
    }  
  }
  
  if(!fnd)
    curr_pos = mpis.insert(curr_pos, mpi);
  else
    delete mpi;
    
  (*curr_pos)->inc_mult();  
  
  return !fnd;
}

bool d2o_main_class::write_files(std::string output_base_name, double n_store, bool dry_run, bool merge_confs)
{
  if(dry_run && (!merge_confs) )
    return true;
  
  std::vector<lbl_order> lo;
  lo = set_lbl_order();
  
  list<map_comp_item *> mpis;
      
  typedef map<int, vector<int> > t_map_comb;
  t_map_comb cur_combs;

  for(int i = 0; i < occup_groups.size(); i++)
  {  
    if( occup_groups[i].is_fixed() )
      continue;
    
    map<int, int> mvc;

    for(int j = 0; j < occup_groups[i].items.size(); j++)
      mvc[j] = occup_groups[i].items[j].num_of_atoms_sc;
    
    mvc[-1] = occup_groups[i].number_of_sites - 
                      occup_groups[i].get_total_num_occup_sites();
    
    cur_combs[i] = create_start_combination(mvc);
  }  
  
  int index = 0;
  bool done;
  
  double tot_comb = total_combinations();
  
  #ifdef USE_FIXED_N_RND
  set<int> rnd_index;
  if(n_store > 0 )
  {
    vector<int> rnd_index_v = get_random_numbers(n_store, tot_comb - 1);
    rnd_index.insert(rnd_index_v.begin(), rnd_index_v.end());
  }
  #else
  boost::mt19937 rnd_gen = create_rnd_gen();
  #endif

  if(!dry_run)
  {
    string del_command = "rm -f " + output_base_name + "*.cif";
    int rc = system(del_command.c_str());
    if( (verbose_level >= 2) && (rc == 0) )
      cout << "Output files was deleted successfully" << endl;
  }
    
    
  int total_index;
  do
  {
    #ifdef USE_FIXED_N_RND
    if( (n_store <= 0) || (rnd_index.count(total_index) > 0) )    
    #else
    if( (n_store <= 0) || get_rnd_value_in_interval(rnd_gen, 0, tot_comb) <= n_store )
    #endif
    {  
      map_comp_item * mpi = new map_comp_item(lo);
      init_atom_change_mol(mpi->mol);

      add_confs_to_mol(mpi->mol, cur_combs);

      if( !merge_confs )
      {  
        if(!dry_run)
        {  
          OBConversion obc;
        
          string index_str = get_index_str(index, tot_comb - 1);

          obc.SetOutFormat("cif");
          obc.WriteFile(mpi->mol, 
                output_base_name + "_ind" + index_str + ".cif");
          
          if( (index % 500 == 0) && (index != 0) && (verbose_level >= 2))
            cout << "Processed " << index << " configuration..." << endl;
        }
        
        delete mpi;
      }
      else if( merge_confs )
      {
        add_to_list(mpis, mpi);
        if( (index % 500 == 0) && (index != 0) && (verbose_level >= 2))
        {  
          cout << "Processed " << index << " configuration merged to " << 
                  mpis.size() << " configuration..."  << endl;
        }  
      }
      index++;    
    }  

    done = true;    
    for(t_map_comb::iterator it  = cur_combs.begin(); 
                             it != cur_combs.end(); it++)
    {
      if( std::next_permutation(it->second.begin(), it->second.end()) )
      {
        done = false;
        break;
      }  
    }
    total_index++;
  }while(!done);
  
  if(merge_confs && (verbose_level >= 1) )
    cout << "Combinations after merge: " << mpis.size() << endl;
  
  if( (!dry_run) && merge_confs )
  {
    int c_num = 0;
    for(list<map_comp_item *>::iterator it = mpis.begin(); it != mpis.end(); it++)
    {  
      OBConversion obc;
      
      string index_str = get_index_str(c_num, mpis.size() - 1);      
      
      obc.SetOutFormat("cif");
      obc.WriteFile((*it)->mol, 
              output_base_name + "_ind" + index_str +
              "w_" + boost::lexical_cast<string>((*it)->get_mult())
              + ".cif");
      delete (*it);
      c_num++;
    }  
  }  
  
  return true;
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

std::vector<lbl_order> d2o_main_class::set_lbl_order()
{
  const double min_tol = 0.01;
  std::vector<lbl_order> result;
  
  for(int i = 0; i < occup_groups.size(); i++)
  {
    // tolerance depricated
    //double tol1 = occup_groups[i].max_dis_within_group + min_tol;
    int64_t m1 = occup_groups[i].get_number_of_combinations();
    for(int j = 0; j < occup_groups[i].items.size(); j++)
    {
      string lbl1 = occup_groups[i].items[j].label;
      for(int k = i; k < occup_groups.size(); k++)
      {
        // tolerance depricated
        //double tol2 = occup_groups[k].max_dis_within_group + min_tol;
        int64_t m2 = occup_groups[k].get_number_of_combinations();
        
        int begp = k == i ? j : 0;
          
        for(int l = begp; l < occup_groups[k].items.size(); l++)
        {
          string lbl2 = occup_groups[k].items[l].label;
          if( m1 * m2 > 1)
          {  
            lbl_order tp(lbl1, lbl2, min_tol, m1 * m2);
            result.push_back(tp);
          }  
        }  
      }  
    }  
  }  
  
  std::sort(result.begin(), result.end());
  std::reverse(result.begin(), result.end());
  
  return result;
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
      rd.min_value = curr_group.number_of_sites;
      rd.max_value = curr_group.number_of_sites;
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
          correct_rms_range(curr_group.number_of_sites, 
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
        overoccup = overoccup || (occup_groups[i].number_of_sites < ocp_t[i]);
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
          underoccup = (occup_groups[i].number_of_sites > ocp_t[i]) &&
                       (abs(1.0 - occup_groups[i].get_total_occup_input()) < occup_tol);
        }          
        
        if((verbose_level >= 5) && underoccup)
          cout << "Under occup" << endl;

        if( underoccup ) break;
      }
      
      if(verbose_level >= 5)
      {
        for(int i = 0; i < occup_groups.size(); i++)
          cout << "Curr occup is " << ocp_t[i] << ". Total is " << occup_groups[i].number_of_sites << endl;
      }  
      
      if(verbose_level >= 5)
        cout << "charge: " << charge << endl;

      
      if( (abs(charge) < charge_tol) && (!overoccup) && (!underoccup))
      {  
        double rms_curr = 0;
        //calculate RMS
        for(int i = 0; i < rc.size(); i++)
        {  
          double group_sites = occup_groups[rc[i].group_index].number_of_sites;
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


bool d2o_main_class::show_groups_information()
{
  cout << "Chemical formula of the supercell: " << get_formula_by_groups() << endl;
  cout << "Minimal distance between groups: " << min_dist_between_groups << endl;
  
  for(int i = 0; i < occup_groups.size(); i++)
  {
    cout << " Group " << i << " has " 
         << occup_groups[i].number_of_sites 
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
                  double(occup_groups[i].number_of_sites) );
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
    
    if( (occup_groups[i].get_total_num_occup_sites() < occup_groups[i].number_of_sites) &&
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
                             bool merge_confs,
                             c_man_atom_prop &manual_properties_v,
                             double n_store,
                             std::string output_base_name)
{
  assert(scs.size() == 3);
  
  r_tolerance = max(tolerance_v, 1.0E-6);
  manual_properties = &manual_properties_v;  
          
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
  
  if(total_combinations() > 5E8)
  {
    cerr << "ERROR: Number of total combinations is too high to work with." << endl;
    return false;
  }
    
  
  if(!write_files(output_base_name, n_store, dry_run, merge_confs))
  {
    cerr << "Write files error." << endl;
    return false;
  }
  
  return true;
}


d2o_main_class::~d2o_main_class()
{
}

