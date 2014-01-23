/* 
 * File:   main.cpp
 * Author: kirill
 *
 * Created on July 18, 2013, 1:49 PM
 */

#include <cstdlib>

#include "boost/program_options.hpp" 
#include "boost/filesystem.hpp"
#include "parse_d2o_input.h"
#include "d2o_main_class.h"
#include <eigen3/Eigen/Dense>

#include "common_types.h"
 
#include <iostream> 
#include <string> 
#include <set>

namespace 
{ 
  const size_t ERROR_IN_COMMAND_LINE = 1; 
  const size_t SUCCESS = 0; 
  const size_t ERROR_UNHANDLED_EXCEPTION = 2; 
  const size_t ERROR_PROCESS_EXECUTION = 3; 
 
} // namespace 

using namespace std;

/*
 * 
 */

int main(int argc, char** argv)
{
 
  try 
  { 
    //std::string appName = boost::filesystem::basename(argv[0]); 

    string input_file;
    string output_file;
    bool dry_run = false;
    bool merge_confs = false;
    bool calc_q = false;
    //double memory_limit;
    int verb_level;
    string cell_size_str;
    string charge_balance_str;
    
    double pos_tol;
    
    vector<string> manual_properties;
    vector<int> supercell_mult;
    supercell_mult.resize(3);
    
 
    /** Define and parse the program options 
     */ 
    namespace po = boost::program_options; 
    po::options_description desc("Options"); 
    desc.add_options() 
      ("help,h", "Print help messages") 
      ("verbose,v", po::value<int>(&verb_level)->default_value(1), 
                  "output data with verbosity level. Default is 1. Suggested for regular users.") 
      ("input,i", po::value<std::string>(&input_file)->required(), 
                  "Input structure")
      ("dry-run,d", "Show information, but not generate structures")
      ("cell-size,s", po::value<std::string>(&cell_size_str)->default_value("1x1x1"), 
                  "Supercell size. Example: -s 2x2x5")
      ("charge-balance,c", po::value<std::string>(&charge_balance_str)->default_value("try"), 
                           (cb_names::get_name(cb_no)     + " - no charge balancing.\n" +
                            cb_names::get_name(cb_try)    + " - Try to charge balance system, " +
                                                           "if initial system is not charged.\n" +
                            cb_names::get_name(cb_input)  + " - Charge balance the system. ").c_str())
      ("property,p", po::value<vector<string> >(&manual_properties), 
                    (string("Set properties of atoms by labels. ") +
                            "For detail description see manual.").c_str())
      ("tolerance,t", po::value<double>(&pos_tol)->default_value(0.75), 
                      "Skip structures with atoms closer than arg Angstrom.")
      ("merge-symmetric,m", "Merge output equivalent (by symmetry) structures.")
      ("coulomb-energy,q", "Calculate Coulomb energy of output structures.")
      ("output,o", po::value<std::string>(&output_file)->default_value("supercell"), 
                   "Output file name base. The extension will be cif. The multiplicity of structure will be added."); 
 
    po::variables_map vm; 
 
    try 
    { 
      po::store(po::command_line_parser(argc, argv).options(desc).run(), vm); 
      // throws on error 
 
      /** --help option 
       */ 
      if ( vm.count("help")  ) 
      { 
        std::cout << "Disordered cell to ordered supercell program." << endl;
        std::cout << "The values in parenthesis are default values." << endl;
        std::cout << desc << endl;
        
        return SUCCESS; 
      } 
 
      po::notify(vm); // throws on error, so do after help in case 
                      // there are any problems 
    } 
    
    catch(boost::program_options::required_option& e) 
    { 
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl; 
      return ERROR_IN_COMMAND_LINE; 
    } 
    catch(boost::program_options::error& e) 
    { 
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl; 
      return ERROR_IN_COMMAND_LINE; 
    } 
    
    //Title
    if(verb_level > 0)
    {  
      cout << "-------------------------------" << endl;
      cout << "-----------Title---------------" << endl;
      cout << "-------------------------------" << endl;
    }  
    
    dry_run = vm.count("dry-run") > 0;
    merge_confs = vm.count("merge-symmetric") > 0;
    calc_q = vm.count("coloumb-energy");
    
    if(!parse_d2o_input::get_supercell_size(cell_size_str, supercell_mult))
    {
      cerr << "Wrong supercell format input." << endl;
      return ERROR_IN_COMMAND_LINE;
    }

    c_man_atom_prop_cli m_prop;
    m_prop.set_verbose(verb_level);
    string param_error;
    
    if(!m_prop.parse_input(manual_properties, param_error))
    {
      cerr << "Error in input parameter " << param_error << endl;
      return false;
    }  
    
    charge_balance cb;
    
    if(!parse_d2o_input::get_charge_balance(charge_balance_str, cb))
    {
      cerr << "Wrong Charge Balance input." << endl;
      return ERROR_IN_COMMAND_LINE;
    }
    
    d2o_main_class mc;
    
    mc.set_verbosity(verb_level);
    
    bool processed = 
    mc.process(input_file, dry_run, supercell_mult, 
               cb, pos_tol, merge_confs, calc_q, 
               m_prop, output_file);
    
    if(!processed)
      return ERROR_PROCESS_EXECUTION;
  } 
  catch(std::exception& e) 
  { 
    std::cerr << "Unhandled Exception reached the top of main: " 
              << e.what() << ", application will now exit" << std::endl; 
    return ERROR_UNHANDLED_EXCEPTION; 
 
  } 
 
  return SUCCESS; 
  
}

/*
  CB_combination::num_map mp;
  CB_combination cbc;
  
  mp[0] = 2;
  mp[1] = 2;
  mp[2] = 2;  
  mp[3] = 2;    
  
  vector<int> fc;

  cout << cbc.create_first_combination(mp, fc) << endl;
  
  int i = 1;
  while(cbc.get_next_combination(fc))
  {  
    for(int j = 0; j < fc.size(); j++)
      cout << fc[j] << " - ";
    cout << endl;    
    i++;
  }

  //sleep(1);  
  
  cout << i << endl;
  
  return 0; 
 
   c_man_atom_prop_cli cit;
  
  cit.regex_test(" c= 0 c = 2 fixed p=4 charge=1.2 ");
  cit.regex_test("c=-1.2");
  //cit.parse_input_item("r(La(    M):p");
  
  return 0;

 
 */