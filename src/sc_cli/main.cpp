/* 
 * File:   main.cpp
 * Author: kirill
 *
 * Created on July 18, 2013, 1:49 PM
 */

#include <cstdlib>
#include <random>
#include <iostream>
#include <string>
#include <set>

#include "boost/program_options.hpp" 
#include "boost/filesystem.hpp"
#include "parse_d2o_input.h"
#include "d2o_main_class.h"
#include <eigen3/Eigen/Dense>

#include "common_types.h"
 

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
    string output_tar_file;
    bool dry_run = false;
    bool merge_confs = false;
    bool calc_q = false;
    bool store_q_all = false;
    //double memory_limit;
    int verb_level;
    string cell_size_str;
    string charge_balance_str;
    
    double pos_tol;
    
    vector<string> manual_properties;
    vector<string> structure_sampling;
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
      ("dry-run,d", "Shows information, but does not generate structures. Always use when trying a new system.")
      ("cell-size,s", po::value<std::string>(&cell_size_str)->default_value("1x1x1"), 
                  "Supercell size. Example: -s 2x2x5. Default is 1x1x1.")
      ("charge-balance,c", po::value<std::string>(&charge_balance_str)->default_value("try"), 
                           (cb_names::get_name(charge_balance::cb_no)     + " - no charge balancing.\n" +
                            cb_names::get_name(charge_balance::cb_try)    + " - Try to charge balance system, " +
                                                           "if initial system is not charged.\n" +
                            cb_names::get_name(charge_balance::cb_yes)  + " - Charge balance the system. ").c_str())
      ("property,p", po::value<vector<string> >(&manual_properties), 
                    (string("Set properties of atoms by labels. ") +
                            "For detailed description see manual.").c_str())
      ("tolerance,t", po::value<double>(&pos_tol)->default_value(0.75), 
                      "Skip structures with atoms closer than <arg> Angstrom.")
      ("merge-symmetric,m", "Merge output structures that are equivalent by symmetry.")
      ("coulomb-energy,q", "Calculate Coulomb energies of output structures.")
      ("coulomb-store,g", "Store electrostatic energy for all structures.")
      ("store-structures,n", po::value<vector<string> >(&structure_sampling),
                            "Generate structures according to certain criteria. See manual for details.")
      ("random-seed", po::value<std::random_device::result_type>(),
                            "Set seed for random structure sampling.")
      ("archive,a", po::value<std::string>(&output_tar_file)->default_value(""),
                   (string("A target archive file for all output files. If empty (default) no packing will be performed. ") +
#ifdef LIBARCHIVE_ENABLED
                           "The option is enabled."   
#else
                           "The option is disabled." 
#endif    
                     ).c_str() )
      ("output,o", po::value<std::string>(&output_file)->default_value("supercell"), 
                   "Output file base name. The extension will be cif. The structure multiplicity will be added."); 
 
    po::variables_map vm; 
 
    try 
    { 
      po::store(po::command_line_parser(argc, argv).options(desc).run(), vm); 
      // throws on error 
 
      /** --help option 
       */
      if (vm.count("help")) {
        std::cout << "----------------------------------------------------" << endl;
        std::cout << "Program generating supercells from crystal structures" << endl;
        std::cout << "with partial occupancies and/or mixed-composition sites." << endl;
        std::cout << "----------------------------------------------------" << endl;
        std::cout << "supercell -i <INPUT_FILE> <OPTIONS>" << endl;
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
      cout << "-----------------------------------------------------" << endl;
      cout << "-            Supercell program (v2.1.1)             -" << endl;
      cout << "-         https://orex.github.io/supercell/         -" << endl;
      cout << "-----------------------------------------------------" << endl;
      cout << "-      Authors:   * Kirill Okhotnikov               -" << endl;
      cout << "-                  (kirill.okhotnikov@gmail.com)    -" << endl;
      cout << "-                 * Sylvian Cadars                  -" << endl;
      cout << "-                  (sylvian.cadars@cnrs-imn.fr)     -" << endl;
      cout << "-                 * Thibault Charpentier            -" << endl;
      cout << "-                  (Thibault.Charpentier@cea.fr)    -" << endl;
      cout << "-----------------------------------------------------" << endl;
      cout << "-  please cite:                                     -" << endl;
      cout << "-    K. Okhotnikov, T. Charpentier and S. Cadars    -" << endl;
      cout << "-    J. Cheminform. 8 (2016) 17 – 33.               -" << endl;
      cout << "-----------------------------------------------------" << endl;
      cout << endl;
    }

    dry_run = vm.count("dry-run") > 0;
    merge_confs = vm.count("merge-symmetric") > 0;
    calc_q = vm.count("coulomb-energy") > 0;
    store_q_all = vm.count("coulomb-store") > 0 && calc_q;

    std::random_device::result_type rnd_seed
        = vm.count("random-seed") == 0 ? std::random_device()()
                                       : vm["random-seed"].as<std::random_device::result_type>();

    if( verb_level > 0 ) {
      std::cout << "Command line: ";
      for(int i = 0; i < argc; i++)
        std::cout << argv[i] << " ";
      std::cout << std::endl;
    }

    if( verb_level > 0 ) {
      cout << "Random SEED: " << rnd_seed << endl;
    }

    if(!parse_d2o_input::get_supercell_size(cell_size_str, supercell_mult))
    {
      cerr << "Wrong supercell format input." << endl;
      return ERROR_IN_COMMAND_LINE;
    }
    
    #ifndef LIBARCHIVE_ENABLED
    if( ! output_tar_file.empty() )
    {
      cerr << "Cannot pack output. LibArchive is not enabled." << endl;
      return ERROR_IN_COMMAND_LINE;
    }  
    #endif

    c_man_atom_prop_cli m_prop;
    m_prop.set_verbose(verb_level);
    string param_error;
    
    if(!m_prop.parse_input(manual_properties, param_error))
    {
      cerr << "Error in manual property input parameter " << param_error << endl;
      return ERROR_IN_COMMAND_LINE;
    }
    
    c_struct_sel sampl_prop;

    if(!parse_d2o_input::parse_sel_input(structure_sampling, sampl_prop, param_error))
    {
      cerr << "Error in sampling input parameter " << param_error << endl;
      return ERROR_IN_COMMAND_LINE;
    }

    if( !merge_confs && sampl_prop.str_weight_limit() > 0 ) {
      cerr << "Weight can be enabled only for merge configuration run." << endl;
      return ERROR_IN_COMMAND_LINE;
    }

    if (!parse_d2o_input::check_adjust_output_name(
            output_file, !output_tar_file.empty(), param_error)) {
      cerr << "Error in output suffix: " << param_error << endl;
      return ERROR_IN_COMMAND_LINE;
    }

    charge_balance cb;
    
    if(!parse_d2o_input::get_charge_balance(charge_balance_str, cb))
    {
      cerr << "Wrong Charge Balance input." << endl;
      return ERROR_IN_COMMAND_LINE;
    }
    
    d2o_main_class mc(rnd_seed);
    
    mc.set_verbosity(verb_level);
    
    bool processed =
        mc.process(input_file, dry_run, supercell_mult,
                   cb, pos_tol, merge_confs, calc_q, store_q_all,
                   m_prop, sampl_prop, output_file,
                   output_tar_file);
    
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