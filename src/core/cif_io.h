//
// Created by kirill on 5/12/21.
//

#ifndef SUPERCELL_SRC_CORE_CIF_IO_H_
#define SUPERCELL_SRC_CORE_CIF_IO_H_

#include <ostream>
#include <sstream>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "common_types.h"

class cif_output {
 private:
  cryst_structure_t cs;
  std::ostream &so;
 public:
  cif_output(std::ostream &output_stream, const cryst_structure_t &cs_struct, const std::string &title);
  void add_atom(int el_num, const std::string &label, const Eigen::Vector3d &pos, double occupancy);
};

bool read_cif_file(const std::string &file_name, cryst_structure_t& result, std::string &msg);

#endif //SUPERCELL_SRC_CORE_CIF_IO_H_
