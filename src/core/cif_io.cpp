//
// Created by kirill on 5/12/21.
//

#include "cif_io.h"
#include <cassert>
#include <sstream>
#include <gemmi/cif.hpp>
#include <gemmi/cifdoc.hpp>
#include <gemmi/unitcell.hpp>
#include <gemmi/symmetry.hpp>
#include <gemmi/elem.hpp>
#include <gemmi/numb.hpp>
#include <gemmi/small.hpp>

#include <cryst_tools/cryst_tools.h>

cif_output::cif_output(std::ostream &output_stream, const cryst_structure_t &cs_struct,
                       const std::string &title, const std::vector<std::pair<std::string, double>> &charges) :
    cs(cs_struct), so(output_stream) {

  auto l = cs.unit_cell.lengths();
  auto a = cs.unit_cell.angles_deg();

  so
      << "# CIF file generated by supercell program. See https://orex.github.io/supercell" << std::endl
      << "# " << title << std::endl
      << "data_" << cs_struct.block_name << std::endl
      << "_chemical_name_common '" << cs_struct.chem_name << "'" << std::endl
      << "_cell_length_a " << l.x() << std::endl
      << "_cell_length_b " << l.y() << std::endl
      << "_cell_length_c " << l.z() << std::endl
      << "_cell_angle_alpha " << a.x() << std::endl
      << "_cell_angle_beta " << a.y() << std::endl
      << "_cell_angle_gamma " << a.z() << std::endl
      << "_space_group_IT_number 1" << std::endl
      << "_space_group_name_H-M_alt 'P 1'" << std::endl
      << "_space_group_name_Hall 'P 1'" << std::endl
      << "loop_" << std::endl
      << "    _space_group_symop_operation_xyz" << std::endl
      << "    x,y,z" << std::endl;
  if( !charges.empty()) {
    char buffer[64];
    so
      << "loop_ " << std::endl
      << "    _atom_type_symbol" << std::endl
      << "    _atom_type_oxidation_number" << std::endl;
    for(const auto &p : charges) {
      snprintf(buffer, sizeof(buffer), "    %-8s%+.3g\n", p.first.c_str(), p.second);
      so << buffer;
    }
  }

  so
      << "loop_" << std::endl
      << "    _atom_site_label" << std::endl
      << "    _atom_site_type_symbol" << std::endl
      << "    _atom_site_fract_x" << std::endl
      << "    _atom_site_fract_y" << std::endl
      << "    _atom_site_fract_z" << std::endl
      << "    _atom_site_occupancy" << std::endl;
}

void cif_output::add_atom(int el_num, const std::string &label, const Eigen::Vector3d &pos, double occupancy) {
  char buffer[256];

  Eigen::Vector3d v = cs.unit_cell.cell().inverse() * pos;
  for(int i = 0; i < v.size(); i++) {
    const double prec = 1e-6;
    v[i] = std::fmod(v[i], 1.0);
    if( v[i] < -0.0 )
      v[i] += 1.0;
    if( v[i] + prec >= 1.0 )
      v[i] -= 1.0;
    if( std::abs(v[i]) < prec )
      v[i] = 0.0;
  }

  int spnr = snprintf(buffer, sizeof(buffer), "    %-8s%-5s%.5f%10.5f%10.5f%8.3f\n",
           label.c_str(),  gemmi::Element(el_num).name(),
           v.x(), v.y(), v.z(), occupancy);

  if( spnr > 0 || spnr < sizeof(buffer) ) {
    so << buffer;
  } else {
    so << "ERROR data write!!!";
  }
}

static
Eigen::Vector3d norm_frac(const Eigen::Vector3d &d) {
  return Eigen::Vector3d(d.x() - round(d.x()),
                         d.y() - round(d.y()),
                         d.z() - round(d.z()));
}

static
bool verify_sg_uc_compatibility(const gemmi::SpaceGroup * sg, const gemmi::UnitCell &uc) {
  const gemmi::Mat33 & c_d = uc.orth.mat;
  const gemmi::Mat33 & c_i = uc.frac.mat;
  for( const auto &op : sg->operations() ) {
    double mult = 1.0 / op.DEN;
    gemmi::Mat33 mf(
        mult * op.rot[0][0], mult * op.rot[0][1], mult * op.rot[0][2],
        mult * op.rot[1][0], mult * op.rot[1][1], mult * op.rot[1][2],
        mult * op.rot[2][0], mult * op.rot[2][1], mult * op.rot[2][2] );
    gemmi::Mat33 mt = c_d.multiply(mf).multiply(c_i);
    if( !mt.transpose().approx(mt.inverse(), 1e-3) )
      return false;
  }
  return true;
}

const gemmi::SpaceGroup * resolve_sg_by_uc(const gemmi::SpaceGroup * sg, const gemmi::UnitCell &uc) {
  if( sg == nullptr )
    return nullptr;
  if( verify_sg_uc_compatibility(sg, uc) )
    return sg;
  for (const auto& st : gemmi::spacegroup_tables::main) {
    if( st.number != sg->number )
      continue;
    if( verify_sg_uc_compatibility(&st, uc) )
      return &st;
  }
  return nullptr;
}

static
bool read_cif_file_gemmi(const std::string &file_name, cryst_structure_t& result, std::string &msg) {
  msg.clear();
  using namespace gemmi;
  using cif::as_number;
  using cif::as_string;

  auto doc = cif::read_file(file_name);

  cif::Block& block = doc.sole_block();
  result.block_name = block.name;

  {
    cif::Table cm_name = block.find({"_chemical_name_common"});
    if( cm_name.ok() )
      result.chem_name = as_string(cm_name.one()[0]);
  }
  // unit cell and symmetry
  gemmi::UnitCell uc_supercell;
  {
    cif::Table cell = block.find("_cell_",
                                 {"length_a", "length_b", "length_c",
                                  "angle_alpha", "angle_beta", "angle_gamma"});
    if (cell.ok()) {
      auto c = cell.one();
      if (!cif::is_null(c[0]) && !cif::is_null(c[1]) && !cif::is_null(c[2])) {
        uc_supercell.set(as_number(c[0]),
                         as_number(c[1]),
                         as_number(c[2]),
                   as_number(c[3]), as_number(c[4]), as_number(c[5]));
      } else {
        msg += "ERROR: Some of cell values are wrong";
        return false;
      }
    }
  }

  std::vector<gemmi::Op> ops;
  /*for (const auto &t : {"_space_group_symop_operation_xyz",
                        "_symmetry_equiv_pos_as_xyz"}) {
    cif::Table op_table = block.find("", {t});
    if (op_table.ok()) {
      for (auto row : op_table) {
        ops.emplace_back(gemmi::parse_triplet(as_string(row[0])));
      }
      break;
    }
  }*/

  const gemmi::SpaceGroup * sg = nullptr;
  {
    for (const char *tag : {"_space_group_name_H-M_alt",
                            "_symmetry_space_group_name_H-M",
                            "_space_group_IT_number"}) {
      if (const std::string *val = block.find_value(tag)) {
        sg = gemmi::find_spacegroup_by_name(as_string(*val), uc_supercell.alpha, uc_supercell.gamma);
        if (sg == nullptr) {
          msg += "ERROR: Group name '" + as_string(*val) + "' in '" + tag + "' is wrong\n";
          return false;
        } else
          break;
      }
    }
    if( sg == nullptr ) {
      msg += "ERROR: H-M tags is not found.\n";
      return false;
    }
    sg = resolve_sg_by_uc(sg, uc_supercell);
    if( sg == nullptr ) {
      msg += "ERROR: Spacegroup is not consistence with unitcell.\n";
      return false;
    }
    if (const std::string *val = block.find_value("_space_group_IT_number") ) {
      if( as_number(*val) != sg->number ) {
        msg += "WARNING: Space group number (" + as_string(*val) +  ") is different from H-M ("
            + std::to_string(sg->number) + ").\n";
      }
    }
    if (const std::string *val = block.find_value("_space_group_name_Hall") ) {
      if( as_string(*val) != sg->hall ) {
        msg += "WARNING: Space group Hall name (" + as_string(*val) +  ") is different from H-M ("
            + sg->hall + ").\n";
      }
    }
  }
  if( !ops.empty() ) {
    msg += "INFO: Using direct symmetries.\n";
  } else {
    msg += "INFO: Using symmetries from space group.\n";
    for(const auto &x : sg->operations())
      ops.emplace_back(x);
  }

  Eigen::Matrix3d eg_s;
  const auto &m = uc_supercell.orth.mat;
  eg_s << m.a[0][0], m.a[0][1], m.a[0][2],
          m.a[1][0], m.a[1][1], m.a[1][2],
          m.a[2][0], m.a[2][1], m.a[2][2];
  result.unit_cell.set(eg_s);

  enum { kTypeName, kChargeValue };
  cif::Table charge_table = block.find("_atom_type_",
                                     {"symbol",
                                      "oxidation_number"});
  std::map<std::string, double> m_charges;
  for (auto row : charge_table) {
    m_charges[as_string(row[kTypeName])] = as_number(row[kChargeValue]);
  }

  enum { kLabel, kSymbol, kX, kY, kZ, kSymmMult, kOcc };
  cif::Table atom_table = block.find("_atom_site_",
                                     {"label",
                                      "?type_symbol",
                                      "fract_x",
                                      "fract_y",
                                      "fract_z",
                                      "?symmetry_multiplicity",
                                      "?occupancy"});
  result.atoms.clear();
  for (auto row : atom_table) {
    atom_t site;
    site.label = as_string(row[kLabel]);
    std::string type_symbol = row.has(kSymbol) ? as_string(row[kSymbol]) : site.label;

    site.fract_pos = Eigen::Vector3d(as_number(row[kX]), as_number(row[kY]), as_number(row[kZ]));
    site.occupancy = row.has(kOcc) ? as_number(row[kOcc], 1.0) : 1.0;
    struct { Element element = El::X; signed char charge = 0; } spl_out;
    split_element_and_charge(type_symbol, &spl_out);
    site.el_num = spl_out.element.atomic_number();
    site.charge = spl_out.charge;
    for(const auto &n : {type_symbol, site.label}) {
      auto it = m_charges.find(n);
      if( it != m_charges.end() )
        site.charge = it->second;
    }
    std::vector<Eigen::Vector3d> smx;
    for (const auto& op : ops) {
      double mult = 1.0 / op.DEN;
      Eigen::Matrix3d rot;
      rot <<
                mult * op.rot[0][0], mult * op.rot[0][1], mult * op.rot[0][2],
                mult * op.rot[1][0], mult * op.rot[1][1], mult * op.rot[1][2],
                mult * op.rot[2][0], mult * op.rot[2][1], mult * op.rot[2][2];
      Eigen::Vector3d tran(mult * op.tran[0], mult * op.tran[1], mult * op.tran[2]);
      Eigen::Vector3d np = rot * site.fract_pos + tran;
      bool dup = false;
      for(int i = 0; i < smx.size() && !dup; i++) {
        dup = dup || norm_frac( smx[i] - np).norm() < 1E-3;
      }
      if( !dup ) {
        smx.push_back(np);
      }
    }
    if (row.has(kSymmMult) && as_number(row[kSymmMult]) != smx.size() ) {
      msg += "WARNING: Site '" + site.label + "' has unexpected '_atom_site_symmetry_multiplicity' \n";
      msg += "         CIF: " + as_string(row[kSymmMult]) + " and actual: " + std::to_string(smx.size()) + "\n";
    }
    for(const auto &a : smx) {
      result.atoms.emplace_back(site);
      result.atoms.back().fract_pos = a;
    }
  }
  return true;
}

bool read_cif_file(const std::string &file_name, cryst_structure_t& result, std::string &msg) {
  return read_cif_file_gemmi(file_name, result, msg);
}
