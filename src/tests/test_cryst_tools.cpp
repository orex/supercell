#define BOOST_TEST_MODULE cryst_tool_test

#include <boost/test/unit_test.hpp>

#include <filesystem>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <gemmi/cif.hpp>
#include <gemmi/cifdoc.hpp>
#include <gemmi/numb.hpp>
#include <gemmi/symmetry.hpp>
#include <gemmi/unitcell.hpp>

#include "cif_io.h"
#include "common_types.h"
#include "cryst_tools/cryst_tools.h"

#ifndef CT_DATA_DIR
#define CT_DATA_DIR "."
#endif

using namespace std;
using namespace cryst_tools;
using namespace Eigen;

namespace {

// Resolve the space group declared in the CIF, mirroring the (read-only) logic
// of read_cif_file_gemmi() in cif_io.cpp.
const gemmi::SpaceGroup *declared_spacegroup(const std::string &file_name)
{
  using namespace gemmi;
  using cif::as_number;
  using cif::as_string;

  auto doc = cif::read_file(file_name);
  cif::Block &block = doc.sole_block();

  UnitCell uc;
  cif::Table cell = block.find("_cell_",
                               {"length_a", "length_b", "length_c",
                                "angle_alpha", "angle_beta", "angle_gamma"});
  if (!cell.ok())
    return nullptr;
  auto c = cell.one();
  uc.set(as_number(c[0]), as_number(c[1]), as_number(c[2]),
         as_number(c[3]), as_number(c[4]), as_number(c[5]));

  for (const char *tag : {"_space_group_name_H-M_alt",
                          "_symmetry_space_group_name_H-M",
                          "_space_group_IT_number"}) {
    if (const std::string *val = block.find_value(tag)) {
      const SpaceGroup *sg =
          find_spacegroup_by_name(as_string(*val), uc.alpha, uc.gamma);
      if (sg)
        return sg;
    }
  }
  return nullptr;
}

Affine3d gemmi_op_to_affine(const gemmi::Op &op)
{
  const double mult = 1.0 / op.DEN;
  Matrix3d rot;
  rot << mult * op.rot[0][0], mult * op.rot[0][1], mult * op.rot[0][2],
         mult * op.rot[1][0], mult * op.rot[1][1], mult * op.rot[1][2],
         mult * op.rot[2][0], mult * op.rot[2][1], mult * op.rot[2][2];
  Affine3d a;
  a.linear() = rot;
  a.translation() = Vector3d(mult * op.tran[0], mult * op.tran[1], mult * op.tran[2]);
  return a;
}

vc_sets group_by_element(const cryst_structure_t &cs)
{
  std::map<int, std::vector<Vector3d>> by_el;
  for (const auto &a : cs.atoms)
    by_el[a.el_num].push_back(a.fract_pos);
  vc_sets out;
  out.reserve(by_el.size());
  for (auto &kv : by_el)
    out.push_back(std::move(kv.second));
  return out;
}

} // namespace


BOOST_AUTO_TEST_SUITE(CrystToolsTest)

BOOST_AUTO_TEST_CASE(Test_first)
{
  namespace fs = std::filesystem;
  fs::path test_data_dir(CT_DATA_DIR);

  BOOST_REQUIRE(fs::exists(test_data_dir));
  BOOST_REQUIRE(fs::is_directory(test_data_dir));

  for (const auto &entry : fs::directory_iterator(test_data_dir)) {
    if (!entry.is_regular_file() || entry.path().extension() != ".cif")
      continue;

    cout << "Checking file: " << entry.path() << endl;
    const std::string path = entry.path().string();

    // P1-expanded atoms + cell, via the production CIF reader.
    cryst_structure_t cs;
    std::string msg;
    BOOST_REQUIRE_MESSAGE(read_cif_file(path, cs, msg),
                          "read_cif_file failed for " + path + ": " + msg);

    vc_sets frc = group_by_element(cs);
    vector_Affine3d syms = get_all_symmetries(cs.unit_cell.cell(), frc, 1E-2);

    // Compare inferred symmetries to the ones declared by the CIF's space group.
    const gemmi::SpaceGroup *sg = declared_spacegroup(path);
    BOOST_REQUIRE_MESSAGE(sg != nullptr,
                          "No space group resolved for " + path);

    int num_declared = 0;
    for (const auto &op : sg->operations()) {
      Affine3d et = gemmi_op_to_affine(op);
      bool exists = false;
      for (const auto &s : syms) {
        if ((s.linear() - et.linear()).norm() < 1E-4 &&
            min_frac(s.translation() - et.translation()).norm() < 1E-4) {
          exists = true;
          break;
        }
      }
      BOOST_CHECK_MESSAGE(exists, "Declared operation missing from inferred set in " + path);
      num_declared++;
    }
    BOOST_CHECK_EQUAL(int(syms.size()), num_declared);
  }
}

BOOST_AUTO_TEST_SUITE_END()
