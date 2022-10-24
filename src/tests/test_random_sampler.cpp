#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_random_sampler

#include <boost/test/unit_test.hpp>

#include <cstdarg>
#include <random>
#include <chrono>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/timer/timer.hpp>

#include "src/core/common_types.h"

BOOST_AUTO_TEST_SUITE(RandomSamplerTest)

BOOST_AUTO_TEST_CASE(Test_time) {
  double cm = 0;
  auto ctime = std::chrono::steady_clock::now();
  for(int i = 0; i < 1000000000; i++) {
    auto mt = std::chrono::steady_clock::now();
    double dt = std::chrono::duration<double>(mt - ctime).count();
    cm += dt;
  }
  BOOST_CHECK_GT(cm, 0.0);
}

BOOST_AUTO_TEST_CASE(Test_matrix) {
  std::array<std::random_device::result_type, 10> seeds = {
      606087199,  1124675372, 1505642234, 638846697, 20542946,
      1424417894, 1629022254, 1616284817, 141993466, 1873233613};

  //std::array<int64_t, 12> total_cmb = {
  //    1, 5, 11, 220, 792, 1001, 10010, 96525, 110110, 168168, 1387386, 8588580};

  std::array<int64_t, 8> total_cmb = {
      1, 5, 11, 220, 792, 1001, 10010, 96525};

  std::array<int, 10> random_count = {
      0, 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 10000000};

  std::array<int, 3> packet_size = {1, 1000, 100000};

  std::array<int, 5> symmetires = {1, 3, 12, 128, 1536};

  std::array<double, 5> fill_coeff = {0, 0.1, 0.5, 0.7, 1.0};

  int f_index = 0;
  for (const auto &tc : total_cmb) {
    for (const auto &rc : random_count) {
      for (const auto &sm : symmetires) {
        for (const auto &psz : packet_size) {
          for (const auto &fc : fill_coeff) {
            int64_t nump =
                std::min<int64_t>(fc * (tc - tc / sm) + tc / sm, tc);
            for (const auto &s : seeds) {
              rnd_indexer_t ri(s);
              ri.set_properties(rc, tc, sm);
              std::vector<struct_info> sc;
              sc.reserve(ri.reserve_size());
              if (ri.get_mode() == rnd_indexer_t::sampling_method_t::ALL) {
                for (int i = 0; i < nump; i++) {
                  struct_info c;
                  c.index = i;
                  sc.emplace_back(c);
                }
              }
              if (ri.get_mode() ==
                  rnd_indexer_t::sampling_method_t::BY_INDEXES) {
                for (int i = 0; i < nump;) {
                  ri.reserve_indexes(i + psz);
                  while (ri.has_next_index() &&
                         ri.get_current_index() < i + psz) {
                    struct_info c;
                    c.index = ri.get_current_index();
                    sc.emplace_back(c);
                    ri.pop_index();
                  }
                  i += psz;
                }
              }
              ri.postprocess_rnd_container(sc);
              BOOST_CHECK_EQUAL(sc.size(), std::min<int64_t>(rc, nump));
              auto cmp = [] (const struct_info &a, const struct_info &b) {
                return a.index < b.index;
              };
              BOOST_CHECK(std::is_sorted(sc.cbegin(), sc.cend(), cmp));
              BOOST_CHECK(std::unique(sc.begin(), sc.end(), cmp) == sc.cend());
            }
          }
        }
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
