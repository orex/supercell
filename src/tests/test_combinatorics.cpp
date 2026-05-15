#define BOOST_TEST_MODULE test_combinatorics

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include "science/combinatorics.h"

#ifndef COMBINATORICS_DATA_DIR
#define COMBINATORICS_DATA_DIR "."
#endif

namespace utf = boost::unit_test;

template <typename T>
void multiplication_check(int v1, int v2, int v3, int v4, bool overflow)
{
  std::vector<int> vm{v1, v2, v3, v4};

  T mult_result;
  bool multiplication_ok = safe_multiplication(vm, mult_result);

  std::string msg = (boost::format("v={%1%, %2%, %3%, %4%}") % v1 % v2 % v3 % v4).str();

  if (!overflow) {
    int64_t result_ref = static_cast<int64_t>(v1) * static_cast<int64_t>(v2) *
                         static_cast<int64_t>(v3) * static_cast<int64_t>(v4);
    BOOST_CHECK_MESSAGE(multiplication_ok, msg);
    BOOST_CHECK_EQUAL(mult_result, result_ref);
    BOOST_CHECK_MESSAGE(mult_result == result_ref, msg);
  } else {
    BOOST_CHECK_MESSAGE(!multiplication_ok, "Should be overflow: " + msg);
  }
}

static int64_t create_combinations(int v0, ...)
{
  va_list ap;
  va_start(ap, v0);
  std::vector<int> vc{v0};
  while (true) {
    int n = va_arg(ap, int);
    if (n < 0)
      break;
    vc.push_back(n);
  }
  va_end(ap);

  int64_t result;
  bool good = num_combinations(vc, result);
  if (good)
    BOOST_CHECK_GE(result, 0);
  else
    result = -1;
  return result;
}

BOOST_AUTO_TEST_SUITE(CombinatoricsTest)

BOOST_AUTO_TEST_CASE(Test_safe_mult_simple)
{
  multiplication_check<int>(1, 1, 1, 0, false);

  multiplication_check<int>(1, 1, 1, 1, false);
  multiplication_check<unsigned int>(1, 1, 1, 1, false);

  multiplication_check<int>(1, -1, 1, 1, false);
  multiplication_check<unsigned int>(1, -1, 1, 1, true);
  multiplication_check<unsigned int>(1, -1, -1, 1, false);

  multiplication_check<int64_t>(std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), 1, 1, false);
  multiplication_check<int64_t>(std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), 2, 1, false);
  multiplication_check<int64_t>(std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), 2, 2, true);

  multiplication_check<int>(32767, 32767, 2, 2, true);
  multiplication_check<unsigned int>(32767, 32767, 2, 2, false);
}

BOOST_AUTO_TEST_CASE(Test_safe_mult_random)
{
  std::mt19937 rnd(std::random_device{}());

  for (int i = 0; i < 200000; i++) {
    int maxn = std::uniform_int_distribution<int>(0, 10)(rnd);
    std::vector<int> vm;
    vm.reserve(maxn);
    bool all_positive = true;
    for (int j = 0; j < maxn; j++) {
      int a = std::uniform_int_distribution<int>(-32767, 32767)(rnd);
      all_positive = all_positive && (a > 0);
      vm.push_back(a);
    }

    int64_t mult_64;
    bool good64 = safe_multiplication(vm, mult_64);

    int mult_32;
    bool good32 = safe_multiplication(vm, mult_32);

    if (!good64 || good32)
      BOOST_CHECK_EQUAL(good32, good64);

    if (good64) {
      int64_t ref = 1;
      for (int j = 0; j < maxn; j++)
        ref *= vm[j];
      if (all_positive)
        BOOST_CHECK_GE(mult_64, 0);
      BOOST_CHECK_EQUAL(mult_64, ref);
      BOOST_CHECK_EQUAL(std::abs(mult_64) <= std::numeric_limits<int>::max(), good32);
      if (good32 && good64)
        BOOST_CHECK_EQUAL(int64_t(mult_32), mult_64);
    }
  }
}

BOOST_AUTO_TEST_CASE(Test_combinatorics_simple)
{
  BOOST_CHECK_EQUAL(create_combinations(0, -1), 0);
  BOOST_CHECK_EQUAL(create_combinations(1, -1), 1);
  BOOST_CHECK_EQUAL(create_combinations(1, 0, -1), 1);
  BOOST_CHECK_EQUAL(create_combinations(1, 1, -1), 2);
  BOOST_CHECK_EQUAL(create_combinations(100, 1, -1), 101);
  BOOST_CHECK_EQUAL(create_combinations(1000, 2, -1), 1002 * 1001 / 2);
}

BOOST_AUTO_TEST_CASE(Test_combinatorics_ref)
{
  const std::string path = std::string(COMBINATORICS_DATA_DIR) + "/ctest.dat";
  std::ifstream dat(path);
  BOOST_REQUIRE_MESSAGE(dat.is_open(), "Cannot open " + path);

  std::string line;
  while (std::getline(dat, line).good()) {
    std::vector<std::string> tks;
    boost::split(tks, line, boost::is_any_of(", "), boost::algorithm::token_compress_on);
    if (tks.size() < 2)
      continue;
    int64_t ref = boost::lexical_cast<int64_t>(tks[0]);
    std::vector<int> vc;
    for (size_t i = 1; i < tks.size(); i++)
      vc.push_back(boost::lexical_cast<int>(tks[i]));

    int64_t result;
    bool good = num_combinations(vc, result);
    if (!good)
      result = -1;
    BOOST_CHECK_MESSAGE(result == ref, std::string("Data: ") + line);
  }
}

BOOST_AUTO_TEST_CASE(Test_next_permutation_it_rnd)
{
  std::mt19937 rnd;
  for (int i = 0; i < 12; i++) {
    std::vector<int> v(i, 0);
    for (size_t j = 0; j < v.size(); j++)
      v[j] = std::uniform_int_distribution<int>(-1, 4)(rnd);
    std::sort(v.begin(), v.end());
    while (true) {
      auto vn = v;
      auto vl = v;
      auto it = next_permutation_it(vn.begin(), vn.end());
      bool b = std::next_permutation(vl.begin(), vl.end());
      BOOST_CHECK_EQUAL(b, (it != vn.end()));
      BOOST_CHECK_EQUAL_COLLECTIONS(vn.begin(), vn.end(), vl.begin(), vl.end());
      if (it == vn.end())
        break;
      int ix = std::distance(vn.begin(), it);
      BOOST_CHECK_NE(v[ix], vn[ix]);
      BOOST_CHECK_EQUAL_COLLECTIONS(vn.begin(), vn.begin() + ix,
                                    vl.begin(), vl.begin() + ix);
      std::swap(v, vn);
    }
  }
}

BOOST_AUTO_TEST_CASE(Test_permutation_index_rnd)
{
  std::mt19937 rnd;
  for (int i = 0; i < 12; i++) {
    std::vector<int> v(i, 0);
    for (size_t j = 0; j < v.size(); j++)
      v[j] = std::uniform_int_distribution<int>(-1, 4)(rnd);
    std::sort(v.begin(), v.end());
    int64_t count = 0;
    do {
      BOOST_CHECK_EQUAL(permutation_index<int64_t>(v.begin(), v.end()), count);
      count++;
    } while (std::next_permutation(v.begin(), v.end()));
  }
}

#define CHECK_NEXT_K_REF(arr, k, ref_array, ref_result)                        \
  {                                                                            \
    auto a = (arr);                                                            \
    auto result = next_k_permutations(a.begin(), a.end(), (k));                \
    BOOST_CHECK_EQUAL((result), (ref_result));                                 \
    BOOST_CHECK_EQUAL_COLLECTIONS(a.begin(), a.end(), (ref_array).begin(),     \
                                  (ref_array).end());                          \
  }

BOOST_AUTO_TEST_CASE(Test_next_k_permutation_trivial)
{
  const std::array<int, 0> empty({});
  const std::array<int, 1> one_e({1});
  const std::array<int, 2> two_same({1, 1});
  const std::array<int, 2> two_different_first({0, 1});
  const std::array<int, 2> two_different_last({1, 0});

  for (int i = 0; i < 12; i++) {
    CHECK_NEXT_K_REF(empty, i, empty, i);
    CHECK_NEXT_K_REF(one_e, i, one_e, i);
    CHECK_NEXT_K_REF(two_same, i, two_same, i);
    if (i % 2 == 0) {
      CHECK_NEXT_K_REF(two_different_first, i, two_different_first, i / 2);
      CHECK_NEXT_K_REF(two_different_last, i, two_different_last, i / 2);
    } else {
      CHECK_NEXT_K_REF(two_different_first, i, two_different_last, i / 2);
      CHECK_NEXT_K_REF(two_different_last, i, two_different_first, i / 2 + 1);
    }
  }
}

BOOST_AUTO_TEST_CASE(Test_next_k_permutation_same)
{
  std::vector<int> data;
  for (int i = 0; i < 12; i++) {
    for (int64_t k : {0LL, 12LL, 456LL, 1237LL, 1276897LL, 1573003353979063832LL}) {
      std::vector<int> a = data;
      BOOST_CHECK_EQUAL(next_k_permutations(a.begin(), a.end(), k), k);
      BOOST_CHECK_EQUAL_COLLECTIONS(a.begin(), a.end(), data.begin(), data.end());
    }
    data.emplace_back(2);
  }
}

BOOST_AUTO_TEST_CASE(Test_next_k_permutation_one_change)
{
  constexpr int64_t length = 12;
  for (int64_t i = 0; i < length; i++) {
    for (int64_t k : {0LL, 5LL, 12LL, 1717LL, 1237LL, 1276897LL, 1573003353979063832LL}) {
      std::vector<int> a(length, 1);
      a[i] = 0;
      int64_t l = next_k_permutations(a.begin(), a.end(), k);
      BOOST_CHECK_EQUAL(std::count(a.begin(), a.end(), 0), 1);
      BOOST_CHECK_EQUAL(std::count(a.begin(), a.end(), 1), length - 1);
      BOOST_CHECK_EQUAL(l, (k + i) / length);
      BOOST_CHECK_EQUAL(a[(k + i) % length], 0);
    }
  }
}

BOOST_AUTO_TEST_CASE(Test_next_k_permutation_extensive)
{
  constexpr int a_size = 6;
  constexpr int cmb_num = 180;
  constexpr int num_rep = 5;
  const std::array<int, a_size> start_p({0, 1, 2, 2, 3, 3});
  std::vector<std::array<int, a_size>> v(1, start_p);
  for (;;) {
    v.emplace_back(v.back());
    if (!std::next_permutation(v.back().begin(), v.back().end()))
      break;
  }
  v.pop_back();
  BOOST_CHECK_EQUAL(v.size(), size_t(cmb_num));

  for (int i = 0; i < cmb_num; i++) {
    for (int j = i; j < cmb_num * num_rep; j++) {
      int k = j - i;
      CHECK_NEXT_K_REF(v[i % cmb_num], k, v[j % cmb_num],
                       (j / cmb_num) - (i / cmb_num));
    }
  }
}

BOOST_AUTO_TEST_CASE(Test_next_k_permutation_benchmark, *utf::disabled())
{
  constexpr int64_t k_step = 50000;
  std::vector<int> start_p(34, 1);
  std::fill(start_p.begin(), start_p.begin() + 17, 0);
  volatile int64_t count_num = 0;
  auto wp = start_p;

  auto t0 = std::chrono::steady_clock::now();
  do {
    count_num++;
  } while (std::next_permutation(wp.begin(), wp.end()));
  auto t1 = std::chrono::steady_clock::now();

  BOOST_CHECK_EQUAL_COLLECTIONS(wp.begin(), wp.end(), start_p.begin(), start_p.end());
  auto et_old = std::chrono::duration<double>(t1 - t0).count();

  t0 = std::chrono::steady_clock::now();
  for (int i = 0; i < count_num / k_step; i++)
    next_k_permutations(wp.begin(), wp.end(), k_step);
  next_k_permutations(wp.begin(), wp.end(), count_num % k_step);
  t1 = std::chrono::steady_clock::now();
  auto et_new = std::chrono::duration<double>(t1 - t0).count();

  BOOST_CHECK_EQUAL_COLLECTIONS(wp.begin(), wp.end(), start_p.begin(), start_p.end());
  std::cout << "Performance boost: " << (et_old / et_new) << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()
