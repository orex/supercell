/* 
 * File:   common_types.h
 * Author: kirill
 *
 * Created on July 29, 2013, 9:56 AM
 */

#ifndef COMMON_TYPES_H
#define	COMMON_TYPES_H

#include <string>
#include <map>
#include <vector>
#include <set>
#include <cmath>
#include <array>

#include <boost/align/aligned_allocator.hpp>
#include <boost/optional.hpp>

#include <Eigen/Core>
#include <random>
#include <deque>

enum class charge_balance {cb_no, cb_yes, cb_try};

class cb_names {
private:
  static constexpr std::size_t __psize = 3;
  static constexpr charge_balance __cbt[__psize] = {charge_balance::cb_no, charge_balance::cb_yes, charge_balance::cb_try};
  static constexpr const char * __cbn[__psize]   = {"no",                  "yes",                  "try"};
public:
  static std::string get_name(const charge_balance cb) {
    for(int i = 0; i < __psize; i++) {
      if( cb == __cbt[i] ) {
        return __cbn[i];
      }
    }
    return "";
  }
  static bool get_cb(const std::string &str, charge_balance &cb) {
    for(int i = 0; i < __psize; i++) {
      if( str == __cbn[i] ) {
        cb = __cbt[i];
        return true;
      }
    }
    return false;
  }
};

class c_man_atom_prop_item
{
public:
  boost::optional<double> charge;
  boost::optional<bool> fixed;
  boost::optional<int> population;
 
  inline void assign(const c_man_atom_prop_item &orig)
  {
    if( orig.charge.is_initialized() )
      charge = orig.charge;
    if( orig.fixed.is_initialized() )
      fixed = orig.fixed;
    if( orig.population.is_initialized() )
      population = orig.population;
  };
};        

class c_man_atom_prop
{
public:
  typedef std::map<std::string, c_man_atom_prop_item> data_type;
protected:
  data_type data_map;
public:
  const data_type &data()
  { return data_map; };
  c_man_atom_prop_item operator[](const std::string &lbl)
  { return data_map[lbl]; };
  bool exist(const std::string &lbl)
  { return data_map.count(lbl) != 0; };
  virtual void convert_properties(const std::set<std::string> &labels) = 0;
};

class c_struct_sel
{
private:
  typedef std::array<int, 26> sel_prop;
  sel_prop sp;
  bool empty;
  inline int str_count(char c) const {
    return sp[c - 'a'];
  };
public:
  c_struct_sel() : sp(), empty(true) {
    std::fill(sp.begin(), sp.end(), 0);
  }
  inline void set_sampling(char c, int count ) {
    if( c >= 'a' && c <= 'z' && count > 0 ) {
      sp[c - 'a'] = count;
      empty = false;
    }
  };
  inline bool save_all() const
  { return empty; };
  inline int str_first_count() const
  { return str_count('f'); };
  inline int str_last_count() const
  { return str_count('a'); };
  inline int str_low_count() const
  { return str_count('l'); };
  inline int str_high_count() const
  { return str_count('h'); };
  inline int str_random_count() const
  { return str_count('r'); };
  inline int str_weight_limit() const
  { return str_count('w'); };

  void assign_base(const c_struct_sel &orig) {
    empty = orig.empty;
    sp = orig.sp;
  };
};

typedef std::uint32_t base_prm_t;

typedef std::vector<base_prm_t, boost::alignment::aligned_allocator<base_prm_t, 32>> t_vec_comb;

struct struct_info_t {
  t_vec_comb cmb;
  double energy;
  //weight should be equal 0 for non-merge run.
  int weight;

  struct_info_t(): energy(0), weight(0) {};
};

struct struct_info_index_t : public struct_info_t {
  int64_t index;
  struct_info_index_t(): struct_info_t(), index(0) {};
  struct_info_index_t(const struct_info_t &base_info,
                      int64_t index_v): struct_info_t(base_info),
                                         index(index_v) {};
};

class unit_cell_t {
 private:
  Eigen::Matrix3d _cell;
  inline static double v_ang(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2) {
    return abs(std::atan2(v1.cross(v2).norm(), v1.dot(v2))) * 180.0 / M_PI;
  }
 public:
  unit_cell_t() = default;
  explicit unit_cell_t(const Eigen::Matrix3d &cell) : _cell(cell) {};
  void set(const Eigen::Matrix3d &cell) { _cell = cell; };
  Eigen::Vector3d lengths() const {
    return Eigen::Vector3d(_cell.col(0).norm(), _cell.col(1).norm(), _cell.col(2).norm());
  }
  Eigen::Vector3d angles_deg() const {
    return Eigen::Vector3d(v_ang(_cell.col(1), _cell.col(2)),
                           v_ang(_cell.col(0), _cell.col(2)),
                           v_ang(_cell.col(0), _cell.col(1)));

  }
  const Eigen::Matrix3d & cell() const {
    return _cell;
  }
};

class rnd_indexer_t {
public:
  enum class sampling_method_t {
    DISABLED,
    BY_INDEXES,
    ALL
  };
private:
  typedef std::mt19937_64 rnd_engine_t;
  static_assert(rnd_engine_t::min() == 0
                    && rnd_engine_t::max() == std::numeric_limits<std::uint64_t>::max(),
                "Random engine error");
private:
  sampling_method_t method;
  std::random_device::result_type initial_seed;
  int random_count;
  int64_t total_comb;
  int64_t total_samples;
  int64_t base_index;
  int64_t index_step;
  std::deque<int64_t> indexes;
  void create_indexes();
public:
  rnd_indexer_t(std::random_device::result_type seed): initial_seed(seed) {};
  void set_properties(int random_count, int64_t total_comb, int symm_op);
  inline sampling_method_t get_mode() const {
    return method;
  }

  inline void reserve_indexes(int64_t index) {
    while( method == sampling_method_t::BY_INDEXES
           && index >= base_index + index_step) {
      base_index += index_step;
      create_indexes();
    }
  }
  inline int64_t get_current_index() const {
    return indexes.empty() ? total_comb : indexes.front();
  }
  inline int64_t has_next_index() const {
    return !indexes.empty();
  }
  inline void pop_index()  {
    if(!indexes.empty())
      indexes.pop_front();
  }
  inline std::size_t reserve_size() {
    switch (method) {
    case sampling_method_t::DISABLED:
      return 0;
    case sampling_method_t::BY_INDEXES:
      return indexes.size();
    case sampling_method_t::ALL:
      return total_comb;
    default:
      assert(false); // "Wrong sampling method"
      return 0;
    };
  }
  void postprocess_rnd_container(std::vector<struct_info_index_t> &rnd_container);
};


struct atom_t {
  int el_num;
  std::string label;
  Eigen::Vector3d fract_pos;
  double occupancy;
  double charge = std::nan("");
};

struct cryst_structure_t {
  std::string block_name;
  std::string chem_name;
  unit_cell_t unit_cell;
  std::vector<atom_t> atoms;
};

#endif	/* COMMON_TYPES_H */

