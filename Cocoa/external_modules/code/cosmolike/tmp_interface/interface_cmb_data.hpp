#include <carma.h>
#include <armadillo>

#ifndef __COSMOLIKE_GENERIC_INTERFACE_CMB_DATA_HPP
#define __COSMOLIKE_GENERIC_INTERFACE_CMB_DATA_HPP

namespace cosmolike_interface
{


class IPCMB
{
public:
  static IPCMB& get_instance()
  {
    static IPCMB instance;
    return instance;
  }
  ~IPCMB() = default;

  bool is_cov_kkkk_term_from_sim() const;
  bool is_bandpower() const;
  bool is_binmat_set() const;
  bool is_offset_set() const;

  int get_nbandpower() const;
  std::array<double,2> get_lrange_bandpower() const;
  int get_ncls() const;
  std::array<double,2> get_lrange_kappacmb() const;

  void set_binning_mat(std::string binned_matrix_filename);
  void set_theory_offset(std::string theory_offset_filename);

  double get_alpha_hartlap_cov_term_kkkk() const;
  double get_binned_matrix_with_correction(const int ci, const int cj) const;
  double get_theory_offset(const int ci) const;

  arma::Mat<double> get_binned_matrix_with_correction() const;
  arma::Mat<double> get_theory_offset() const;

private: 
  bool is_cov_kkkk_term_from_sim_ = true;
  bool is_bandpower_  = true;  
  bool is_binmat_set_ = false;
  bool is_offset_set_ = false;
  int nbandpower_ = -1;
  int ncl_ = -1;
  double alpha_hartlap_cov_kkkk_;
  std::array<double,2> lrange_kappacmb_;
  std::array<double,2> lrange_bandpower_;
  std::string binmat_filename_;               // LSS x CMB
  std::string offset_filename_;               // LSS x CMB
  arma::Col<double> theory_offset_;           // Eqn35 Planck-2018-VIII
  arma::Mat<double> binmat_with_correction_;  // Eqn35 Planck-2018-VIII
  
  IPCMB() = default;
  IPCMB(IP const&) = delete;

  friend void init_cmb_bandpower(const int,const int,const double);
  friend void init_binning_cmb_bandpower(const int, const int, const int)
};


} // namespace cosmolike_interface