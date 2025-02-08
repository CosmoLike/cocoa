#include <carma.h>
#include <armadillo>
#include <map>

// Python Binding
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>

#ifndef __COSMOLIKE_COSMO2D_WRAPPER_HPP
#define __COSMOLIKE_COSMO2D_WRAPPER_HPP

namespace cosmolike_interface
{

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

arma::Col<double> get_binning_real_space();

pybind11::tuple xi_pm_tomo_cpp();

arma::Col<double> w_gammat_tomo_cpp();

arma::Col<double> w_gg_tomo_cpp();

arma::Col<double> w_gk_tomo_cpp();

arma::Col<double> w_ks_tomo_cpp();

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

py::tuple C_ss_tomo_limber_cpp(
    const double l, 
    const int ni, 
    const int nj
  );

py::tuple C_ss_tomo_limber_cpp(
    const arma::Col<double> l
  );

py::tuple int_for_C_ss_tomo_limber_cpp(
    const double a, 
    const double l, 
    const int ni, 
    const int nj
  );

py::tuple int_for_C_ss_tomo_limber_cpp(
    const arma::Col<double> a, 
    const arma::Col<double> l
  );

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

arma::Mat<double> gs_bins();

double C_gs_tomo_limber_cpp(
    const double l, 
    const int ni, 
    const int nj
  );

arma::Mat<double> C_gs_tomo_limber_cpp(
    const arma::Col<double> l
  );

double int_for_C_gs_tomo_limber_cpp(
    const double a, 
    const double l, 
    const int nl, 
    const int ns
  );

arma::Cube<double> int_for_C_gs_tomo_limber_cpp(
    const arma::Col<double> a, 
    const arma::Col<double> l
  );



// ---------------------------------------------------------------------------

arma::Mat<double> C_gg_tomo_cpp(
    const arma::Col<double> l
  );

// ---------------------------------------------------------------------------

double C_gk_tomo_limber_cpp(
    const double l, 
    const int ni
  );

arma::Mat<double> C_gk_tomo_limber_cpp(
    const arma::Col<double> l
  );

double int_for_C_gk_tomo_limber_cpp(
    const double a, 
    const double l, 
    const int nz
  );

arma::Cube<double> int_for_C_gk_tomo_limber_cpp(
    const arma::Col<double> a, 
    const arma::Col<double> l
  );

// ---------------------------------------------------------------------------

double C_ks_tomo_limber_cpp(
    const double l, 
    const int ni
  );

arma::Mat<double> C_ks_tomo_limber_cpp(
    const arma::Col<double> l
  );

double int_for_C_ks_tomo_limber_cpp(
    const double a, 
    const double l, 
    const int nz
  );

arma::Cube<double> int_for_C_ks_tomo_limber_cpp(
    const arma::Col<double> a, 
    const arma::Col<double> l
  );

// ---------------------------------------------------------------------------

double C_kk_limber_cpp(
    const double l
  );

arma::Col<double> C_kk_limber_cpp(
    const arma::Col<double> l
  );

double int_for_C_kk_limber_cpp(
    const double a, 
    const double l
  );

arma::Mat<double> int_for_C_kk_limber_cpp(
    const arma::Col<double> a, 
    const arma::Col<double> l
  );

// ---------------------------------------------------------------------------

}  // namespace cosmolike_interface
#endif // HEADER GUARD