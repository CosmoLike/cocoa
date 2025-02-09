#include <string>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <stdexcept>
#include <array>
#include <random>
#include <variant>
#include <cmath> 

// SPDLOG
#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/cfg/env.h>

// ARMADILLO LIB AND PYBIND WRAPPER (CARMA)
#include <carma.h>
#include <armadillo>

// Python Binding
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/pytypes.h>
namespace py = pybind11;

// cosmolike
#include "cosmolike/bias.h"
#include "cosmolike/IA.h"
#include "cosmolike/cosmo2D.h"
#include "cosmolike/redshift_spline.h"
#include "cosmolike/structs.h"

using vector = arma::Col<double>;
using matrix = arma::Mat<double>;
using cube = arma::Cube<double>;

namespace cosmolike_interface
{

static int has_b2_galaxies()
{
  int res = 0;
  for (int i=0; i<redshift.clustering_nbin; i++) 
    if (nuisance.gb[1][i])
      res = 1;
  return res;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

arma::Col<double> get_binning_real_space()
{  
  arma::Col<double> result(Ntable.Ntheta, arma::fill::none);
  
  const double logdt=(std::log(Ntable.vtmax)-std::log(Ntable.vtmin))/Ntable.Ntheta;
  for (int i = 0; i < Ntable.Ntheta; i++)
  {  
    const double thetamin = std::exp(log(Ntable.vtmin) + (i + 0.0) * logdt);
    const double thetamax = std::exp(log(Ntable.vtmin) + (i + 1.0) * logdt);
    const double theta = (2./ 3.) * (std::pow(thetamax,3) - std::pow(thetamin,3)) /
                                    (thetamax*thetamax    - thetamin*thetamin);
    result(i) = theta / 2.90888208665721580e-4; 
  }
  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

py::tuple xi_pm_tomo_cpp()
{ 
  arma::Mat<double> xp(Ntable.Ntheta,tomo.shear_Npowerspectra,arma::fill::none);
  arma::Mat<double> xm(Ntable.Ntheta,tomo.shear_Npowerspectra,arma::fill::none);

  for (int nz=0; nz<tomo.shear_Npowerspectra; nz++)
  {    
    for (int i=0; i<Ntable.Ntheta; i++)
    {
      const int z1 = Z1(nz);
      const int z2 = Z2(nz);
      xp(i,nz) = xi_pm_tomo(1, i, z1, z2, 1);
      xm(i,nz) = xi_pm_tomo(-1, i, z1, z2, 1);
    }
  }
  return py::make_tuple(carma::mat_to_arr(xp), carma::mat_to_arr(xm));
}

arma::Cube<double> w_gammat_tomo_cpp()
{  
  arma::Cube<double> result(Ntable.Ntheta,
                            redshift.clustering_nbin, 
                            redshift.shear_nbin,
                            arma::fill::zeros);
  
  for (int nz=0; nz<tomo.ggl_Npowerspectra; nz++)
    for (int i=0; i<Ntable.Ntheta; i++)
      result(i,ZL(nz),ZS(nz)) = w_gammat_tomo(i, ZL(nz), ZS(nz), 1);
  return result;
}

arma::Cube<double> w_gg_tomo_cpp()
{
  arma::Cube<double> result(Ntable.Ntheta,
                            redshift.clustering_nbin,
                            redshift.clustering_nbin,
                            arma::fill::zeros);

  for (int nz=0; nz<tomo.clustering_Npowerspectra; nz++)
    for (int i=0; i<Ntable.Ntheta; i++)
      result(i, nz, nz) = w_gg_tomo(i, nz, nz, 0);
  return result;
}

/*

arma::Col<double> w_gk_tomo_cpp()
{
  arma::Col<double> result(Ntable.Ntheta*redshift.clustering_nbin,arma::fill::none);
  for (int nz=0; nz<redshift.clustering_nbin; nz++)
    for (int i=0; i<Ntable.Ntheta; i++)
      result(Ntable.Ntheta*nz+i) = w_gk_tomo(i, nz, 1);
  return result;
}

arma::Col<double> w_ks_tomo_cpp()
{
  arma::Col<double> result(Ntable.Ntheta*redshift.shear_nbin,arma::fill::none);
  for (int nz=0; nz<redshift.clustering_nbin; nz++)
    for (int i=0; i<Ntable.Ntheta; i++)
      result(Ntable.Ntheta*nz+i) = w_ks_tomo(i, nz, 1);
  return result;
}
*/

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

py::tuple C_ss_tomo_limber_cpp(const double l, const int ni, const int nj)
{
  return py::make_tuple(
    C_ss_tomo_limber_nointerp(l, ni, nj, 1, 0),
    C_ss_tomo_limber_nointerp(l, ni, nj, 0, 0) 
  );
}

py::tuple C_ss_tomo_limber_cpp(const arma::Col<double> l)
{
  if (!(l.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    exit(1);
  }
  
  arma::Mat<double> EE(l.n_elem,tomo.shear_Npowerspectra,arma::fill::zeros);
  arma::Mat<double> BB(l.n_elem,tomo.shear_Npowerspectra,arma::fill::zeros);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  { // init static variables
    const int ni = Z1(0);
    const int nj = Z2(0);
    double trash = C_ss_tomo_limber_nointerp(l(0), ni, nj, 1, 1);
    trash = C_ss_tomo_limber_nointerp(l(0), ni, nj, 0, 1);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for collapse(2)
  for (int nz=0; nz<tomo.shear_Npowerspectra; nz++)
  {
    for (int i=0; i<static_cast<int>(l.n_elem); i++)
    {
      const int ni = Z1(nz);
      const int nj = Z2(nz);
      EE(i, nz) = C_ss_tomo_limber_nointerp(l(i), ni, nj, 1, 0);
      BB(i, nz) = C_ss_tomo_limber_nointerp(l(i), ni, nj, 0, 0);
    }
  }
  return py::make_tuple(carma::mat_to_arr(EE), carma::mat_to_arr(BB));
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

arma::Mat<double> gs_bins()
{
  arma::Mat<double> result(tomo.ggl_Npowerspectra, 2);
  for (int nz=0; nz<tomo.ggl_Npowerspectra; nz++)
  {
    result(nz,0) = ZL(nz);
    result(nz,1) = ZS(nz);
  }
  return result;
}

double C_gs_tomo_limber_cpp(const double l, const int ni, const int nj)
{
  return C_gs_tomo_limber_nointerp(l, ni, nj, 0);
}

arma::Cube<double> C_gs_tomo_limber_cpp(const arma::Col<double> l)
{
  if (!(l.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    exit(1);
  }

  arma::Cube<double> result(l.n_elem,
                            redshift.clustering_nbin, 
                            redshift.shear_nbin,
                            arma::fill::zeros);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<tomo.ggl_Npowerspectra; nz++)
  { // init static variables
    double tmp = C_gs_tomo_limber_nointerp(l(0),ZL(nz),ZS(nz),1);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for collapse(2)
  for (int nz=0; nz<tomo.ggl_Npowerspectra; nz++)
    for (int i=0; i<static_cast<int>(l.n_elem); i++)
      result(i,ZL(nz),ZS(nz))=C_gs_tomo_limber_nointerp(l(i),ZL(nz),ZS(nz),0);
  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double C_gg_tomo_limber_cpp(const double l, const int nz)
{
  return C_gg_tomo_limber_nointerp(l, nz, nz, 0);
}

arma::Cube<double> C_gg_tomo_limber_cpp(const arma::Col<double> l)
{
  if (!(l.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    exit(1);
  }

  arma::Cube<double> result(l.n_elem, 
                            redshift.clustering_nbin,
                            redshift.clustering_nbin,
                            arma::fill::zeros);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<redshift.clustering_nbin; nz++)
  { // init static variables
    double tmp = C_gg_tomo_limber_nointerp(l(0), 0, 0, 1);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for collapse(2)
  for (int nz=0; nz<redshift.clustering_nbin; nz++)
    for (int i=0; i<static_cast<int>(l.n_elem); i++)
      result(i, nz, nz) = C_gg_tomo_limber_nointerp(l(i), nz, nz, 0);
  return result;
}

arma::Cube<double> C_gg_tomo_cpp(const arma::Col<double> l)
{
  if (!(l.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    exit(1);
  }

  arma::Cube<double> result = C_gg_tomo_limber_cpp(l);

  for (int nz=0; nz<redshift.clustering_nbin; nz++)
  {
    arma::uvec idxs = arma::find(l<limits.LMAX_NOLIMBER);

    if (idxs.n_elem > 0)
    {
      const int L = 1;
      const double tolerance = 0.01;     // required fractional accuracy in C(l)
      const double dev = 10. * tolerance; // will be diff  exact vs Limber init to

      arma::Col<double> Cl(limits.LMAX_NOLIMBER+1);
      C_cl_tomo(L, nz, nz, Cl.memptr(), dev, tolerance);
      
      for (int i=0; i<static_cast<int>(idxs.n_elem); i++)
        result(idxs(i), nz, nz) = Cl(static_cast<int>(l(idxs(i))+1e-13));
    }
  }
  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double C_gk_tomo_limber_cpp(const double l, const int ni)
{
  return C_gk_tomo_limber_nointerp(l, ni, 0);
}

arma::Mat<double> C_gk_tomo_limber_cpp(const arma::Col<double> l)
{
  if (!(l.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    exit(1);
  }

  arma::Mat<double> result(l.n_elem, redshift.clustering_nbin);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<redshift.clustering_nbin; nz++)
  { // init static variables
    double trash = C_gk_tomo_limber_nointerp(l(0), nz, 1);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for collapse(2)
  for (int nz=0; nz<redshift.clustering_nbin; nz++)
    for (int i=0; i<static_cast<int>(l.n_elem); i++)
      result(i, nz) = C_gk_tomo_limber_nointerp(l(i), nz, 0);
  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double C_ks_tomo_limber_cpp(const double l, const int ni)
{
  return C_ks_tomo_limber_nointerp(l, ni, 0);
}

arma::Mat<double> C_ks_tomo_limber_cpp(const arma::Col<double> l)
{
  if (!(l.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    exit(1);
  }

  arma::Mat<double> result(l.n_elem, redshift.shear_nbin);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<redshift.shear_nbin; nz++)
  { // init static variables
    double tmp = C_ks_tomo_limber_nointerp(l(0), nz, 1);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for collapse(2)
  for (int nz=0; nz<redshift.shear_nbin; nz++)
    for (int i=0; i<static_cast<int>(l.n_elem); i++)
      result(i, nz) = C_ks_tomo_limber_nointerp(l(i), nz, 0);
  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*
double C_gy_tomo_limber_cpp(const double l, const int ni)
{
  return C_gy_tomo_limber_nointerp(l, ni, 0, 0);
}

arma::Mat<double> C_gy_tomo_limber_cpp(const arma::Col<double> l)
{
  if (!(l.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    exit(1);
  }

  arma::Mat<double> result(l.n_elem, redshift.clustering_nbin);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<redshift.clustering_nbin; nz++)
  { // init static variables
    double tmp = C_gy_tomo_limber_nointerp(l(0), nz, 0, 1);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for collapse(2)
  for (int nz=0; nz<redshift.clustering_nbin; nz++)
    for (int i=0; i<static_cast<int>(l.n_elem); i++)
      result(i, nz) = C_gy_tomo_limber_nointerp(l(i), nz, 0, 0);
  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double C_ys_tomo_limber_cpp(const double l, const int ni)
{
  return C_ys_tomo_limber_nointerp(l, ni, 0, 0);
}

arma::Mat<double> C_ys_tomo_limber_cpp(const arma::Col<double> l)
{
  if (!(l.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    exit(1);
  }

  arma::Mat<double> result(l.n_elem, redshift.shear_nbin);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<redshift.shear_nbin; nz++)
  { // init static variables
    double tmp = C_ys_tomo_limber_nointerp(l(0), nz, 0, 1);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for collapse(2)
  for (int nz=0; nz<redshift.shear_nbin; nz++)
    for (int i=0; i<l.n_elem; i++)
      result(i, nz) = C_ys_tomo_limber_nointerp(l(i), nz, 0, 0);
  return result;
}
*/

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double C_kk_limber_cpp(const double l)
{
  return C_kk_limber_nointerp(l, 0);
}

arma::Col<double> C_kk_limber_cpp(const arma::Col<double> l)
{
  if (!(l.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    exit(1);
  }

  arma::Col<double> result(l.n_elem);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  { // init static variables
    double tmp = C_kk_limber_nointerp(l(0), 1);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for
  for (int i=0; i<l.n_elem; i++)
    result(i) = C_kk_limber_nointerp(l(i), 0);
  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*
double C_ky_limber_cpp(const double l)
{
  return C_ky_limber_nointerp(l, 0, 0);
}

arma::Col<double> C_ky_limber_nointerp_cpp(const arma::Col<double> l)
{
  if (!(l.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    exit(1);
  }

  arma::Col<double> result(l.n_elem);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  { // init static variables
    double tmp = C_ky_limber_nointerp(l(0), 0, 1);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for
  for (int i=0; i<l.n_elem; i++)
    result(i) = C_ky_limber_nointerp(l(i), 0, 0);
  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double C_yy_limber_cpp(double l)
{
  return C_yy_limber_nointerp(l, 0, 0);
}

arma::Col<double> C_yy_limber_nointerp_cpp(const arma::Col<double> l)
{
  if (!(l.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    exit(1);
  }

  arma::Col<double> result(l.n_elem);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  { // init static variables
    double tmp = C_yy_limber_nointerp(l(0), 0, 1);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for
  for (int i=0; i<l.n_elem; i++)
    result(i) = C_yy_limber_nointerp(l(i), 0, 0);
  return result;
}
*/
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_ss_EE_tomo_limber_cpp(
    const double a, 
    const double l, 
    const int ni, 
    const int nj
  )
{
  double ar[4] = {(double) ni, (double) nj, l, 1};
  return int_for_C_ss_tomo_limber(a, (void*) ar); 
}

arma::Cube<double> int_for_C_ss_EE_tomo_limber_cpp(
    arma::Col<double> a, 
    arma::Col<double> l
  )
{
  if (!(l.n_elem > 0 && a.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    spdlog::critical("\x1b[90m{}\x1b[0m: a array size = {}", a.n_elem);
    exit(1);
  }

  arma::Cube<double> result(a.n_elem, l.n_elem, tomo.shear_Npowerspectra);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<tomo.shear_Npowerspectra; nz++)
  { // init static variables
    double tmp = int_for_C_ss_EE_tomo_limber_cpp(a(0),l(0),Z1(nz),Z2(nz));
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for collapse(3)
  for (int nz=0; nz<tomo.shear_Npowerspectra; nz++)
    for (int i=0; i<l.n_elem; i++)
      for (int j=0; j<a.n_elem; j++)
        result(j,i,nz) = int_for_C_ss_EE_tomo_limber_cpp(a(j),l(i),Z1(nz),Z2(nz));
  return result;
}

double int_for_C_ss_BB_tomo_limber_cpp(
    const double a, 
    const double l, 
    const int ni, 
    const int nj
  )
{
  double ar[4] = {(double) ni, (double) nj, l, 0};
  return int_for_C_ss_tomo_limber(a, (void*) ar); 
}

arma::Cube<double> int_for_C_ss_BB_tomo_limber_cpp(
    arma::Col<double> a, 
    arma::Col<double> l
  )
{
  if (!(l.n_elem > 0 && a.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    spdlog::critical("\x1b[90m{}\x1b[0m: a array size = {}", a.n_elem);
    exit(1);
  }

  arma::Cube<double> result(a.n_elem, l.n_elem, tomo.shear_Npowerspectra);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<tomo.shear_Npowerspectra; nz++)
  { // init static variables
    double tmp = int_for_C_ss_BB_tomo_limber_cpp(a(0), l(0), Z1(nz), Z2(nz));
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for collapse(3)
  for (int nz=0; nz<tomo.shear_Npowerspectra; nz++)
    for (int i=0; i<l.n_elem; i++)
      for (int j=0; j<a.n_elem; j++)
        result(j,i,nz) = int_for_C_ss_BB_tomo_limber_cpp(a(j),l(i),Z1(nz),Z2(nz));
  return result;
}

py::tuple int_for_C_ss_tomo_limber_cpp(
    const double a, 
    const double l, 
    const int ni, 
    const int nj
  )
{
  return py::make_tuple(
    int_for_C_ss_EE_tomo_limber_cpp(a, l, ni, nj),
    int_for_C_ss_BB_tomo_limber_cpp(a, l, ni, nj)
  );
}

py::tuple int_for_C_ss_tomo_limber_cpp(
  const arma::Col<double> a, const arma::Col<double> l)
{
  if (!(l.n_elem > 0 && a.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    spdlog::critical("\x1b[90m{}\x1b[0m: a array size = {}", a.n_elem);
    exit(1);
  }

  arma::Cube<double> EE(a.n_elem, l.n_elem, tomo.shear_Npowerspectra);
  arma::Cube<double> BB(a.n_elem, l.n_elem, tomo.shear_Npowerspectra);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  { // init static variables
    double tmp = int_for_C_ss_EE_tomo_limber_cpp(a(0), l(0), Z1(0), Z2(0));
    tmp = int_for_C_ss_BB_tomo_limber_cpp(a(0), l(0), Z1(0), Z2(0));
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for collapse(3)
  for (int nz=0; nz<tomo.shear_Npowerspectra; nz++)
  {
    for (int i=0; i<l.n_elem; i++)
    {
      for (int j=0; j<a.n_elem; j++)
      {
        EE(j,i,nz) = int_for_C_ss_EE_tomo_limber_cpp(a(j),l(i),Z1(nz),Z2(nz));
        BB(j,i,nz) = int_for_C_ss_BB_tomo_limber_cpp(a(j),l(i),Z1(nz),Z2(nz));
      }
    }
  }
  return py::make_tuple(carma::cube_to_arr(EE), carma::cube_to_arr(BB));
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_gs_tomo_limber_cpp(
    const double a, 
    const double l, 
    const int nl, 
    const int ns
  )
{
  double ar[4] = {(double) nl, (double) ns, l, (double) has_b2_galaxies()};
  return int_for_C_gs_tomo_limber(a, (void*) ar);
}

arma::Cube<double> int_for_C_gs_tomo_limber_cpp(
    const arma::Col<double> a, 
    const arma::Col<double> l
  )
{
  if (!(l.n_elem > 0 && a.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    spdlog::critical("\x1b[90m{}\x1b[0m: a array size = {}", a.n_elem);
    exit(1);
  }

  arma::Cube<double> result(a.n_elem, l.n_elem, tomo.ggl_Npowerspectra);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<tomo.ggl_Npowerspectra; nz++)
  { // init static variables
    double tmp = int_for_C_gs_tomo_limber_cpp(a(0),l(0),ZL(nz),ZS(nz));
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for collapse(3)
  for (int nz=0; nz<tomo.ggl_Npowerspectra; nz++)
    for (int i=0; i<l.n_elem; i++)
      for (int j=0; j<a.n_elem; j++)
        result(j,i,nz) = int_for_C_gs_tomo_limber_cpp(a(j),l(i),ZL(nz),ZS(nz));
  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_gg_tomo_limber_cpp(
    const double a, 
    const double l, 
    const int ni, 
    const int nj
  )
{
  double ar[5] = {(double) ni, (double) nj, l, 0, (double) has_b2_galaxies()};
  return int_for_C_gg_tomo_limber(a, (void*) ar);
}

arma::Cube<double> int_for_C_gg_tomo_limber_cpp(
    const arma::Col<double> a, 
    const arma::Col<double> l
  )
{
  if (!(l.n_elem > 0 && a.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    spdlog::critical("\x1b[90m{}\x1b[0m: a array size = {}", a.n_elem);
    exit(1);
  }

  arma::Cube<double> result(a.n_elem, l.n_elem, redshift.clustering_nbin);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<redshift.clustering_nbin; nz++)
  { // init static variablesf
    double tmp = int_for_C_gg_tomo_limber_cpp(a(0), l(0), nz, nz);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for collapse(3)
  for (int nz=0; nz<redshift.clustering_nbin; nz++)
    for (int i=0; i<l.n_elem; i++)
      for (int j=0; j<a.n_elem; j++)
        result(j, i, nz) = int_for_C_gg_tomo_limber_cpp(a(j), l(i), nz, nz);
  return result;
}


/*
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_gk_tomo_limber_cpp(
    const double a, 
    const double l, 
    const int nz
  )
{
  double ar[3] = {(double) nz, l, (double) 0.0};

  double res = 0.0;

  if (has_b2_galaxies())
    res = int_for_C_gk_tomo_limber_withb2(a, (void*) ar);
  else
    res = int_for_C_gk_tomo_limber(a, (void*) ar);

  return res;
}

arma::Cube<double> int_for_C_gk_tomo_limber_cpp(
    const arma::Col<double> a, 
    const arma::Col<double> l
  )
{
  if (!(l.n_elem > 0 && a.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    spdlog::critical("\x1b[90m{}\x1b[0m: a array size = {}", a.n_elem);
    exit(1);
  }

  arma::Cube<double> result(a.n_elem, l.n_elem, redshift.clustering_nbin);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<redshift.clustering_nbin; nz++)
  { // init static variables
    double tmp = int_for_C_gk_tomo_limber_cpp(a(0), l(0), nz);
  }
  #pragma GCC diagnostic pop

  if (has_b2_galaxies())
  {
    #pragma omp parallel for collapse(3)
    for (int nz=0; nz<redshift.clustering_nbin; nz++)
    {
      for (int i=0; i<l.n_elem; i++)
      {
        for (int j=0; j<a.n_elem; j++)
        {
          double ar[3] = {(double) nz, l(i), (double) 0.0};
          result(j, i, nz) = int_for_C_gk_tomo_limber_withb2(a(j), (void*) ar);
        }
      }
    }
  }
  else
  {
    #pragma omp parallel for collapse(3)
    for (int nz=0; nz<redshift.clustering_nbin; nz++)
    {
      for (int i=0; i<l.n_elem; i++)
      {
        for (int j=0; j<a.n_elem; j++)
        {
          double ar[3] = {(double) nz, l(i), (double) 0.0};
          result(j, i, nz) = int_for_C_gk_tomo_limber(a(j), (void*) ar);
        }
      }
    }
  }

  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*
double int_for_C_gy_tomo_limber_cpp(
    const double a, 
    const double l,
    const int nz
  )
{
  double ar[3] = {(double) nz, l, (double) 0.0};
  
  double res;

  if (has_b2_galaxies())
  {
    spdlog::critical("b2 not supported in C_gy_nointerp");
    exit(1);
  }
  else
  {
    res = int_for_C_gy_tomo_limber(a, (void*) ar);
  }

  return res;
}

arma::Cube<double> int_for_C_gy_tomo_limber_cpp(
    const arma::Col<double> a, 
    const arma::Col<double> l
  )
{
  if (!(l.n_elem > 0 && a.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    spdlog::critical("\x1b[90m{}\x1b[0m: a array size = {}", a.n_elem);
    exit(1);
  }

  arma::Cube<double> result(a.n_elem, l.n_elem, redshift.clustering_nbin);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<redshift.clustering_nbin; nz++)
  { // init static variables
    double tmp = int_for_C_gy_tomo_limber_cpp(a(0), l(0), nz);
  }
  #pragma GCC diagnostic pop

  if (has_b2_galaxies())
  {
    spdlog::critical("b2 not supported in C_gy_nointerp");
    exit(1);
  }
  else
  {
    #pragma omp parallel for collapse(3)
    for (int nz=0; nz<redshift.clustering_nbin; nz++)
    {
      for (int i=0; i<l.n_elem; i++)
      {
        for (int j=0; j<a.n_elem; j++)
        {
          double ar[3] = {(double) nz, l(i), (double) 0.0};
          result(j, i, nz) = int_for_C_gy_tomo_limber(a(j), (void*) ar);
        }
      }
    }
  }

  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_ks_tomo_limber_cpp(
    const double a, 
    const double l, 
    const int nz
  )
{
  double ar[3] = {(double) nz, l, (double) 0.0};

  double res = 0.0;

  switch(nuisance.IA_MODEL)
  {
    case IA_MODEL_NLA:
    {
      res = int_for_C_ks_tomo_limber(a, (void*) ar);
    }
    default:
    {
      spdlog::critical("nuisance.IA_MODEL = {} not supported", nuisance.IA_MODEL);
      exit(1);
    }
  }

  return res;
}

arma::Cube<double> int_for_C_ks_tomo_limber_cpp(
    const arma::Col<double> a, 
    const arma::Col<double> l
  )
{
  if (!(l.n_elem > 0 && a.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    spdlog::critical("\x1b[90m{}\x1b[0m: a array size = {}", a.n_elem);
    exit(1);
  }

  arma::Cube<double> result(a.n_elem, l.n_elem, redshift.shear_nbin);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<redshift.shear_nbin; nz++)
  { // init static variables
    double tmp = int_for_C_ks_tomo_limber_cpp(a(0), l(0), nz);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for collapse(3)
  for (int nz=0; nz<redshift.shear_nbin; nz++)
  {
    for (int i=0; i<l.n_elem; i++)
    {
      for (int j=0; j<a.n_elem; j++)
      {
        result(j, i, nz) = int_for_C_ks_tomo_limber_cpp(a(j), l(i), nz);
      }
    }
  }

  return result;
}

/*
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_ys_tomo_limber_cpp(
    const double a, 
    const double l, 
    const int nz
  )
{
  double ar[3] = {(double) nz, l, (double) 0.0};

  double res = 0.0;

  switch(nuisance.IA_MODEL)
  {
    case IA_MODEL_NLA:
    {
      res = int_for_C_ys_tomo_limber(a, (void*) ar); 
    }
    default:
    {
      spdlog::critical("nuisance.IA_MODEL = {} not supported", nuisance.IA_MODEL);
      exit(1);
    }
  }

  return res;
}

arma::Cube<double> int_for_C_ys_tomo_limber_cpp(
    const arma::Col<double> a, 
    const arma::Col<double> l
  )
{
  if (!(l.n_elem > 0 && a.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    spdlog::critical("\x1b[90m{}\x1b[0m: a array size = {}", a.n_elem);
    exit(1);
  }

  arma::Cube<double> result(a.n_elem, l.n_elem, redshift.shear_nbin);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<redshift.shear_nbin; nz++)
  { // init static variables
    double tmp = int_for_C_ys_tomo_limber_cpp(a(0), l(0), nz);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for collapse(3)
  for (int nz=0; nz<redshift.shear_nbin; nz++)
  {
    for (int i=0; i<l.n_elem; i++)
    {
      for (int j=0; j<a.n_elem; j++)
      {
        result(j, i, nz) = int_for_C_ys_tomo_limber_cpp(a(j), l(i), nz);
      }
    }
  }

  return result;
}
*/

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_kk_limber_cpp(const double a, const double l)
{
  double ar[2] = {l, (double) 0.0}; 
  return int_for_C_kk_limber(a, (void*) ar);
}

arma::Mat<double> int_for_C_kk_limber_cpp(
    const arma::Col<double> a, 
    const arma::Col<double> l
  )
{
  if (!(l.n_elem > 0 && a.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    spdlog::critical("\x1b[90m{}\x1b[0m: a array size = {}", a.n_elem);
    exit(1);
  }

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  { // init static variables
    double tmp = int_for_C_kk_limber_cpp(a(0), l(0));
  }
  #pragma GCC diagnostic pop

  arma::Mat<double> result(a.n_elem, l.n_elem);

  #pragma omp parallel for collapse(2)
  for (int i=0; i<l.n_elem; i++)
    for (int j=0; j<a.n_elem; j++)
      result(j, i) = int_for_C_kk_limber_cpp(a(j), l(i));
  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*
double int_for_C_ky_limber_cpp(const double a, const double l)
{
  double ar[2] = {l, (double) 0.0}; 
  return int_for_C_ky_limber(a, (void*) ar);
}

arma::Mat<double> int_for_C_ky_limber_cpp(
    const arma::Col<double> a, 
    const arma::Col<double> l
  )
{
  if (!(l.n_elem > 0 && a.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    spdlog::critical("\x1b[90m{}\x1b[0m: a array size = {}", a.n_elem);
    exit(1);
  }

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  { // init static variables
    double tmp = int_for_C_ky_limber_cpp(a(0), l(0));
  }
  #pragma GCC diagnostic pop

  arma::Mat<double> result(a.n_elem, l.n_elem);

  #pragma omp parallel for collapse(2)
  for (int i=0; i<l.n_elem; i++)
  {
    for (int j=0; j<a.n_elem; j++)
    {
      result(j, i) = int_for_C_ky_limber_cpp(a(j), l(i));
    }
  }

  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_yy_limber_cpp(const double a, const double l)
{
  double ar[2] = {l, (double) 0.0}; 
  return int_for_C_yy_limber(a, (void*) ar);
}

arma::Mat<double> int_for_C_yy_limber_cpp(
    const arma::Col<double> a, 
    const arma::Col<double> l
  )
{
  if (!(l.n_elem > 0 && a.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    spdlog::critical("\x1b[90m{}\x1b[0m: a array size = {}", a.n_elem);
    exit(1);
  }

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  { // init static variables
    double tmp = int_for_C_yy_limber_cpp(a(0), l(0));
  }
  #pragma GCC diagnostic pop

  arma::Mat<double> result(a.n_elem, l.n_elem);

  #pragma omp parallel for collapse(2)
  for (int i=0; i<l.n_elem; i++)
  {
    for (int j=0; j<a.n_elem; j++)
    {
      result(j, i) = int_for_C_yy_limber_cpp(a(j), l(i));
    }
  }

  return result;
}

*/

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

} // end namespace cosmolike_interface

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------