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
#include "cosmolike/cosmo2D.h"
#include "cosmolike/redshift_spline.h"
#include "cosmolike/structs.h"

using Vector = arma::Col<double>;
using Matrix = arma::Mat<double>;
using Cube = arma::Cube<double>;

namespace cosmolike
{

static int has_b2_galaxies()
{
  int res = 0;
  for (int i=0; i<tomo.clustering_Nbin; i++) 
  {
    if (gbias.b2[i])
    {
      res = 1;
    }
  }
  return res;
}

static int nell()
{
  int res = 0;
  switch(like.IA_MODEL)
  {
    case IA_MODEL_TATT:
    { 
      res = Ntable.N_ell_TATT;
      break;
    }
    case IA_MODEL_NLA:
    {
      res = Ntable.N_ell;
      break;
    }
    default:
    {
      spdlog::critical("like.IA_MODEL = {} not supported", like.IA_MODEL);
      exit(1);
    }
  }
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

py::tuple xi_pm_tomo_cpp()
{
  if (tomo.shear_Nbin == 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} = 0 is invalid",
      "compute_data_vector_masked", 
      "shear_Nbin"
    );
    exit(1);
  }
  if (like.Ntheta == 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} = 0 is invalid",
      "compute_data_vector_masked", 
      "Ntheta"
    );
    exit(1);
  }
 
  Vector xp(like.Ntheta*tomo.shear_Npowerspectra, arma::fill::none);
  Vector xm(like.Ntheta*tomo.shear_Npowerspectra, arma::fill::none);

  for (int nz=0; nz<tomo.shear_Npowerspectra; nz++)
  {    
    for (int i=0; i<like.Ntheta; i++)
    {
      const int z1 = Z1(nz);
      const int z2 = Z2(nz);
      xp(like.Ntheta*nz+i) = xi_pm_tomo(1, i, z1, z2, 1);
      xm(like.Ntheta*nz+i) = xi_pm_tomo(-1, i, z1, z2, 1);
    }
  }

  return py::make_tuple(carma::col_to_arr(xp), carma::col_to_arr(xm));
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

Vector w_gammat_tomo_cpp()
{
  if (tomo.shear_Nbin == 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} = 0 is invalid",
      "compute_data_vector_masked", 
      "shear_Nbin"
    );
    exit(1);
  }

  if (tomo.clustering_Nbin == 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} = 0 is invalid",
      "compute_data_vector_masked", 
      "clustering_Nbin"
    );
    exit(1);
  }

  if (like.Ntheta == 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} = 0 is invalid",
      "compute_data_vector_masked", 
      "Ntheta"
    );
    exit(1);
  }
  
  Vector result(like.Ntheta*tomo.ggl_Npowerspectra, arma::fill::none);

  for (int nz=0; nz<tomo.ggl_Npowerspectra; nz++)
  {
    for (int i=0; i<like.Ntheta; i++)
    {
      result(like.Ntheta*nz+i) = w_gammat_tomo(i, ZL(nz), ZS(nz), 1);
    }
  }

  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

Vector w_gg_tomo_cpp()
{
  if (tomo.clustering_Nbin == 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} = 0 is invalid",
      "compute_data_vector_masked", 
      "clustering_Nbin"
    );
    exit(1);
  }

  if (like.Ntheta == 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} = 0 is invalid",
      "compute_data_vector_masked", 
      "Ntheta"
    );
    exit(1);
  }

  Vector result(like.Ntheta*tomo.clustering_Nbin, arma::fill::none);

  for (int nz=0; nz<tomo.clustering_Npowerspectra; nz++)
  {
    for (int i=0; i<like.Ntheta; i++)
    {
      result(like.Ntheta*nz + i) = w_gg_tomo(i, nz, nz, 0);
    }
  }

  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

Vector w_gk_tomo_cpp()
{
  if (tomo.clustering_Nbin == 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} = 0 is invalid",
      "compute_data_vector_masked", 
      "clustering_Nbin"
    );
    exit(1);
  }

  if (like.Ntheta == 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} = 0 is invalid",
      "compute_data_vector_masked", 
      "Ntheta"
    );
    exit(1);
  }

  Vector result(like.Ntheta*tomo.clustering_Nbin, arma::fill::none);

  for (int nz=0; nz<tomo.clustering_Nbin; nz++)
  {
    for (int i=0; i<like.Ntheta; i++)
    {
      result(like.Ntheta*nz + i) = w_gk_tomo(i, nz, 1);
    }
  }

  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

Vector w_ks_tomo_cpp()
{
  if (tomo.shear_Nbin == 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} = 0 is invalid",
      "compute_data_vector_masked", 
      "clustering_Nbin"
    );
    exit(1);
  }
  if (like.Ntheta == 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} = 0 is invalid",
      "compute_data_vector_masked", 
      "Ntheta"
    );
    exit(1);
  }

  Vector result(like.Ntheta*tomo.shear_Nbin, arma::fill::none);

  for (int nz=0; nz<tomo.clustering_Nbin; nz++)
  {
    for (int i=0; i<like.Ntheta; i++)
    {
      result(like.Ntheta*nz + i) = w_ks_tomo(i, nz, 1);
    }
  }

  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double C_ss_NLA_tomo_limber_cpp(const double l, const int ni, const int nj)
{
  return C_ss_tomo_limber_nointerp(l, ni, nj, 0, 0);
}

Matrix C_ss_NLA_tomo_limber_cpp(const Vector l)
{
  if (!(l.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    exit(1);
  }

  Matrix result(l.n_elem, tomo.shear_Npowerspectra);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<tomo.shear_Npowerspectra; nz++)
  { // init static variables
    const int ni = Z1(nz);
    const int nj = Z2(nz);
    double trash = C_ss_tomo_limber_nointerp(l(0), ni, nj, 0, 1);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for collapse(2)
  for (int nz=0; nz<tomo.shear_Npowerspectra; nz++)
  {
    for (int i=0; i<static_cast<int>(l.n_elem); i++)
    {
      const int ni = Z1(nz);
      const int nj = Z2(nz);
      result(i, nz) = C_ss_tomo_limber_nointerp(l(i), ni, nj, 0, 0);
    }
  }

  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

py::tuple C_ss_TATT_tomo_limber_cpp(const double l, const int ni, const int nj)
{
  return py::make_tuple(
    C_ss_TATT_EE_tomo_limber_nointerp(l, ni, nj, 0),
    C_ss_TATT_BB_tomo_limber_nointerp(l, ni, nj, 0) 
  );
}

py::tuple C_ss_TATT_tomo_limber_cpp(const Vector l)
{
  if (!(l.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    exit(1);
  }
  
  Matrix EE(l.n_elem, tomo.shear_Npowerspectra);
  Matrix BB(l.n_elem, tomo.shear_Npowerspectra);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<tomo.shear_Npowerspectra; nz++)
  { // init static variables
    const int ni = Z1(nz);
    const int nj = Z2(nz);
    double trash = C_ss_TATT_EE_tomo_limber_nointerp(l(0), ni, nj, 1);
    trash = C_ss_TATT_BB_tomo_limber_nointerp(l(0), ni, nj, 1);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for collapse(2)
  for (int nz=0; nz<tomo.shear_Npowerspectra; nz++)
  {
    for (int i=0; i<static_cast<int>(l.n_elem); i++)
    {
      const int ni = Z1(nz);
      const int nj = Z2(nz);
      EE(i, nz) = C_ss_TATT_EE_tomo_limber_nointerp(l(i), ni, nj, 0);
      BB(i, nz) = C_ss_TATT_BB_tomo_limber_nointerp(l(i), ni, nj, 0);
    }
  }

  return py::make_tuple(carma::mat_to_arr(EE), carma::mat_to_arr(BB));
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

py::tuple C_ss_tomo_limber_cpp(const double l, const int ni, const int nj)
{
  switch(like.IA_MODEL)
  {
    case IA_MODEL_TATT:
    {
      return C_ss_TATT_tomo_limber_cpp(l, ni, nj);
    }
    case IA_MODEL_NLA:
    {
      return py::make_tuple(C_ss_NLA_tomo_limber_cpp(l, ni, nj), 0.0);
    }
    default:
    {
      spdlog::critical("like.IA_MODEL = {} not supported", like.IA_MODEL);
      exit(1);
    }
  }
}

py::tuple C_ss_tomo_limber_cpp(Vector l)
{
  switch(like.IA_MODEL)
  {
    case IA_MODEL_TATT:
    {
      return C_ss_TATT_tomo_limber_cpp(l);
    }
    case IA_MODEL_NLA:
    {
      py::make_tuple(
        carma::mat_to_arr(C_ss_NLA_tomo_limber_cpp(l)),
        carma::mat_to_arr(Matrix{l.n_elem, tomo.shear_Npowerspectra, 
          arma::fill::zeros})
      );
    }
    default:
    {
      spdlog::critical("like.IA_MODEL = {} not supported", like.IA_MODEL);
      exit(1);
    }
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double C_gs_tomo_limber_cpp(const double l, const int ni, const int nj)
{
  return C_gs_tomo_limber_nointerp(l, ni, nj, 0, 0);
}

Matrix C_gs_tomo_limber_cpp(const Vector l)
{
  if (!(l.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    exit(1);
  }

  Matrix result(l.n_elem, tomo.ggl_Npowerspectra);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<tomo.ggl_Npowerspectra; nz++)
  { // init static variables
    double tmp = C_gs_tomo_limber_nointerp(l(0), ZL(nz), ZS(nz), 0, 1);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for collapse(2)
  for (int nz=0; nz<tomo.ggl_Npowerspectra; nz++)
  {
    for (int i=0; i<static_cast<int>(l.n_elem); i++)
    {
      result(i, nz) = C_gs_tomo_limber_nointerp(l(i), ZL(nz), ZS(nz), 0, 0);
    }
  }

  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double C_gg_tomo_limber_cpp(const double l, const int nz)
{
  return C_gg_tomo_limber_nointerp(l, nz, nz, 0, 0);
}

Matrix C_gg_tomo_limber_cpp(const Vector l)
{
  if (!(l.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    exit(1);
  }

  Matrix result(l.n_elem, tomo.clustering_Nbin);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<tomo.clustering_Nbin; nz++)
  { // init static variables
    double tmp = C_gg_tomo_limber_nointerp(l(0), nz, nz, 0, 1);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for collapse(2)
  for (int nz=0; nz<tomo.clustering_Nbin; nz++)
  {
    for (int i=0; i<static_cast<int>(l.n_elem); i++)
    {
      result(i, nz) = C_gg_tomo_limber_nointerp(l(i), nz, nz, 0, 0);
    }
  }

  return result;
}

Matrix C_gg_tomo_cpp(const Vector l)
{
  if (!(l.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    exit(1);
  }

  Matrix result = C_gg_tomo_limber_cpp(l);

  for (int nz=0; nz<tomo.clustering_Nbin; nz++)
  {
    arma::uvec idxs = arma::find(l < limits.LMAX_NOLIMBER);

    if (idxs.n_elem > 0)
    {
      const int L = 1;
      const double tolerance = 0.0075;    
      const double dev = 10. * tolerance;

      Vector Cl(limits.LMAX_NOLIMBER+1);
      //C_cl_tomo(1, nz, nz, Cl.memptr(), dev, tolerance);
      
      for (int i=0; i<static_cast<int>(idxs.n_elem); i++)
      {
        result(idxs(i), nz) = Cl(static_cast<int>(l(idxs(i))+1e-13));
      }
    }
  }

  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double C_gk_tomo_limber_cpp(const double l, const int ni)
{
  return C_gk_tomo_limber_nointerp(l, ni, 0, 0);
}

Matrix C_gk_tomo_limber_cpp(const Vector l)
{
  if (!(l.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    exit(1);
  }

  Matrix result(l.n_elem, tomo.clustering_Nbin);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<tomo.clustering_Nbin; nz++)
  { // init static variables
    double trash = C_gk_tomo_limber_nointerp(l(0), nz, 0, 1);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for collapse(2)
  for (int nz=0; nz<tomo.clustering_Nbin; nz++)
  {
    for (int i=0; i<static_cast<int>(l.n_elem); i++)
    {
      result(i, nz) = C_gk_tomo_limber_nointerp(l(i), nz, 0, 0);
    }
  }

  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double C_ks_tomo_limber_cpp(const double l, const int ni)
{
  return C_ks_tomo_limber_nointerp(l, ni, 0, 0);
}

Matrix C_ks_tomo_limber_cpp(const Vector l)
{
  if (!(l.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    exit(1);
  }

  Matrix result(l.n_elem, tomo.shear_Nbin);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<tomo.shear_Nbin; nz++)
  { // init static variables
    double tmp = C_ks_tomo_limber_nointerp(l(0), nz, 0, 1);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for collapse(2)
  for (int nz=0; nz<tomo.shear_Nbin; nz++)
  {
    for (int i=0; i<static_cast<int>(l.n_elem); i++)
    {
      result(i, nz) = C_ks_tomo_limber_nointerp(l(i), nz, 0, 0);
    }
  }

  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double C_gy_tomo_limber_cpp(const double l, const int ni)
{
  return C_gy_tomo_limber_nointerp(l, ni, 0, 0);
}

Matrix C_gy_tomo_limber_cpp(const Vector l)
{
  if (!(l.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    exit(1);
  }

  Matrix result(l.n_elem, tomo.clustering_Nbin);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<tomo.clustering_Nbin; nz++)
  { // init static variables
    double tmp = C_gy_tomo_limber_nointerp(l(0), nz, 0, 1);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for collapse(2)
  for (int nz=0; nz<tomo.clustering_Nbin; nz++)
  {
    for (int i=0; i<static_cast<int>(l.n_elem); i++)
    {
      const int ni = nz;
      result(i, nz) = C_gy_tomo_limber_nointerp(l(i), nz, 0, 0);
    }
  }

  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double C_ys_tomo_limber_cpp(const double l, const int ni)
{
  return C_ys_tomo_limber_nointerp(l, ni, 0, 0);
}

Matrix C_ys_tomo_limber_cpp(const Vector l)
{
  if (!(l.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    exit(1);
  }

  Matrix result(l.n_elem, tomo.shear_Nbin);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<tomo.shear_Nbin; nz++)
  { // init static variables
    double tmp = C_ys_tomo_limber_nointerp(l(0), nz, 0, 1);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for collapse(2)
  for (int nz=0; nz<tomo.shear_Nbin; nz++)
  {
    for (int i=0; i<l.n_elem; i++)
    {
      result(i, nz) = C_ys_tomo_limber_nointerp(l(i), nz, 0, 0);
    }
  }

  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double C_kk_limber_cpp(const double l)
{
  return C_kk_limber_nointerp(l, 0, 0);
}

Vector C_kk_limber_cpp(const Vector l)
{
  if (!(l.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    exit(1);
  }

  Vector result(l.n_elem);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  { // init static variables
    double tmp = C_kk_limber_nointerp(l(0), 0, 1);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for
  for (int i=0; i<l.n_elem; i++)
  {
    result(i) = C_kk_limber_nointerp(l(i), 0, 0);
  }

  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double C_ky_limber_cpp(const double l)
{
  return C_ky_limber_nointerp(l, 0, 0);
}

Vector C_ky_limber_nointerp_cpp(const Vector l)
{
  if (!(l.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    exit(1);
  }

  Vector result(l.n_elem);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  { // init static variables
    double tmp = C_ky_limber_nointerp(l(0), 0, 1);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for
  for (int i=0; i<l.n_elem; i++)
  {
    result(i) = C_ky_limber_nointerp(l(i), 0, 0);
  }

  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double C_yy_limber_cpp(double l)
{
  return C_yy_limber_nointerp(l, 0, 0);
}

Vector C_yy_limber_nointerp_cpp(const Vector l)
{
  if (!(l.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    exit(1);
  }

  Vector result(l.n_elem);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  { // init static variables
    double tmp = C_yy_limber_nointerp(l(0), 0, 1);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for
  for (int i=0; i<l.n_elem; i++)
  {
    result(i) = C_yy_limber_nointerp(l(i), 0, 0);
  }

  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_ss_NLA_tomo_limber_cpp(
  const double a, const double l, const int ni, const int nj)
{
  double ar[4] = {(double) ni, (double) nj, l, (double) 0.0};
  return int_for_C_ss_tomo_limber(a, (void*) ar); 
}

Cube int_for_C_ss_NLA_tomo_limber_cpp(Vector a, Vector l)
{
  if (!(l.n_elem > 0 && a.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    spdlog::critical("\x1b[90m{}\x1b[0m: a array size = {}", a.n_elem);
    exit(1);
  }

  Cube result(a.n_elem, l.n_elem, tomo.shear_Npowerspectra);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<tomo.shear_Npowerspectra; nz++)
  { // init static variables
    double tmp = int_for_C_ss_NLA_tomo_limber_cpp(a(0), l(0), Z1(nz), Z2(nz));
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for collapse(3)
  for (int nz=0; nz<tomo.shear_Npowerspectra; nz++)
  {
    for (int i=0; i<l.n_elem; i++)
    {
      for (int j=0; j<a.n_elem; j++)
      {
        result(j, i, nz) = int_for_C_ss_NLA_tomo_limber_cpp(a(j), l(i), Z1(nz), Z2(nz));
      }
    }
  }

  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

py::tuple int_for_C_ss_TATT_tomo_limber_cpp(const double a, const double l, const int ni, const int nj)
{
  double ar[4] = {(double) ni, (double) nj, l, (double) 0.0};

  return py::make_tuple(
    int_for_C_ss_TATT_EE_tomo_limber(a, (void*) ar),
    int_for_C_ss_TATT_BB_tomo_limber(a, (void*) ar)
  );
}

py::tuple int_for_C_ss_TATT_tomo_limber_cpp(const Vector a, const Vector l)
{
  if (!(l.n_elem > 0 && a.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    spdlog::critical("\x1b[90m{}\x1b[0m: a array size = {}", a.n_elem);
    exit(1);
  }

  Cube EE(a.n_elem, l.n_elem, tomo.shear_Npowerspectra);
  Cube BB(a.n_elem, l.n_elem, tomo.shear_Npowerspectra);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<tomo.shear_Npowerspectra; nz++)
  { // init static variables
    double ar[4] = {(double) Z1(nz), (double) Z2(nz), l(0), (double) 0.0};
    double tmp = int_for_C_ss_TATT_EE_tomo_limber(a(0), (void*) ar);
    tmp = int_for_C_ss_TATT_BB_tomo_limber(a(0), (void*) ar);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for collapse(3)
  for (int nz=0; nz<tomo.shear_Npowerspectra; nz++)
  {
    for (int i=0; i<l.n_elem; i++)
    {
      for (int j=0; j<a.n_elem; j++)
      {
        double ar[4] = {(double) Z1(nz), (double) Z2(nz), l(i), (double) 0.0}; 
        EE(j, i, nz) = int_for_C_ss_TATT_EE_tomo_limber(a(j), (void*) ar);
        BB(j, i, nz) = int_for_C_ss_TATT_BB_tomo_limber(a(j), (void*) ar);
      }
    }
  }

  return py::make_tuple(carma::cube_to_arr(EE), carma::cube_to_arr(BB));
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

py::tuple int_for_C_ss_tomo_limber_cpp(const double a, const double l, const int ni, const int nj)
{
  switch(like.IA_MODEL)
  {
    case IA_MODEL_TATT:
    {
      return int_for_C_ss_TATT_tomo_limber_cpp(a, l, ni, nj);
    }
    case IA_MODEL_NLA:
    {
      return py::make_tuple(int_for_C_ss_NLA_tomo_limber_cpp(a, l, ni, nj), 0.0);
    }
    default:
    {
      spdlog::critical("like.IA_MODEL = {} not supported", like.IA_MODEL);
      exit(1);
    }
  }
}

py::tuple int_for_C_ss_tomo_limber_cpp(const Vector a, const Vector l)
{
  switch(like.IA_MODEL)
  {
    case IA_MODEL_TATT:
    {
      return int_for_C_ss_TATT_tomo_limber_cpp(a, l);
    }
    case IA_MODEL_NLA:
    {
      return py::make_tuple(
        carma::cube_to_arr(int_for_C_ss_NLA_tomo_limber_cpp(a, l)),
        carma::cube_to_arr(
          Cube{a.n_elem, l.n_elem, tomo.shear_Npowerspectra, arma::fill::zeros})
      );
    }
    default:
    {
      spdlog::critical("like.IA_MODEL = {} not supported", like.IA_MODEL);
      exit(1);
    }
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_gs_tomo_limber_cpp(const double a, const double l, const int nl, const int ns)
{
  double ar[4] = {(double) nl, (double) ns, l, (double) 0.0};

  double res;

  switch(like.IA_MODEL)
  {
    case IA_MODEL_TATT:
    {
      res = int_for_C_gs_TATT_tomo_limber(a, (void*) ar);
      break;
    }
    case IA_MODEL_NLA:
    {
      if (has_b2_galaxies())
      {
        res = int_for_C_gs_tomo_limber_withb2(a, (void*) ar);
      }
      else 
      {
        res = int_for_C_gs_tomo_limber(a, (void*) ar);
      }
      break;
    }
    default:
    {
      spdlog::critical("like.IA_MODEL = {} not supported", like.IA_MODEL);
      exit(1);
    }
  }

  return res;
}

Cube int_for_C_gs_tomo_limber_cpp(const Vector a, const Vector l)
{
  if (!(l.n_elem > 0 && a.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    spdlog::critical("\x1b[90m{}\x1b[0m: a array size = {}", a.n_elem);
    exit(1);
  }

  Cube result(a.n_elem, l.n_elem, tomo.ggl_Npowerspectra);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<tomo.ggl_Npowerspectra; nz++)
  { // init static variables
    double tmp = int_for_C_gs_tomo_limber_cpp(a(0), l(0), ZL(nz), ZS(nz));
  }
  #pragma GCC diagnostic pop

  switch(like.IA_MODEL)
  {
    case IA_MODEL_TATT:
    {
      #pragma omp parallel for collapse(3)
      for (int nz=0; nz<tomo.ggl_Npowerspectra; nz++)
      {
        for (int i=0; i<l.n_elem; i++)
        {
          for (int j=0; j<a.n_elem; j++)
          {
            double ar[4] = {(double) ZL(nz), (double) ZS(nz), l(i), (double) 0};
            result(j, i, nz) = int_for_C_gs_TATT_tomo_limber(a(j), (void*) ar);
          }
        }
      }
      break;
    }
    case IA_MODEL_NLA:
    {
      if (has_b2_galaxies())
      {
        #pragma omp parallel for collapse(3)
        for (int nz=0; nz<tomo.ggl_Npowerspectra; nz++)
        {
          for (int i=0; i<l.n_elem; i++)
          {
            for (int j=0; j<a.n_elem; j++)
            {
              double ar[4] = {(double) ZL(nz), (double) ZS(nz), l(i), (double) 0};
              result(j, i, nz) = int_for_C_gs_tomo_limber_withb2(a(j), (void*) ar);
            }
          }
        }
      }
      else 
      {
        #pragma omp parallel for collapse(3)
        for (int nz=0; nz<tomo.ggl_Npowerspectra; nz++)
        {
          for (int i=0; i<l.n_elem; i++)
          {
            for (int j=0; j<a.n_elem; j++)
            {
              double ar[4] = {(double) ZL(nz), (double) ZS(nz), l(i), (double) 0};
              result(j, i, nz) = int_for_C_gs_tomo_limber(a(j), (void*) ar);
            }
          }
        }
      }
      break;
    }
    default:
    {
      spdlog::critical("like.IA_MODEL = {} not supported", like.IA_MODEL);
      exit(1);
    }
  }

  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_gg_tomo_limber_cpp(const double a, const double l, 
    const int ni, const int nj
  )
{
  double ar[4] = {(double) ni, (double) nj, l, (double) 0};

  double res = 0.0;
  
  if (has_b2_galaxies())
  {
    res = int_for_C_gg_tomo_limber_withb2(a, (void*) ar);
  }
  else
  {
    res = int_for_C_gg_tomo_limber(a, (void*) ar);
  }

  return res;
}

Cube int_for_C_gg_tomo_limber_cpp(const Vector a, const Vector l)
{
  if (!(l.n_elem > 0 && a.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    spdlog::critical("\x1b[90m{}\x1b[0m: a array size = {}", a.n_elem);
    exit(1);
  }

  const int NSIZE = tomo.clustering_Nbin;
  Cube result(a.n_elem, l.n_elem, NSIZE);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<NSIZE; nz++)
  { // init static variables
    double tmp = int_for_C_gg_tomo_limber_cpp(a(0), l(0), nz, nz);
  }
  #pragma GCC diagnostic pop

  if (has_b2_galaxies())
  {
    #pragma omp parallel for collapse(3)
    for (int nz=0; nz<NSIZE; nz++)
    {
      for (int i=0; i<l.n_elem; i++)
      {
        for (int j=0; j<a.n_elem; j++)
        {
          const int ni = nz;
          const int nj = nz;
          double ar[4] = {(double) ni, (double) nj, l(i), (double) 0};
          result(j, i, nz) = int_for_C_gg_tomo_limber_withb2(a(j), (void*) ar);
        }
      }
    }
  }
  else
  {
    #pragma omp parallel for collapse(3)
    for (int nz=0; nz<NSIZE; nz++)
    {
      for (int i=0; i<l.n_elem; i++)
      {
        for (int j=0; j<a.n_elem; j++)
        {
          const int ni = nz;
          const int nj = nz;
          double ar[4] = {(double) ni, (double) nj, l(i), (double) 0};
          result(j, i, nz) = int_for_C_gg_tomo_limber(a(j), (void*) ar);
        }
      }
    }
  }

  return result;
}

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
  {
    res = int_for_C_gk_tomo_limber_withb2(a, (void*) ar);
  }
  else
  {
    res = int_for_C_gk_tomo_limber(a, (void*) ar);
  }

  return res;
}

Cube int_for_C_gk_tomo_limber_cpp(
    const Vector a, 
    const Vector l
  )
{
  if (!(l.n_elem > 0 && a.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    spdlog::critical("\x1b[90m{}\x1b[0m: a array size = {}", a.n_elem);
    exit(1);
  }

  Cube result(a.n_elem, l.n_elem, tomo.clustering_Nbin);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<tomo.clustering_Nbin; nz++)
  { // init static variables
    double tmp = int_for_C_gk_tomo_limber_cpp(a(0), l(0), nz);
  }
  #pragma GCC diagnostic pop

  if (has_b2_galaxies())
  {
    #pragma omp parallel for collapse(3)
    for (int nz=0; nz<tomo.clustering_Nbin; nz++)
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
    for (int nz=0; nz<tomo.clustering_Nbin; nz++)
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

Cube int_for_C_gy_tomo_limber_cpp(
    const Vector a, 
    const Vector l
  )
{
  if (!(l.n_elem > 0 && a.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    spdlog::critical("\x1b[90m{}\x1b[0m: a array size = {}", a.n_elem);
    exit(1);
  }

  Cube result(a.n_elem, l.n_elem, tomo.clustering_Nbin);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<tomo.clustering_Nbin; nz++)
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
    for (int nz=0; nz<tomo.clustering_Nbin; nz++)
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

  switch(like.IA_MODEL)
  {
    case IA_MODEL_NLA:
    {
      res = int_for_C_ks_tomo_limber(a, (void*) ar);
    }
    default:
    {
      spdlog::critical("like.IA_MODEL = {} not supported", like.IA_MODEL);
      exit(1);
    }
  }

  return res;
}

Cube int_for_C_ks_tomo_limber_cpp(
    const Vector a, 
    const Vector l
  )
{
  if (!(l.n_elem > 0 && a.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    spdlog::critical("\x1b[90m{}\x1b[0m: a array size = {}", a.n_elem);
    exit(1);
  }

  Cube result(a.n_elem, l.n_elem, tomo.shear_Nbin);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<tomo.shear_Nbin; nz++)
  { // init static variables
    double tmp = int_for_C_ks_tomo_limber_cpp(a(0), l(0), nz);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for collapse(3)
  for (int nz=0; nz<tomo.shear_Nbin; nz++)
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

  switch(like.IA_MODEL)
  {
    case IA_MODEL_NLA:
    {
      res = int_for_C_ys_tomo_limber(a, (void*) ar); 
    }
    default:
    {
      spdlog::critical("like.IA_MODEL = {} not supported", like.IA_MODEL);
      exit(1);
    }
  }

  return res;
}

Cube int_for_C_ys_tomo_limber_cpp(
    const Vector a, 
    const Vector l
  )
{
  if (!(l.n_elem > 0 && a.n_elem > 0))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: l array size = {}", l.n_elem);
    spdlog::critical("\x1b[90m{}\x1b[0m: a array size = {}", a.n_elem);
    exit(1);
  }

  Cube result(a.n_elem, l.n_elem, tomo.shear_Nbin);

  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  for (int nz=0; nz<tomo.shear_Nbin; nz++)
  { // init static variables
    double tmp = int_for_C_ys_tomo_limber_cpp(a(0), l(0), nz);
  }
  #pragma GCC diagnostic pop

  #pragma omp parallel for collapse(3)
  for (int nz=0; nz<tomo.shear_Nbin; nz++)
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

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_kk_limber_cpp(const double a, const double l)
{
  double ar[2] = {l, (double) 0.0}; 
  return int_for_C_kk_limber(a, (void*) ar);
}

Matrix int_for_C_kk_limber_cpp(const Vector a, const Vector l)
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

  Matrix result(a.n_elem, l.n_elem);

  #pragma omp parallel for collapse(2)
  for (int i=0; i<l.n_elem; i++)
  {
    for (int j=0; j<a.n_elem; j++)
    {
      result(j, i) = int_for_C_kk_limber_cpp(a(j), l(i));
    }
  }

  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_ky_limber_cpp(const double a, const double l)
{
  double ar[2] = {l, (double) 0.0}; 
  return int_for_C_ky_limber(a, (void*) ar);
}

Matrix int_for_C_ky_limber_cpp(const Vector a, const Vector l)
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

  Matrix result(a.n_elem, l.n_elem);

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

Matrix int_for_C_yy_limber_cpp(const Vector a, const Vector l)
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

  Matrix result(a.n_elem, l.n_elem);

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

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

} // end namespace cosmolike

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------