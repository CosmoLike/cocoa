#include <string>
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
#include <tuple>

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

// boost library
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/replace.hpp>

// COSMOLIKE
#include "cosmolike/basics.h"
#include "cosmolike/baryons.h"
#include "cosmolike/cosmo3D.h"
#include "cosmolike/interface_aux.hpp"
#include "cosmolike/interface_baryons.hpp"
#include "cosmolike/interface_cosmology.hpp"
#include "cosmolike/interface_dv_3x2pt.hpp"
#include "cosmolike/interface_masked_data.hpp"
#include "cosmolike/redshift_spline.h"
#include "cosmolike/structs.h"

namespace cosmolike_interface
{

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void init_data_vector_size_3x2pt_real_space()
{
  spdlog::debug("{}: Begins", "init_data_vector_size_3x2pt_real_space");

  if (tomo.shear_Npowerspectra == 0)
  {
    spdlog::critical(
      "{}: {} not set prior to this function call",
      "init_data_vector_size_3x2pt_real_space", "tomo.shear_Npowerspectra");
    exit(1);
  }
  if (tomo.ggl_Npowerspectra == 0)
  {
    spdlog::critical(
      "{}: {} not set prior to this function call",
      "init_data_vector_size_3x2pt_real_space", "tomo.ggl_Npowerspectra");
    exit(1);
  }
  if (tomo.clustering_Npowerspectra == 0)
  {
    spdlog::critical(
      "{}: {} not set prior to this function call",
      "init_data_vector_size_3x2pt_real_space", "tomo.clustering_Npowerspectra");
    exit(1);
  }
  if (Ntable.Ntheta == 0) 
  {
    spdlog::critical(
      "{}: {} not set prior to this function call",
      "init_data_vector_size_3x2pt_real_space", "Ntable.Ntheta");
    exit(1);
  }

  like.Ndata =   Ntable.Ntheta*2*tomo.shear_Npowerspectra
               + Ntable.Ntheta*tomo.ggl_Npowerspectra
               + Ntable.Ntheta*tomo.clustering_Npowerspectra;

  spdlog::debug("{}: {} = {} selected.", 
    "init_data_vector_size_3x2pt_real_space", "Ndata", like.Ndata);
  
  spdlog::debug("{}: Ends", "init_data_vector_size_3x2pt_real_space");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void init_data_vector_size_3x2pt_fourier_space()
{
  spdlog::debug("{}: Begins", "init_data_vector_size_3x2pt_fourier_space");

  if (tomo.shear_Npowerspectra == 0)
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "init_data_vector_size_3x2pt_fourier_space", "tomo.shear_Npowerspectra");
    exit(1);
  }
  if (tomo.ggl_Npowerspectra == 0)
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "init_data_vector_size_3x2pt_fourier_space", "tomo.ggl_Npowerspectra");
    exit(1);
  }
  if (tomo.clustering_Npowerspectra == 0)
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "init_data_vector_size_3x2pt_fourier_space", 
      "tomo.clustering_Npowerspectra");
    exit(1);
  }
  if (like.Ncl == 0) 
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "init_data_vector_size_3x2pt_fourier_space", 
      "like.Ncl");
    exit(1);
  }

  like.Ndata =  like.Ncl*tomo.shear_Npowerspectra
                + like.Ncl*tomo.ggl_Npowerspectra
                + like.Ncl*tomo.clustering_Npowerspectra;
  
  spdlog::debug("{}: {} = {} selected.", 
    "init_data_vector_size_3x2pt_fourier_space", 
    "Ndata", like.Ndata);
  
  spdlog::debug("{}: Ends", "init_data_vector_size_3x2pt_fourier_space");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double compute_pm(const int zl, const int zs, const double theta)
{
  return PointMass::get_instance().get_pm(zl, zs, theta);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void compute_ss_real_masked(arma::Col<double>& data_vector, const int start)
{
  if (like.shear_shear == 1)
  {
    for (int nz=0; nz<tomo.shear_Npowerspectra; nz++)
    {
      const int z1 = Z1(nz);
      const int z2 = Z2(nz);

      for (int i=0; i<Ntable.Ntheta; i++)
      {
        int index = start + Ntable.Ntheta*nz + i;
        
        if (IP::get_instance().get_mask(index))
        {
          data_vector(index) = xi_pm_tomo(1, i, z1, z2, 1)*
            (1.0 + nuisance.shear_calibration_m[z1])*
            (1.0 + nuisance.shear_calibration_m[z2]);
        }
        
        index += Ntable.Ntheta*tomo.shear_Npowerspectra;
        
        if (IP::get_instance().get_mask(index))
        {
          data_vector(index) = xi_pm_tomo(-1, i, z1, z2, 1)*
            (1.0 + nuisance.shear_calibration_m[z1])*
            (1.0 + nuisance.shear_calibration_m[z2]);
        }
      }
    }
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void compute_gs_real_masked(arma::Col<double>& data_vector, const int start)
{
  if (like.shear_pos == 1)
  {
    for (int nz=0; nz<tomo.ggl_Npowerspectra; nz++)
    {
      const int zl = ZL(nz);
      const int zs = ZS(nz);
      
      for (int i=0; i<Ntable.Ntheta; i++)
      {
        const int index = start + Ntable.Ntheta*nz + i;
        if (IP::get_instance().get_mask(index))
        {
          const double logdt=(std::log(Ntable.vtmax)-std::log(Ntable.vtmin))/Ntable.Ntheta;
          const double thetamin = std::exp(log(Ntable.vtmin) + (i+0.) * logdt);
          const double thetamax = std::exp(log(Ntable.vtmin) + (i+1.) * logdt);
          const double theta = (2./3.) * (std::pow(thetamax,3) - std::pow(thetamin,3)) /
                                         (thetamax*thetamax    - thetamin*thetamin);
          data_vector(index) = (w_gammat_tomo(i,zl,zs,1)+
            compute_pm(zl,zs,theta))*(1.0+nuisance.shear_calibration_m[zs]);
        }
      }
    }
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void compute_gg_real_masked(arma::Col<double>& data_vector, const int start)
{
  if (like.pos_pos == 1)
  {
    for (int nz=0; nz<tomo.clustering_Npowerspectra; nz++)
    {
      for (int i=0; i<Ntable.Ntheta; i++)
      {
        const int index = start + Ntable.Ntheta*nz + i;
        if (IP::get_instance().get_mask(index))
          data_vector(index) = w_gg_tomo(i, nz, nz, like.adopt_limber_gg);
      }
    }
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

vector compute_data_vector_3x2pt_real_masked_any_order(arma::Col<int>::fixed<3> order)
{ // order = (1,2,3) => Cosmic Shear, ggl, gg
  spdlog::debug("{}: Begins", "compute_data_vector_3x2pt_real_masked_any_order");

  if (redshift.shear_nbin == 0)
  {
    spdlog::critical("{}: {} = 0 is invalid",
      "compute_data_vector_3x2pt_real_masked_any_order", "shear_Nbin");
    exit(1);
  }
  if (Ntable.Ntheta == 0)
  {
    spdlog::critical("{}: {} = 0 is invalid",
      "compute_data_vector_3x2pt_real_masked_any_order", "Ntheta");
    exit(1);
  }
  if (!IP::get_instance().is_mask_set())
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "compute_data_vector_3x2pt_real_masked_any_order", "mask");
    exit(1);
  }

  constexpr int sz = 3;

  auto indices = arma::conv_to<arma::Col<int>>::from(
      arma::stable_sort_index(order, "ascend")
    );
  arma::Col<int>::fixed<sz> sizes =
    {
      2*Ntable.Ntheta*tomo.shear_Npowerspectra,
      Ntable.Ntheta*tomo.ggl_Npowerspectra,
      Ntable.Ntheta*tomo.clustering_Npowerspectra
    };

  arma::Col<int>::fixed<sz> start = {0,0,0};

  for(int i=0; i<sz; i++)
    for(int j=0; j<indices(i); j++)
      start(i) += sizes(indices(j));
  
  arma::Col<double> data_vector(like.Ndata, arma::fill::zeros);
  
  compute_ss_real_masked(data_vector, start(0));

  compute_gs_real_masked(data_vector, start(1));

  compute_gg_real_masked(data_vector, start(2));

  spdlog::debug("{}: Ends", "compute_data_vector_3x2pt_real_masked_any_order");

  return data_vector;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

arma::Col<double> compute_data_vector_3x2pt_real_masked_any_order(
    arma::Col<double> Q,                // PC amplitudes
    arma::Col<int>::fixed<3> order
  )
{
  if (!BaryonScenario::get_instance().is_pcs_set())
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "compute_data_vector_3x2pt_real_masked_any_order", "baryon PCs");
    exit(1);
  }
  if (BaryonScenario::get_instance().get_pcs().row(0).n_elem < Q.n_elem)
  {
    spdlog::critical("{}: {} invalid PC amplitude vector or PC eigenvectors",
      "compute_data_vector_3x2pt_real_masked_any_order");
    exit(1);
  }

  arma::Col<double> dv = compute_data_vector_3x2pt_real_masked_any_order(order);
  
  if (BaryonScenario::get_instance().get_pcs().col(0).n_elem != dv.n_elem)
  {
    spdlog::critical("{}: {} invalid datavector or PC eigenvectors",
      "compute_data_vector_3x2pt_real_masked_any_order");
    exit(1);
  }

  #pragma omp parallel for
  for (int j=0; j<dv.n_elem; j++)
    for (int i=0; i<Q.n_elem; i++)
      if (IP::get_instance().get_mask(j))
        dv(j) += Q(i) * BaryonScenario::get_instance().get_pcs(j, i);
  return dv;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void compute_ss_fourier_masked(arma::Col<double>& data_vector, const int start)
{ // HERE WE ASSUME NLA EE TODO CHANGE NAMES
  if (like.shear_shear == 1)
  {
    for (int nz=0; nz<tomo.shear_Npowerspectra; nz++)
    {
      const int z1 = Z1(nz);
      const int z2 = Z2(nz);

      for (int i=0; i<like.Ncl; i++)
      {
        const int index = start + like.Ncl*nz + i;
        
        if (like.ell[i] < like.lmax_shear)
        {
          data_vector(index) = C_ss_tomo_limber(like.ell[i], z1, z2, 1)*
            (1.0 + nuisance.shear_calibration_m[z1])*
            (1.0 + nuisance.shear_calibration_m[z2]);
        }
        else
        {
          data_vector(index) = 0.0;
        }
      }
    }
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void compute_gs_fourier_masked(arma::Col<double>& data_vector, const int start)
{
  if (like.shear_pos == 1)
  {
    for (int nz=0; nz<tomo.shear_Npowerspectra; nz++)
    {
      const int zl = ZL(nz);
      const int zs = ZS(nz);

      for (int i=0; i<like.Ncl; i++)
      {
        const int index = start + like.Ncl*nz + i;
        
        if (test_kmax(like.ell[i], zl))
        {
          data_vector(index) = C_gs_tomo_limber(like.ell[i], zl, zs)*
              (1.0 + nuisance.shear_calibration_m[zs]);
        }
        else
        {
          data_vector(index) = 0.0;
        }
      }
    }
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void compute_gg_fourier_masked(arma::Col<double>& data_vector, const int start)
{
  if (like.pos_pos == 1)
  {
    for (int nz=0; nz<tomo.clustering_Npowerspectra; nz++)
    {
      for (int i=0; i<like.Ncl; i++)
      {
        const int index = start + Ntable.Ntheta*nz + i;

        if (test_kmax(like.ell[i], nz))
          data_vector(index) = C_gg_tomo_limber(like.ell[i], nz, nz);
        else
          data_vector(index) = 0.0;
      }
    }
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

vector compute_data_vector_3x2pt_fourier_masked_any_order(arma::Col<int>::fixed<3> order)
{ // order = (1,2,3) => Cosmic Shear, ggl, gg
  spdlog::debug("{}: Begins","compute_data_vector_3x2pt_fourier_masked_any_order");

  if (redshift.shear_nbin == 0)
  {
    spdlog::critical("{}: {} = 0 is invalid",
      "compute_data_vector_3x2pt_fourier_masked_any_order", "shear_Nbin");
    exit(1);
  }
  if (like.Ncl == 0)
  {
    spdlog::critical("{}: {} = 0 is invalid",
      "compute_data_vector_3x2pt_fourier_masked_any_order", "Ncl");
    exit(1);
  }

  // check if there are not equal elements
  constexpr int sz = 3;

  auto indices = arma::conv_to<arma::Col<int>>::from(
      arma::stable_sort_index(order, "ascend")
    );

  arma::Col<int>::fixed<sz> sizes =
    {
      like.Ncl*tomo.shear_Npowerspectra,
      like.Ncl*tomo.ggl_Npowerspectra,
      like.Ncl*tomo.clustering_Npowerspectra
    };

  arma::Col<int>::fixed<sz> start = {0,0,0};

  for(int i=0; i<sz; i++)
  {
    for(int j=0; j<indices(i); j++)
    {
      start(i) += sizes(indices(j));
    }
  }
  
  vector data_vector(like.Ndata, arma::fill::zeros);
  
  compute_ss_fourier_masked(data_vector, start(0));

  compute_gs_fourier_masked(data_vector, start(1));

  compute_gg_fourier_masked(data_vector, start(2));

  spdlog::debug("{}: Ends", "compute_data_vector_3x2pt_fourier_masked_any_order");

  return data_vector;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void PointMass::set_pm_vector(arma::Col<double> pm)
{
  this->pm_ = pm;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

arma::Col<double> PointMass::get_pm_vector() const
{
  return this->pm_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double PointMass::get_pm(
    const int zl, 
    const int zs, 
    const double theta
  ) const
{
  constexpr double G_over_c2 = 1.6e-23;
  
  const double a_lens = 1.0/(1.0 + zmean(zl));
  
  const double chi_lens = chi(a_lens);

  return 4*G_over_c2*this->pm_[zl]*1.e+13*
    g_tomo(a_lens, zs)/(theta*theta)/(chi_lens*a_lens);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
} // namespace cosmolike_interface
