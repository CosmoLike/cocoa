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
#include "cosmolike/interface_cosmology.hpp"
#include "cosmolike/structs.h"

namespace cosmolike_interface
{

void set_cosmological_parameters(const double omega_matter, const double hubble)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "set_cosmological_parameters");

  // Cosmolike should not need parameters from inflation or dark energy.
  // Cobaya provides P(k,z), H(z), D(z), Chi(z)...
  // It may require H0 to set scales and \Omega_M to set the halo model

  int cache_update = 0;
  if (fdiff(cosmology.Omega_m, omega_matter) ||
      fdiff(cosmology.h0, hubble/100.0)) // assuming H0 in km/s/Mpc
  {
    cache_update = 1;
  }
  if (1 == cache_update)
  {
    cosmology.Omega_m = omega_matter;
    cosmology.Omega_v = 1.0-omega_matter;
    // Cosmolike only needs to know that there are massive neutrinos (>0)
    cosmology.Omega_nu = 0.1;
    cosmology.h0 = hubble/100.0; 
    cosmology.MGSigma = 0.0;
    cosmology.MGmu = 0.0;
    cosmology.random = cosmolike_interface::RandomNumber::get_instance().get();
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "set_cosmological_parameters");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void set_distances(arma::Col<double> io_z, arma::Col<double> io_chi)
{
  spdlog::debug("{}: Begins", "set_distances");

  bool debug_fail = false;
  if (io_z.n_elem != io_chi.n_elem)
    debug_fail = true;
  else
    if (io_z.n_elem == 0)
      debug_fail = true;
  
  if (debug_fail)
  {
    spdlog::critical("{}: incompatible input w/ z.size = {} and G.size = {}",
      "set_distances", io_z.n_elem, io_chi.n_elem);
    exit(1);
  }
  if(io_z.n_elem < 5)
  {
    spdlog::critical("{}: bad input w/ z.size = {} and chi.size = {}"
      "set_distances", io_z.n_elem, io_chi.n_elem);
    exit(1);
  }

  int cache_update = 0;
  if (cosmology.chi_nz != static_cast<int>(io_z.n_elem) || NULL == cosmology.chi)
    cache_update = 1;
  else
  {
    for (int i=0; i<cosmology.chi_nz; i++) 
    {
      if (fdiff(cosmology.chi[0][i], io_z(i)) ||
          fdiff(cosmology.chi[1][i], io_chi(i))) 
      {
        cache_update = 1; 
        break; 
      }    
    }
  }
  if (1 == cache_update)
  {
    cosmology.chi_nz = static_cast<int>(io_z.n_elem);

    if (cosmology.chi != NULL) free(cosmology.chi);
    cosmology.chi = (double**) malloc2d(2, cosmology.chi_nz);

    #pragma omp parallel for
    for (int i=0; i<cosmology.chi_nz; i++)
    {
      cosmology.chi[0][i] = io_z(i);
      cosmology.chi[1][i] = io_chi(i);
    }

    cosmology.random = RandomNumber::get_instance().get();
  }

  spdlog::debug("{}: Ends", "set_distances");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void set_growth(arma::Col<double> io_z, arma::Col<double> io_G)
{ // Growth: D = G * a
  spdlog::debug("{}: Begins", "set_growth");

  bool debug_fail = false;
  if (io_z.n_elem != io_G.n_elem)
    debug_fail = true;
  else
    if (io_z.n_elem == 0)
      debug_fail = true;
  
  if (debug_fail)
  {
    spdlog::critical("{}: incompatible input w/ z.size = {} and G.size = {}",
      "set_growth", io_z.n_elem, io_G.n_elem);
    exit(1);
  }
  if(io_z.n_elem < 5)
  {
    spdlog::critical("{}: bad input w/ z.size = {} and G.size = {}"
      "set_growth", io_z.n_elem, io_G.n_elem);
    exit(1);
  }

  int cache_update = 0;
  if (cosmology.G_nz != static_cast<int>(io_z.n_elem) || NULL == cosmology.G)
    cache_update = 1;
  else
  {
    for (int i=0; i<cosmology.G_nz; i++) 
    {
      if (fdiff(cosmology.G[0][i], io_z(i)) ||
          fdiff(cosmology.G[1][i], io_G(i))) 
      {
        cache_update = 1; 
        break;
      }    
    }
  }
  if (1 == cache_update)
  {
    cosmology.G_nz = static_cast<int>(io_z.n_elem);

    if (cosmology.G != NULL) free(cosmology.G);
    cosmology.G = (double**) malloc2d(2, cosmology.G_nz);

    #pragma omp parallel for
    for (int i=0; i<cosmology.G_nz; i++)
    {
      cosmology.G[0][i] = io_z(i);
      cosmology.G[1][i] = io_G(i);
    }

    cosmology.random = RandomNumber::get_instance().get();
  }

  spdlog::debug("{}: Ends", "set_growth");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void set_linear_power_spectrum(arma::Col<double> io_log10k, 
  arma::Col<double> io_z, arma::Col<double> io_lnP)
{
  spdlog::debug("{}: Begins", "set_linear_power_spectrum");

  bool debug_fail = false;
  if (io_z.n_elem*io_log10k.n_elem != io_lnP.n_elem)
    debug_fail = true;
  else
    if (io_z.n_elem == 0 || io_log10k.n_elem == 0)
      debug_fail = true;
  
  if (debug_fail)
  {
    spdlog::critical("{}: incompatible input w/ k.size = {}, z.size = {}, "
      "and lnP.size = {}", "set_linear_power_spectrum", 
      io_log10k.n_elem, io_z.n_elem, io_lnP.n_elem);
    exit(1);
  }
  if(io_z.n_elem < 5 || io_log10k.n_elem < 5)
  {
    spdlog::critical("{}: bad input w/ k.size = {}, z.size = {}, "
      "and lnP.size = {}", "set_linear_power_spectrum", 
      io_log10k.n_elem, io_z.n_elem, io_lnP.n_elem);
    exit(1);
  }

  int cache_update = 0;
  if (cosmology.lnPL_nk != static_cast<int>(io_log10k.n_elem) ||
      cosmology.lnPL_nz != static_cast<int>(io_z.n_elem) || 
      NULL == cosmology.lnPL)
    cache_update = 1;
  else
  {
    for (int i=0; i<cosmology.lnPL_nk; i++)
    {
      for (int j=0; j<cosmology.lnPL_nz; j++)
      {
        if (fdiff(cosmology.lnPL[i][j], io_lnP(i*cosmology.lnPL_nz+j))) 
        {
          cache_update = 1; 
          goto jump;
        }
      }
    }
    for (int i=0; i<cosmology.lnPL_nk; i++)
    {
      if (fdiff(cosmology.lnPL[i][cosmology.lnPL_nz], io_log10k(i))) 
      {
        cache_update = 1; 
        break;
      }
    }
    for (int j=0; j<cosmology.lnPL_nz; j++)
    {
      if (fdiff(cosmology.lnPL[cosmology.lnPL_nk][j], io_z(j))) 
      {
        cache_update = 1; 
        break;
      }
    }
  }

  jump:

  if (1 == cache_update)
  {
    cosmology.lnPL_nk = static_cast<int>(io_log10k.n_elem);
    cosmology.lnPL_nz = static_cast<int>(io_z.n_elem);

    if (NULL != cosmology.lnPL) free(cosmology.lnPL);
    cosmology.lnPL = (double**) malloc2d(cosmology.lnPL_nk+1,cosmology.lnPL_nz+1);

    #pragma omp parallel for
    for (int i=0; i<cosmology.lnPL_nk; i++)
      cosmology.lnPL[i][cosmology.lnPL_nz] = io_log10k(i);
    #pragma omp parallel for
    for (int j=0; j<cosmology.lnPL_nz; j++)
      cosmology.lnPL[cosmology.lnPL_nk][j] = io_z(j);
    #pragma omp parallel for collapse(2)
    for (int i=0; i<cosmology.lnPL_nk; i++)
      for (int j=0; j<cosmology.lnPL_nz; j++)
        cosmology.lnPL[i][j] = io_lnP(i*cosmology.lnPL_nz+j);

    cosmology.random = RandomNumber::get_instance().get();
  }

  spdlog::debug("{}: Ends", "set_linear_power_spectrum");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void set_non_linear_power_spectrum(arma::Col<double> io_log10k, 
  arma::Col<double> io_z, arma::Col<double> io_lnP)
{
  spdlog::debug("{}: Begins", "set_non_linear_power_spectrum");

  bool debug_fail = false;
  
  if (io_z.n_elem*io_log10k.n_elem != io_lnP.n_elem)
    debug_fail = true;
  else
    if (io_z.n_elem == 0)
      debug_fail = true;
  
  if (debug_fail)
  {
    spdlog::critical("{}: incompatible input w/ k.size = {}, z.size = {}, "
      "and lnP.size = {}", "set_non_linear_power_spectrum", 
      io_log10k.n_elem, io_z.n_elem, io_lnP.n_elem);
    exit(1);
  }
  if (io_z.n_elem < 5 || io_log10k.n_elem < 5)
  {
    spdlog::critical("{}: bad input w/ k.size = {}, z.size = {}, "
      "and lnP.size = {}", "set_non_linear_power_spectrum", 
      io_log10k.n_elem, io_z.n_elem, io_lnP.n_elem);
    exit(1);
  }

  int cache_update = 0;
  if (cosmology.lnP_nk != static_cast<int>(io_log10k.n_elem) ||
      cosmology.lnP_nz != static_cast<int>(io_z.n_elem) || 
      NULL == cosmology.lnP)
    cache_update = 1;
  else
  {
    for (int i=0; i<cosmology.lnP_nk; i++)
    {
      for (int j=0; j<cosmology.lnP_nz; j++)
      {
        if (fdiff(cosmology.lnP[i][j], io_lnP(i*cosmology.lnP_nz+j))) 
        {
          cache_update = 1; 
          goto jump;
        }
      }
    }
    for (int i=0; i<cosmology.lnP_nk; i++)
    {
      if (fdiff(cosmology.lnP[i][cosmology.lnP_nz], io_log10k(i))) 
      {
        cache_update = 1; 
        break;
      }
    }
    for (int j=0; j<cosmology.lnP_nz; j++)
    {
      if (fdiff(cosmology.lnP[cosmology.lnP_nk][j], io_z(j))) 
      {
        cache_update = 1; 
        break;
      }
    }
  }

  jump:

  if (1 == cache_update)
  {
    cosmology.lnP_nk = static_cast<int>(io_log10k.n_elem);
    cosmology.lnP_nz = static_cast<int>(io_z.n_elem);

    if (cosmology.lnP != NULL) free(cosmology.lnP);
    cosmology.lnP = (double**) malloc2d(cosmology.lnP_nk+1,cosmology.lnP_nz+1);

    #pragma omp parallel for
    for (int i=0; i<cosmology.lnP_nk; i++)
      cosmology.lnP[i][cosmology.lnP_nz] = io_log10k(i);
    #pragma omp parallel for
    for (int j=0; j<cosmology.lnP_nz; j++)
      cosmology.lnP[cosmology.lnP_nk][j] = io_z(j);
    #pragma omp parallel for collapse(2)
    for (int i=0; i<cosmology.lnP_nk; i++)
      for (int j=0; j<cosmology.lnP_nz; j++)
        cosmology.lnP[i][j] = io_lnP(i*cosmology.lnP_nz+j);

    cosmology.random = RandomNumber::get_instance().get();
  }

  spdlog::debug("{}: Ends", "set_non_linear_power_spectrum");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
} // namespace cosmolike_interface

