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
#include "cosmolike/bias.h"
#include "cosmolike/IA.h"
#include "cosmolike/interface_aux.hpp"
#include "cosmolike/interface_cosmology.hpp"
#include "cosmolike/interface_nuisance.hpp"
#include "cosmolike/redshift_spline.h"
#include "cosmolike/structs.h"

namespace cosmolike_interface
{

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void init_bias(arma::Col<double> bias_z_evol_model)
{
  spdlog::debug("{}: Begins", "init_bias");
  
  if (MAX_SIZE_ARRAYS < static_cast<int>(bias_z_evol_model.n_elem))
  {
    spdlog::critical("{}: incompatible input {} size = {} (> {})", "init_bias", 
      "bias_z_evol_model", bias_z_evol_model.n_elem, MAX_SIZE_ARRAYS);
    exit(1);
  }

  /*
  int galaxy_bias_model[MAX_SIZE_ARRAYS]; // [0] = b1, 
                                          // [1] = b2, 
                                          // [2] = bs2, 
                                          // [3] = b3, 
                                          // [4] = bmag 
  */
  for(int i=0; i<bias_z_evol_model.n_elem; i++)
  {
    like.galaxy_bias_model[i] = bias_z_evol_model(i);
    
    spdlog::debug("{}: {}[{}] = {} selected.", "init_bias", 
      "like.galaxy_bias_model", i, bias_z_evol_model(i));
  }

  spdlog::debug("{}: Ends", "init_bias");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void init_IA(const int IA_MODEL, const int IA_REDSHIFT_EVOL)
{
  spdlog::debug("{}: Begins", "init_IA");

  spdlog::debug("{}: {} = {} selected.", 
    "init_IA", "IA MODEL", IA_MODEL, 
    "IA REDSHIFT EVOLUTION", IA_REDSHIFT_EVOL);

  if (IA_MODEL == 0 || IA_MODEL == 1)
  {
    like.IA_MODEL = IA_MODEL;
  }
  else
  {
    spdlog::critical("{}: {} = {} not supported", 
      "init_IA", "like.IA_MODEL", IA_MODEL);
    exit(1);
  }

  if (IA_REDSHIFT_EVOL == NO_IA                   || 
      IA_REDSHIFT_EVOL == IA_NLA_LF               ||
      IA_REDSHIFT_EVOL == IA_REDSHIFT_BINNING     || 
      IA_REDSHIFT_EVOL == IA_REDSHIFT_EVOLUTION)
  {
    like.IA = IA_REDSHIFT_EVOL;
  }
  else
  {
    spdlog::critical("{}: {} = {} not supported", 
      "init_IA", "like.IA", IA_REDSHIFT_EVOL);
    exit(1);
  }

  spdlog::debug("{}: Ends", "init_IA");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void set_nuisance_shear_calib(arma::Col<double> M)
{
  spdlog::debug("{}: Begins", "set_nuisance_shear_calib");

  if (redshift.shear_nbin == 0)
  {
    spdlog::critical("{}: {} = 0 is invalid", 
      "set_nuisance_shear_calib", "shear_Nbin");
    exit(1);
  }
  if (redshift.shear_nbin != static_cast<int>(M.n_elem))
  {
    spdlog::critical(
      "{}: incompatible input w/ size = {} (!= {})",
      "set_nuisance_shear_calib", M.n_elem, redshift.shear_nbin);
    exit(1);
  }

  for (int i=0; i<redshift.shear_nbin; i++)
    nuisance.shear_calibration_m[i] = M(i);
  
   spdlog::debug("{}: Ends", "set_nuisance_shear_calib");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void set_nuisance_linear_bias(arma::Col<double> B1)
{
  spdlog::debug("{}: Begins", "set_nuisance_linear_bias");

  if (redshift.clustering_nbin == 0)
  {
    spdlog::critical("{}: {} = 0 is invalid",
      "set_nuisance_linear_bias", "clustering_Nbin");
    exit(1);
  }
  if (redshift.clustering_nbin != static_cast<int>(B1.n_elem))
  {
    spdlog::critical(
      "{}: incompatible input w/ size = {} (!= {})",
      "set_nuisance_linear_bias", B1.n_elem, redshift.clustering_nbin);
    exit(1);
  }

  // GALAXY BIAS ------------------------------------------
  // 1st index: b[0][i] = linear galaxy bias in clustering bin i (b1)
  //            b[1][i] = linear galaxy bias in clustering bin i (b2)
  //            b[2][i] = leading order tidal bias in clustering bin i (b3)
  //            b[3][i] = leading order tidal bias in clustering bin i
  int cache_update = 0;
  for (int i=0; i<redshift.clustering_nbin; i++)
  {
    if(fdiff(nuisance.gb[0][i], B1(i)))
    {
      cache_update = 1;
      nuisance.gb[0][i] = B1(i);
    } 
  }

  if(1 == cache_update)
    nuisance.random_galaxy_bias = RandomNumber::get_instance().get();

  spdlog::debug("{}: Ends", "set_nuisance_linear_bias");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void set_nuisance_nonlinear_bias(arma::Col<double> B1, arma::Col<double> B2)
{
  spdlog::debug("{}: Begins", "set_nuisance_nonlinear_bias");

  if (redshift.clustering_nbin == 0)
  {
    spdlog::critical("{}: {} = 0 is invalid",
      "set_nuisance_nonlinear_bias", "clustering_Nbin");
    exit(1);
  }
  if (redshift.clustering_nbin != static_cast<int>(B1.n_elem) ||
      redshift.clustering_nbin != static_cast<int>(B2.n_elem))
  {
    spdlog::critical(
      "{}: incompatible input w/ sizes = {} and {} (!= {})",
      "set_nuisance_nonlinear_bias", 
      B1.n_elem, B2.n_elem, redshift.clustering_nbin);
    exit(1);
  }

  // GALAXY BIAS ------------------------------------------
  // 1st index: b[0][i]: linear galaxy bias in clustering bin i
  //            b[1][i]: nonlinear b2 galaxy bias in clustering bin i
  //            b[2][i]: leading order tidal bs2 galaxy bias in clustering bin i
  //            b[3][i]: nonlinear b3 galaxy bias  in clustering bin i 
  //            b[4][i]: amplitude of magnification bias in clustering bin i 
  int cache_update = 0;
  for (int i=0; i<redshift.clustering_nbin; i++)
  {
    if(fdiff(nuisance.gb[1][i], B2(i)))
    {
      cache_update = 1;
      nuisance.gb[1][i] = B2(i);
      nuisance.gb[2][i] = almost_equal(B2(i), 0.) ? 0 : (-4./7.)*(B1(i)-1.0);
    }
  }
  if(1 == cache_update)
    nuisance.random_galaxy_bias = RandomNumber::get_instance().get();

  spdlog::debug("{}: Ends", "set_nuisance_nonlinear_bias");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void set_nuisance_magnification_bias(arma::Col<double> B_MAG)
{
  spdlog::debug("{}: Begins", "set_nuisance_magnification_bias");

  if (redshift.clustering_nbin == 0)
  {
    spdlog::critical("{}: {} = 0 is invalid",
      "set_nuisance_magnification_bias", "clustering_Nbin");
    exit(1);
  }
  if (redshift.clustering_nbin != static_cast<int>(B_MAG.n_elem))
  {
    spdlog::critical(
      "{}: incompatible input w/ size = {} (!= {})",
      "set_nuisance_magnification_bias", 
      B_MAG.n_elem, redshift.clustering_nbin);
    exit(1);
  }

  // GALAXY BIAS ------------------------------------------
  // 1st index: b[0][i]: linear galaxy bias in clustering bin i
  //            b[1][i]: nonlinear b2 galaxy bias in clustering bin i
  //            b[2][i]: leading order tidal bs2 galaxy bias in clustering bin i
  //            b[3][i]: nonlinear b3 galaxy bias  in clustering bin i 
  //            b[4][i]: amplitude of magnification bias in clustering bin i
  int cache_update = 0;
  for (int i=0; i<redshift.clustering_nbin; i++)
  {
    if(fdiff(nuisance.gb[4][i], B_MAG(i)))
    {
      cache_update = 1;
      nuisance.gb[4][i] = B_MAG(i);
    }
  }
  if(1 == cache_update)
    nuisance.random_galaxy_bias = RandomNumber::get_instance().get();

  spdlog::debug("{}: Ends", "set_nuisance_magnification_bias");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void set_nuisance_bias(arma::Col<double> B1, arma::Col<double> B2, 
  arma::Col<double> B_MAG)
{
  set_nuisance_linear_bias(B1);
  
  set_nuisance_nonlinear_bias(B1, B2);
  
  set_nuisance_magnification_bias(B_MAG);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void set_nuisance_IA(arma::Col<double> A1, arma::Col<double> A2, 
  arma::Col<double> BTA)
{
  spdlog::debug("{}: Begins", "set_nuisance_IA");

  if (redshift.shear_nbin == 0)
  {
    spdlog::critical("{}: {} = 0 is invalid",
      "set_nuisance_IA", "shear_Nbin");
    exit(1);
  }
  if (redshift.shear_nbin > static_cast<int>(A1.n_elem) ||
      redshift.shear_nbin > static_cast<int>(A2.n_elem) ||
      redshift.shear_nbin > static_cast<int>(BTA.n_elem))
  {
    spdlog::critical(
      "{}: incompatible input w/ sizes = {}, {} and {} (!= {})",
      "set_nuisance_IA", A1.n_elem, A2.n_elem, BTA.n_elem, redshift.shear_nbin);
    exit(1);
  }

  // INTRINSIC ALIGMENT ------------------------------------------  
  // ia[0][0] = A_ia          if(IA_NLA_LF || IA_REDSHIFT_EVOLUTION)
  // ia[0][1] = eta_ia        if(IA_NLA_LF || IA_REDSHIFT_EVOLUTION)
  // ia[0][2] = eta_ia_highz  if(IA_NLA_LF, Joachimi2012)
  // ia[0][3] = beta_ia       if(IA_NLA_LF, Joachimi2012)
  // ia[0][4] = LF_alpha      if(IA_NLA_LF, Joachimi2012)
  // ia[0][5] = LF_P          if(IA_NLA_LF, Joachimi2012)
  // ia[0][6] = LF_Q          if(IA_NLA_LF, Joachimi2012)
  // ia[0][7] = LF_red_alpha  if(IA_NLA_LF, Joachimi2012)
  // ia[0][8] = LF_red_P      if(IA_NLA_LF, Joachimi2012)
  // ia[0][9] = LF_red_Q      if(IA_NLA_LF, Joachimi2012)
  // ------------------
  // ia[1][0] = A2_ia        if IA_REDSHIFT_EVOLUTION
  // ia[1][1] = eta_ia_tt    if IA_REDSHIFT_EVOLUTION
  // ------------------
  // ia[2][MAX_SIZE_ARRAYS] = b_ta_z[MAX_SIZE_ARRAYS]

  int cache_update = 0;
  nuisance.c1rhocrit_ia = 0.01389;
  
  if (like.IA == IA_REDSHIFT_BINNING)
  {
    for (int i=0; i<redshift.shear_nbin; i++)
    {
      if (fdiff(nuisance.ia[0][i],A1(i)) ||
          fdiff(nuisance.ia[1][i],A2(i)) ||
          fdiff(nuisance.ia[2][i],A2(i)))
      {
        nuisance.ia[0][i] = A1(i);
        nuisance.ia[1][i] = A2(i);
        nuisance.ia[2][i] = BTA(i);
        cache_update = 1;
      }
    }
  }
  else if (like.IA == IA_REDSHIFT_EVOLUTION)
  {
    nuisance.oneplusz0_ia = 1.62;
    if (fdiff(nuisance.ia[0][0],A1(0)) ||
        fdiff(nuisance.ia[0][1],A1(1)) ||
        fdiff(nuisance.ia[1][0],A2(0)) ||
        fdiff(nuisance.ia[1][1],A2(1)) ||
        fdiff(nuisance.ia[2][0],BTA(0)))
    {
      nuisance.ia[0][0] = A1(0);
      nuisance.ia[0][1] = A1(1);
      nuisance.ia[1][0] = A2(0);
      nuisance.ia[1][1] = A2(1);
      nuisance.ia[2][0] = BTA(0);
      cache_update = 1;
    }
  }
  if(1 == cache_update)
    nuisance.random_ia = RandomNumber::get_instance().get();

  spdlog::debug("{}: Ends", "set_nuisance_ia");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
} // namespace cosmolike_interface