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

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void init_cmb_bandpower(const int is_bandpower, 
    const int is_cmb_kkkk_cov_from_sim, const double alpha)
{
  spdlog::debug("{}: Begins", "init_cmb_bandpower");

  IPCMB::get_instance().is_bandpower_     = is_bandpower;
  IPCMB::get_instance().is_cov_kkkk_term_from_sim_ = is_cmb_kkkk_cov_from_sim;
  IPCMB::get_instance().alpha_hartlap_cov_kkkk_ = alpha;

  spdlog::debug("{}: Ends", "init_cmb_bandpower");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void init_binning_cmb_bandpower(const int nbp, const int lmin, const int lmax)
{
  spdlog::debug("{}: Begins", "init_binning_cmb_bandpower");

  if (!(nbp > 0))
  {
    spdlog::critical("{}: {} = {} not supported", 
      "init_binning_cmb_bandpower", "Number of Band Powers", nbp);
    exit(1);
  }

  spdlog::debug("{}: {} = {} selected.","init_binning_cmb_bandpower","NBandPower",nbp);
  spdlog::debug("{}: {} = {} selected.","init_binning_cmb_bandpower","lmin",lmin);
  spdlog::debug("{}: {} = {} selected.","init_binning_cmb_bandpower","lmax",lmax);

  IPCMB::get_instance().nbandpower_ = nbandpower;
  IPCMB::get_instance().lrange_bandpower_[0] = lmin;
  IPCMB::get_instance().lrange_bandpower_[1] = lmax;
  
  spdlog::debug("{}: Ends", "init_binning_cmb_bandpower");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void init_cmb(const double lmin_kappa_cmb, const double lmax_kappa_cmb, 
    const double fwhm, std::string pixwin_file) 
{
  spdlog::debug("{}: Begins", "init_cmb");

  like.lmin_kappacmb = lmin_kappa_cmb;
  
  like.lmax_kappacmb = lmax_kappa_cmb;
  
  // fwhm = beam size in arcmin - cmb.fwhm = beam size in rad
  cmb.fwhm = fwhm*2.90888208665721580e-4;
  
  memcpy(cmb.pathHealpixWinFunc, pixwin_file.c_str(), pixwin_file.size()+1);
  
  spdlog::debug("{}: Ends", "init_cmb");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void init_cmb_bandpower_data(std::string binning_matrix, std::string theory_offset)
{
  spdlog::debug("{}: Begins", "init_cmb_bandpower_data");

  IPCMB::get_instance().set_binning_mat(binning_matrix);
  IPCMB::get_instance().set_theory_offset(theory_offset);
  
  spdlog::debug("{}: Ends", "init_cmb_bandpower_data");

  return;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
} // namespace cosmolike_interface