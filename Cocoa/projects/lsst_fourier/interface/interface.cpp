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

#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/cfg/env.h>

#include <boost/algorithm/string.hpp>

#ifdef PYBIND11
// Python Binding
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace py = pybind11;
#endif

#include "cosmolike/basics.h"
#include "cosmolike/bias.h"
#include "cosmolike/cosmo2D_fourier.h"
#include "cosmolike/cosmo2D_exact_fft.h"
#include "cosmolike/cosmo3D.h"
#include "cosmolike/halo.h"
#include "cosmolike/radial_weights.h"
#include "cosmolike/recompute.h"
#include "cosmolike/pt_cfastpt.h"
#include "cosmolike/redshift_spline.h"
#include "cosmolike/structs.h"

#include "interface.hpp"

namespace ima = interface_mpp_aux;

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// init functions
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

void cpp_initial_setup() {
  spdlog::cfg::load_env_levels();
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "initial_setup");

  // restart variables to 0 so error check can flag bad initialization
  tomo.shear_Nbin = 0;
  tomo.clustering_Nbin = 0;

  like.shearcalib = 0; // on cpp_compute_data_vector I erased shearcalib flag
  like.shear_shear = 0;
  like.shear_pos = 0;
  like.pos_pos = 0;

  // bias
  gbias.b1_function = &b1_per_bin;

  // no priors
  like.clusterN = 0;
  like.clusterWL = 0;
  like.clusterCG = 0;
  like.clusterCC = 0;
  like.SN_WFIRST = 0;
  like.BAO = 0;
  like.wlphotoz = 0;
  like.clphotoz = 0;
  like.Planck15_BAO_H070p6_JLA_w0wa = 0;
  like.Planck18_BAO_Riess18_Pantheon_w0wa = 0;
  like.Planck18_BAO_w0wa = 0;
  like.Planck18_w0 = 0;

  // reset bias
  for (int i = 0; i < MAX_SIZE_ARRAYS; i++) {
    gbias.b[i] = 0.0;
    gbias.b2[i] = 0.0;
    gbias.b_mag[i] = 0.0;
  }
  // reset IA
  for (int i = 0; i < MAX_SIZE_ARRAYS; i++) {
    nuisance.A_z[i] = 0.0;
    nuisance.A2_z[i] = 0.0;
    nuisance.b_ta_z[i] = 0.0;
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "initial_setup");
}

void cpp_init_probes(std::string possible_probes) {
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_probes");

  if (possible_probes.compare("ss") == 0) { // cosmolike c interface
    like.shear_shear = 1;
    spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected",
      "init_probes", "possible_probes", "ss");
  } else if (possible_probes.compare("gg") == 0) {
    like.pos_pos = 1;
    spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected",
      "init_probes", "possible_probes", "gg");
  } else if (possible_probes.compare("gs") == 0) {
    like.shear_pos = 1;
    spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected",
      "init_probes", "possible_probes", "gs");
  } else if (possible_probes.compare("2x2pt") == 0) {
    like.shear_pos = 1;
    like.pos_pos = 1;
    spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected",
      "init_probes", "possible_probes", "2x2pt");
  } else if (possible_probes.compare("3x2pt") == 0) {
    like.shear_shear = 1;
    like.shear_pos = 1;
    like.pos_pos = 1;
    spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected", "init_probes", "possible_probes",
      "3x2pt");
  } else if (possible_probes.compare("ss_gs") == 0) {
    like.shear_shear = 1;
    like.shear_pos = 1;
    spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected", "init_probes", "possible_probes",
      "ss + gs (2x2pt)");
  } else {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = {} probe not supported",
      "init_probes", "possible_probes", possible_probes);
    exit(1);
  }
  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_probes");
}

void cpp_init_survey(std::string surveyname, double area, double sigma_e) {
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_survey");
#ifdef DEBUG
  if (surveyname.size() > CHAR_MAX_SIZE - 1) {

    exit(1);
  }
  if (!(surveyname.size()>0)) {
    spdlog::critical("{}: incompatible input", "init_survey");
    exit(1);
  }
#endif

  memcpy(survey.name, surveyname.c_str(), surveyname.size()+1);
  survey.area = area;
  survey.sigma_e = sigma_e;

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_survey");
}

void cpp_init_cosmo_runmode(const bool is_linear) {
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_cosmo_runmode");

  std::string mode = is_linear ? "linear" : "Halofit";
  const size_t size = mode.size();
  memcpy(pdeltaparams.runmode, mode.c_str(), size+1);

  spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected", "init_cosmo_runmode", "runmode", mode);
  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_cosmo_runmode");
}

void cpp_init_IA(int N) {
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_IA");
  spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected.", "init_IA", "IA", N);

  if (N == 3 || N == 4 || N == 5 || N == 6) {
    like.IA = N;
  } else {
    spdlog::critical("{}: {} = {} not supported", "init_IA", "like.IA", N);
    exit(1);
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_IA");
}


void cpp_init_binning(const int Ncl, const double lmin,
const double lmax) {
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_binning");
#ifdef DEBUG
  if (!(Ncl > 0)) {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = {} not supported",
      "init_binning", "like.Ncl", Ncl);
    exit(1);
  }
#endif
  spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected.", "init_binning",
    "Ncl", Ncl);
  spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected.", "init_binning",
    "l_min", lmin);
  spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected.", "init_binning",
    "l_max", lmax);

  like.Ncl = Ncl;
  like.lmin = lmin;
  like.lmax = lmax;
  const double logdl = (std::log(like.lmax) - std::log(like.lmin))/like.Ncl;
  like.ell = create_double_vector(0, like.Ncl-1);
  for (int i = 0; i < like.Ncl; i++) {
    like.ell[i] = std::exp(std::log(like.lmin)+(i+0.5)*logdl);
    spdlog::debug(
     "\x1b[90m{}\x1b[0m: Bin {:d} - {} = {:.4e}, {} = {:.4e} and {} = {:.4e}",
     "init_binning", i, "lmin", lmin, "ell",
     like.ell[i], "lmax", lmax);
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_binning");
}

void cpp_init_lens_sample(std::string multihisto_file, const int Ntomo,
const double ggl_cut) {
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_lens_sample");
#ifdef DEBUG
  if (tomo.shear_Nbin == 0) {
    spdlog::critical("{}: {} not set prior to this function call",
      "init_lens_sample", "tomo.shear_Nbin");
    exit(1);
  }
  if (multihisto_file.size()>CHAR_MAX_SIZE-1) {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: insufficient pre-allocated char memory (max = {}) for"
      "the string: {}", "init_lens_sample", CHAR_MAX_SIZE-1, multihisto_file);
    exit(1);
  }
  if (!(multihisto_file.size() > 0)) {
    spdlog::critical("\x1b[90m{}\x1b[0m: empty {} string not supported",
      "init_lens_sample", "multihisto_file");
    exit(1);
  }
  if (!(Ntomo > 0) || Ntomo > MAX_SIZE_ARRAYS) {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = {} not supported (max = {})",
      "init_lens_sample", "Ntomo", Ntomo, MAX_SIZE_ARRAYS);
    exit(1);
  }
#endif
  memcpy(redshift.clustering_REDSHIFT_FILE, multihisto_file.c_str(),
    multihisto_file.size()+1);

  redshift.clustering_photoz = 4;
  tomo.clustering_Nbin = Ntomo;
  tomo.clustering_Npowerspectra = tomo.clustering_Nbin;

  spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected.", "init_lens_sample",
    "clustering_REDSHIFT_FILE", multihisto_file);

  spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected.", "init_lens_sample",
    "clustering_Nbin", Ntomo);

  if (ggl_cut > 0) {
    survey.ggl_overlap_cut = ggl_cut;
  }
  else {
    survey.ggl_overlap_cut = 0.0;
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected.", "init_lens_sample",
    "survey.ggl_overlap_cut", survey.ggl_overlap_cut);

  pf_photoz(0.1, 0);
  {
    int n = 0;
    for (int i = 0; i < tomo.clustering_Nbin; i++) {
      for (int j = 0; j < tomo.shear_Nbin; j++) {
        n += test_zoverlap(i, j);
      }
    }
    tomo.ggl_Npowerspectra = n;
    spdlog::debug("\x1b[90m{}\x1b[0m: tomo.ggl_Npowerspectra = {}",
      "init_lens_sample", tomo.ggl_Npowerspectra);
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_lens_sample");
}

void cpp_init_source_sample(std::string multihisto_file, const int Ntomo) {
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_source_sample");
#ifdef DEBUG
  if (multihisto_file.size() > CHAR_MAX_SIZE - 1) {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: insufficient pre-allocated char memory (max = {}) for"
      "the string: {}", "init_source_sample", CHAR_MAX_SIZE-1, multihisto_file);
    exit(1);
  }
  if (!(multihisto_file.size() > 0)) {
    spdlog::critical("\x1b[90m{}\x1b[0m: empty {} string not supported",
      "init_source_sample", "multihisto_file");
    exit(1);
  }
  if (!(Ntomo > 0) || Ntomo > MAX_SIZE_ARRAYS) {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = {} not supported (max = {})",
      "init_source_sample", "Ntomo", Ntomo, MAX_SIZE_ARRAYS);
    exit(1);
  }
#endif
  // convert std::string to char*
  memcpy(redshift.shear_REDSHIFT_FILE, multihisto_file.c_str(),
    multihisto_file.size()+1);

  redshift.shear_photoz = 4;
  tomo.shear_Nbin = Ntomo;
  tomo.shear_Npowerspectra = tomo.shear_Nbin * (tomo.shear_Nbin + 1) / 2;

  spdlog::debug("\x1b[90m{}\x1b[0m: tomo.shear_Npowerspectra = {}",
    "init_source_sample", tomo.shear_Npowerspectra);

  for (int i = 0; i < tomo.shear_Nbin; i++) {
    nuisance.bias_zphot_shear[i] = 0.0;
    spdlog::info("\x1b[90m{}\x1b[0m: bin {} - {} = {}.",
      "init_source_sample", i, "<z_s>", zmean_source(i));
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected.", "init_source_sample",
    "shear_REDSHIFT_FILE", multihisto_file);
  spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected.", "init_source_sample",
    "shear_Nbin", Ntomo);
  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_source_sample");
}

void cpp_init_size_data_vector() {
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_size_data_vector");
#ifdef DEBUG
  if (tomo.shear_Nbin == 0) {
    spdlog::critical("{}: {} not set prior to this function call",
      "init_size_data_vector", "tomo.shear_Nbin");
    exit(1);
  }
  if (tomo.clustering_Nbin == 0) {
    spdlog::critical("{}: {} not set prior to this function call",
      "init_size_data_vector", "tomo.clustering_Nbin");
    exit(1);
  }
  if (like.Ncl == 0) {
    spdlog::critical("{}: {} not set prior to this function call",
      "init_size_data_vector", "like.Ncl");
    exit(1);
  }
#endif

  like.Ndata = like.Ncl*(tomo.shear_Npowerspectra +
                         tomo.ggl_Npowerspectra + tomo.clustering_Npowerspectra);

  spdlog::debug("\x1b[90m{}\x1b[0m: {} = {} selected.", "init_size_data_vector", "Ndata",
    like.Ndata);
  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_size_data_vector");
}

void cpp_init_linear_power_spectrum(std::vector<double> io_log10k,
std::vector<double> io_z, std::vector<double> io_lnP) {
#ifdef DEBUG
  {
    bool debug_fail = false;
    if (io_z.size()*io_log10k.size() != io_lnP.size()) {
      debug_fail = true;
    }
    else{
      if (io_z.size() == 0) {
        debug_fail = true;
      }
    }
    if (debug_fail) {
      spdlog::critical("\x1b[90m{}\x1b[0m: incompatible input w/ k.size = {}, z.size = {}, "
        "and lnP.size = {}", "init_linear_power_spectrum", io_log10k.size(),
        io_z.size(), io_lnP.size());
      exit(1);
    }
  }
#endif
  int nlog10k = static_cast<int>(io_log10k.size());
  int nz = static_cast<int>(io_z.size());
  double* log10k = io_log10k.data();
  double* z = io_z.data();
  double* lnP = io_lnP.data();
  setup_p_lin(&nlog10k, &nz, &log10k, &z, &lnP, 1);
  // force initialization - imp to avoid seg fault when openmp is on
  const double io_a = 1.0;
  const double io_k = 0.1*cosmology.coverH0;
  p_lin(io_k, io_a);
}

void cpp_init_non_linear_power_spectrum(std::vector<double> io_log10k,
  std::vector<double> io_z, std::vector<double> io_lnP) {
#ifdef DEBUG
  {
    bool debug_fail = false;
    if (io_z.size()*io_log10k.size() != io_lnP.size()) {
      debug_fail = true;
    }
    else{
      if (io_z.size() == 0) {
        debug_fail = true;
      }
    }
    if (debug_fail) {
      spdlog::critical("\x1b[90m{}\x1b[0m: incompatible input w/ k.size = {}, z.size = {}, "
        "and lnP.size = {}", "init_non_linear_power_spectrum", io_log10k.size(),
        io_z.size(), io_lnP.size());
      exit(1);
    }
  }
#endif
  int nlog10k = static_cast<int>(io_log10k.size());
  int nz = static_cast<int>(io_z.size());
  double* log10k = io_log10k.data();
  double* z = io_z.data();
  double* lnP = io_lnP.data();
  setup_p_nonlin(&nlog10k, &nz, &log10k, &z, &lnP, 1);
  // force initialization - imp to avoid seg fault when openmp is on
  const double io_a = 1.0;
  const double io_k = 0.1*cosmology.coverH0;
  p_nonlin(io_k, io_a);
}

// Growth: D = G * a
void cpp_init_growth(std::vector<double> io_z, std::vector<double> io_G) {
#ifdef DEBUG
  {
    bool debug_fail = false;
    if (io_z.size() != io_G.size()) {
      debug_fail = true;
    }
    else{
      if (io_z.size() == 0) {
        debug_fail = true;
      }
    }
    if (debug_fail) {
      spdlog::critical("\x1b[90m{}\x1b[0m: incompatible input w/ z.size = {} and G.size = {}",
        "init_growth", io_z.size(), io_G.size());
      exit(1);
    }
  }
#endif
  int nz = static_cast<int>(io_z.size());
  double* z = io_z.data();
  double* G = io_G.data();
  setup_growth(&nz, &z, &G, 1);
  // force initialization - imp to avoid seg fault when openmp is on
  const double io_a = 1.0;
  const double zz = 0.0;
  f_growth(zz);
  growfac_all(io_a);
  growfac(io_a);
  return;
}

void cpp_init_distances(std::vector<double> io_z, std::vector<double> io_chi ) {
#ifdef DEBUG
  {
    bool debug_fail = false;
    if (io_z.size() != io_chi.size()) {
      debug_fail = true;
    }
    else{
      if (io_z.size() == 0) {
        debug_fail = true;
      }
    }
    if (debug_fail) {
      spdlog::critical("\x1b[90m{}\x1b[0m: incompatible input w/ z.size = {} and G.size = {}",
        "init_distances", io_z.size(), io_chi.size());
      exit(1);
    }
  }
#endif
  int nz = static_cast<int>(io_z.size());
  double* vz = io_z.data();
  double* vchi = io_chi.data();
  setup_chi(&nz, &vz, &vchi, 1);
  // force initialization - imp to avoid seg fault when openmp is on
  const double io_a = 1.0;
  chi(io_a);
}


void cpp_init_data_real(std::string COV, std::string MASK, std::string DATA) {
  ima::RealData& instance = ima::RealData::get_instance();
  instance.set_data(DATA);
  instance.set_mask(MASK);
  instance.set_inv_cov(COV);
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// SET PARAM FUNCTIONS
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

void cpp_set_cosmological_parameters(const double omega_matter,
const double hubble, const bool is_cached_cosmology) {
  if(!is_cached_cosmology) {
    // Cosmolike should not need parameters from inflation or dark energy.
    // because Cobaya provides P(k,z), H(z), D(z), Chi(z)...
    // It may require H0 to set scales and \Omega_M to set the halo model

    // cosmolike c interface
    cosmology.Omega_m = omega_matter;
    cosmology.Omega_v = 1.0-omega_matter;
    // Cosmolike only needs to know that there are massive neutrinos (>0)
    cosmology.Omega_nu = 0.1;
    cosmology.h0 = hubble/100.0; // assuming H0 in km/s/Mpc
    cosmology.MGSigma = 0.0;
    cosmology.MGmu = 0.0;

    // Technical Problem: we want Cosmolike to calculate the data vector when
    // Cobaya request (no cache). To avoid cache in Cosmolike, we use a
    // random number generators to set cosmology.random
    cosmology.random = ima::RandomNumber::get_instance().get();
    cosmology.is_cached = 0;
  } else {
    cosmology.is_cached = 1;
  }
}

void cpp_set_nuisance_shear_calib(std::vector<double> M) {
#ifdef DEBUG
  if (tomo.shear_Nbin == 0) {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = 0 is invalid", "set_nuisance_shear_calib",
      "shear_Nbin");
    exit(1);
  }
  if (tomo.shear_Nbin != static_cast<int>(M.size())) {
    spdlog::critical("\x1b[90m{}\x1b[0m: incompatible input w/ size = {} (!= {})",
        "set_nuisance_shear_calib", M.size(), tomo.shear_Nbin);
    exit(1);
  }
#endif
  for (int i = 0; i < tomo.shear_Nbin; i++) {
    nuisance.shear_calibration_m[i] = M[i];
  }
}

void cpp_set_nuisance_shear_photoz(std::vector<double> SP) {
#ifdef DEBUG
  if (tomo.shear_Nbin == 0) {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = 0 is invalid", "set_nuisance_shear_photoz",
      "shear_Nbin");
    exit(1);
  }
  if (tomo.shear_Nbin != static_cast<int>(SP.size())) {
    spdlog::critical("\x1b[90m{}\x1b[0m: incompatible input w/ size = {} (!= {})",
        "set_nuisance_shear_photoz", SP.size(), tomo.shear_Nbin);
    exit(1);
  }
#endif
  for (int i = 0; i < tomo.shear_Nbin; i++) {
    nuisance.bias_zphot_shear[i] = SP[i];
  }
}

void cpp_set_nuisance_clustering_photoz(std::vector<double> CP) {
#ifdef DEBUG
  if (tomo.clustering_Nbin == 0) {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = 0 is invalid", "set_nuisance_clustering_photoz",
      "clustering_Nbin");
    exit(1);
  }
  if (tomo.clustering_Nbin != static_cast<int>(CP.size())) {
    spdlog::critical("\x1b[90m{}\x1b[0m: incompatible input w/ size = {} (!= {})",
        "set_nuisance_clustering_photoz", CP.size(), tomo.clustering_Nbin);
    exit(1);
  }
#endif
  for (int i = 0; i < tomo.clustering_Nbin; i++) {
    nuisance.bias_zphot_clustering[i] = CP[i];
  }
}

void cpp_set_nuisance_linear_bias(std::vector<double> B1) {
#ifdef DEBUG
  if (tomo.clustering_Nbin == 0) {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = 0 is invalid", "set_nuisance_linear_bias",
      "clustering_Nbin");
    exit(1);
  }
  if (tomo.clustering_Nbin != static_cast<int>(B1.size())) {
    spdlog::critical("\x1b[90m{}\x1b[0m: incompatible input w/ size = {} (!= {})",
        "set_nuisance_linear_bias", B1.size(), tomo.clustering_Nbin);
    exit(1);
  }
#endif
  for (int i = 0; i < tomo.clustering_Nbin; i++) {
    gbias.b[i] = B1[i];
  }
}

void cpp_set_nuisance_nonlinear_bias(std::vector<double> B1, std::vector<double> B2) {
#ifdef DEBUG
  if (tomo.clustering_Nbin == 0) {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = 0 is invalid", "set_nuisance_nonlinear_bias",
      "clustering_Nbin");
    exit(1);
  }
  if (tomo.clustering_Nbin != static_cast<int>(B1.size()) ||
      tomo.clustering_Nbin != static_cast<int>(B2.size())) {
    spdlog::critical("\x1b[90m{}\x1b[0m: incompatible input w/ sizes = {} and {} (!= {})",
        "set_nuisance_nonlinear_bias", B1.size(), B2.size(),
        tomo.clustering_Nbin);
    exit(1);
  }
#endif
  constexpr double tmp = -4./7.;
  for (int i = 0; i < tomo.clustering_Nbin; i++) {
    gbias.b2[i] = B2[i];
    gbias.bs2[i] = ima::almost_equal(B2[i], 0.) ? 0 : tmp*(B1[i]-1.0);
  }
}

void cpp_set_nuisance_magnification_bias(std::vector<double> B_MAG) {
#ifdef DEBUG
  if (tomo.clustering_Nbin == 0) {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = 0 is invalid", "set_nuisance_magnification_bias",
      "clustering_Nbin");
    exit(1);
  }
  if (tomo.clustering_Nbin != static_cast<int>(B_MAG.size())) {
    spdlog::critical("\x1b[90m{}\x1b[0m: incompatible input w/ size = {} (!= {})",
        "set_nuisance_magnification_bias", B_MAG.size(), tomo.clustering_Nbin);
    exit(1);
  }
#endif
  for (int i = 0; i < tomo.clustering_Nbin; i++) {
    gbias.b_mag[i] = B_MAG[i];
  }
}

void cpp_set_nuisance_bias(std::vector<double> B1, std::vector<double> B2,
  std::vector<double> B_MAG) {
  cpp_set_nuisance_linear_bias(B1);
  cpp_set_nuisance_nonlinear_bias(B1, B2);
  cpp_set_nuisance_magnification_bias(B_MAG);
}

void cpp_set_nuisance_ia_mpp(std::vector<double> A1, std::vector<double> A2,
  std::vector<double> B_TA) {
#ifdef DEBUG
  if (tomo.shear_Nbin == 0) {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = 0 is invalid",
      "set_nuisance_ia_mpp", "shear_Nbin");
    exit(1);
  }
  if (tomo.shear_Nbin != static_cast<int>(A1.size()) ||
      tomo.shear_Nbin != static_cast<int>(A2.size()) ||
  		tomo.shear_Nbin != static_cast<int>(B_TA.size())) {
    spdlog::critical("\x1b[90m{}\x1b[0m: incompatible input w/ sizes = {}, {} and {} (!= {})",
        "set_nuisance_ia_mpp", A1.size(), A2.size(), B_TA.size(),
        tomo.shear_Nbin);
    exit(1);
  }
#endif
  nuisance.c1rhocrit_ia = 0.01389;
  if (like.IA == 3 || like.IA == 5) {
    for (int i = 0; i < tomo.shear_Nbin; i++) {
      nuisance.A_z[i] = A1[i];
      nuisance.A2_z[i] = A2[i];
      nuisance.b_ta_z[i] = B_TA[i];
    }
  }
  else if (like.IA == 4 || like.IA == 6) {
    nuisance.A_ia = A1[0];
    nuisance.eta_ia = A1[1];
    nuisance.oneplusz0_ia = 1.62;
    nuisance.A2_ia = A2[0];
    nuisance.eta_ia_tt = A2[1];
    nuisance.b_ta_z[0] = B_TA[0];
#ifdef DEBUG
    for (int i = 2; i < tomo.shear_Nbin; i++) {
      if ( !(ima::almost_equal(A1[i],0.)) ||
           !(ima::almost_equal(A2[i],0.)) ||
           !(ima::almost_equal(B_TA[i],0.))) {
        spdlog::critical(
        	"set_nuisance_ia_mpp: one of nuisance.A_z[{}]={}, nuisance.A2_z[{}]="
          "{}, nuisance.b_ta[{}]={} was specified w/ power-law evolution\n",
          i, nuisance.A_z[i], i, nuisance.A2_z[i], i, nuisance.b_ta_z[i]
        );
        exit(1);
      }
    }
#endif
  }
}


// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// COMPUTE FUNCTIONS
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double cpp_compute_chi2(std::vector<double> datavector) {
  ima::RealData& instance = ima::RealData::get_instance();
  return instance.get_chi2(datavector);
}

int cpp_compute_mask(const int i) {
  ima::RealData& instance = ima::RealData::get_instance();
  return instance.get_mask(i);
}

std::vector<double> cpp_compute_data_vector() {
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "compute_data_vector");
#ifdef DEBUG
  if (tomo.shear_Nbin == 0) {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = 0 is invalid",
      "compute_data_vector", "shear_Nbin");
    exit(1);
  }
  if (tomo.clustering_Nbin == 0) {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = 0 is invalid",
      "compute_data_vector", "clustering_Nbin");
    exit(1);
  }
  if (like.Ncl == 0) {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} = 0 is invalid",
      "compute_data_vector", "Ncl");
    exit(1);
  }
  if (!ima::RealData::get_instance().is_mask_set()) {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} not set prior to this function call",
      "compute_data_vector", "mask");
    exit(1);
  }
  if (!ima::RealData::get_instance().is_data_set()) {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} not set prior to this function call",
      "compute_data_vector", "data_vector");
    exit(1);
  }
  if (!ima::RealData::get_instance().is_inv_cov_set()) {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} not set prior to this function call",
      "compute_data_vector", "inv_cov");
    exit(1);
  }
#endif
  std::vector<double> data_vector(like.Ndata,0.0);
  int start = 0;
  if (like.shear_shear == 1) {
    for (int nz = 0; nz<tomo.shear_Npowerspectra; nz++) {
      const int z1 = Z1(nz);
      const int z2 = Z2(nz);
      for (int i = 0; i<like.Ncl; i++) {
        if (cpp_compute_mask(like.Ncl*nz+i)) {
          data_vector[like.Ncl*nz+i] =
            C_shear_shear_IA_tab(like.ell[i],z1,z2)*
            (1.0 + nuisance.shear_calibration_m[z1])*
            (1.0 + nuisance.shear_calibration_m[z2]);
        }
      }
    }
  }
  start = start + like.Ncl*tomo.shear_Npowerspectra;
  if (like.shear_pos == 1) {
    for (int nz=0; nz<tomo.ggl_Npowerspectra; nz++) {
      const int zl = ZL(nz);
      const int zs = ZS(nz);
      for (int i=0; i<like.Ncl; i++) {
        if (cpp_compute_mask(start+(like.Ncl*nz)+i)) {
          data_vector[start+(like.Ncl*nz)+i] = C_ggl_IA_tab(like.ell[i],zl,zs)
            *(1.0+nuisance.shear_calibration_m[zs]);
        }
      }
    }
  }
  start = start + like.Ncl*tomo.ggl_Npowerspectra;
  if (like.pos_pos == 1) {
    {
      const int nz = 0;
      {
        const int i = 0;
        if (cpp_compute_mask(start+(like.Ncl*nz)+i)) {
          data_vector[start+(like.Ncl*nz)+i] = C_cl_tomo_nointerp(like.ell[i],nz,nz);
        }
      }
      #pragma omp parallel for
      for (int i=1; i<like.Ncl; i++) {
        if (cpp_compute_mask(start+(like.Ncl*nz)+i)) {
          data_vector[start+(like.Ncl*nz)+i] = C_cl_tomo_nointerp(like.ell[i],nz,nz);
        }
      }
    }
    #pragma omp parallel for
    for (int nz=1; nz<tomo.clustering_Npowerspectra; nz++) {
      for (int i=0; i<like.Ncl; i++) {
        if (cpp_compute_mask(start+(like.Ncl*nz)+i)) {
          data_vector[start+(like.Ncl*nz)+i] = C_cl_tomo_nointerp(like.ell[i],nz,nz);
        }
      }
    }
  }
  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "compute_data_vector");
  return data_vector;
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// CLASS RealData MEMBER FUNCTIONS (& RELATED) - READ MASK, COV..
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

arma::Mat<double> ima::read_table(const std::string file_name) {
  std::ifstream input_file(file_name);
#ifdef DEBUG
  if (!input_file.is_open()) {
    spdlog::critical("\x1b[90m{}\x1b[0m: file {} cannot be opened", "read_table", file_name);
    exit(1);
  }
#endif
  // Read the entire file into memory
  std::string tmp;
  input_file.seekg(0,std::ios::end);
  tmp.resize(static_cast<size_t>(input_file.tellg()));
  input_file.seekg(0,std::ios::beg);
  input_file.read(&tmp[0],tmp.size());
  input_file.close();
#ifdef DEBUG
  if(tmp.empty()) {
    spdlog::critical("\x1b[90m{}\x1b[0m: file {} is empty", "read_table", file_name);
    exit(1);
  }
#endif
  std::vector<std::string> lines;
  lines.reserve(50000);
  // Second: Split file into lines
  boost::trim_if(tmp,boost::is_any_of("\t "));
  boost::trim_if(tmp,boost::is_any_of("\n"));
  boost::split(lines, tmp,boost::is_any_of("\n"), boost::token_compress_on);
  // Erase comment/blank lines
  auto check = [](std::string mystr) -> bool {
    return boost::starts_with(mystr, "#");
  };
  lines.erase(std::remove_if(lines.begin(), lines.end(), check), lines.end());
  // Third: Split line into words
  arma::Mat<double> result;
  size_t ncols = 0;
  { // first line
    std::vector<std::string> words;
    words.reserve(100);
    boost::split(words,lines[0],boost::is_any_of(" \t"), boost::token_compress_on);
    ncols = words.size();
    result.set_size(lines.size(), ncols);
    for (size_t j=0; j<ncols; j++) {
      result(0,j) = std::stod(words[j]);
    }
  }
  #pragma omp parallel for
  for (size_t i=1; i<lines.size(); i++) {
    std::vector<std::string> words;
    boost::split(words, lines[i], boost::is_any_of(" \t"),
      boost::token_compress_on);
    if (words.size() != ncols) {
      spdlog::critical("\x1b[90m{}\x1b[0m: file {} is not well formatted"
      " (regular table required)", "read_table", file_name);
      exit(1);
    }
    for (size_t j=0; j<ncols; j++) {
      result(i,j) = std::stod(words[j]);
    }
  };
  return result;
}

void ima::RealData::set_data(std::string DATA) {
#ifdef DEBUG
  if(!(like.Ndata>0)) {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} not set prior to this function call",
      "set_data", "like.Ndata");
    exit(1);
  }
#endif
  this->data_.set_size(like.Ndata);
  arma::Mat<double> table = ima::read_table(DATA);
  for (int i=0; i<like.Ndata; i++) {
    this->data_(i) = table(i,1);
  }
  this->data_filename_ = DATA;
  this->is_data_set_ = true;
}

void ima::RealData::set_mask(std::string MASK) {
#ifdef DEBUG
  if (!(like.Ndata>0)) {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} not set prior to this function call",
      "set_mask", "like.Ndata");
    exit(1);
  }
#endif
  this->mask_.set_size(like.Ndata);
  arma::Mat<double> table = ima::read_table(MASK);
  for (int i=0; i<like.Ndata; i++) {
    this->mask_(i) = table(i,1);
  }
  if (like.shear_shear==0) {
    const int M = like.Ncl*tomo.shear_Npowerspectra;
    for (int i=0; i<M; i++) {
      this->mask_(i) = 0.0;
    }
  }
  if (like.shear_pos==0) {
    const int N = like.Ncl*tomo.shear_Npowerspectra;
    const int M = N + like.Ncl*tomo.ggl_Npowerspectra;
    for (int i=N; i<M; i++) {
      this->mask_(i) = 0.0;
    }
  }
  if (like.pos_pos==0) {
    const int N = like.Ncl*(tomo.shear_Npowerspectra +
      tomo.ggl_Npowerspectra);
    const int M = N + like.Ncl*tomo.clustering_Npowerspectra;
    for (int i=N; i<M; i++) {
      this->mask_(i) = 0.0;
    }
  }
#ifdef DEBUG
  if(!(std::accumulate(this->mask_.begin(),this->mask_.end(),0.0)>0)) {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: mask file {} left no data points after masking",
      "set_mask",
      MASK
    );
      exit(1);
  }
  spdlog::info(
    "\x1b[90m{}\x1b[0m: mask file {} left {} non-masked elements after masking",
    "set_mask",
    MASK,
    arma::accu(this->mask_)
  );
#endif
  this->mask_filename_ = MASK;
  this->is_mask_set_ = true;
}


void ima::RealData::set_inv_cov(std::string COV) {
#ifdef DEBUG
  if (!(like.Ndata>0)) {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} not set prior to this function call",
      "set_inv_cov", "like.Ndata");
    exit(1);
  }
  if (!(this->is_mask_set_)) {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} not set prior to this function call",
      "set_inv_cov", "mask");
    exit(1);
  }
#endif
  arma::Mat<double> table = ima::read_table(COV); // this reads cov!
  this->inv_cov_mask_.set_size(like.Ndata,like.Ndata);
  this->inv_cov_mask_.zeros();
  switch (table.n_cols) {
    case 3:
    {
      for (int i=0; i<static_cast<int>(table.n_rows); i++) {
        const int j = static_cast<int>(table(i,0));
        const int k = static_cast<int>(table(i,1));
        this->inv_cov_mask_(j,k) = table(i,2);
        if (j!=k) {
          // apply mask to off-diagonal covariance elements
          this->inv_cov_mask_(j,k) *= this->get_mask(j);
          this->inv_cov_mask_(j,k) *= this->get_mask(k);
          // m(i,j) = m(j,i)
          this->inv_cov_mask_(k,j) = this->inv_cov_mask_(j,k);
        }
      };
      break;
    }
    case 4:
    {
      for (int i=0; i<static_cast<int>(table.n_rows); i++) {
        const int j = static_cast<int>(table(i,0));
        const int k = static_cast<int>(table(i,1));
        this->inv_cov_mask_(j,k) = table(i,2) + table(i,3);
        this->inv_cov_mask_(k,j) = this->inv_cov_mask_(j,k);
        if (j!=k) {
          // apply mask to off-diagonal covariance elements
          this->inv_cov_mask_(j,k) *= this->get_mask(j);
          this->inv_cov_mask_(j,k) *= this->get_mask(k);
          // m(i,j) = m(j,i)
          this->inv_cov_mask_(k,j) = this->inv_cov_mask_(j,k);
        }
      };
      break;
    }
    case 10:
    {
      for (int i=0; i<static_cast<int>(table.n_rows); i++) {
        const int j = static_cast<int>(table(i,0));
        const int k = static_cast<int>(table(i,1));
        this->inv_cov_mask_(j,k) = table(i,8) + table(i,9);
        this->inv_cov_mask_(k,j) = this->inv_cov_mask_(j,k);
        this->inv_cov_mask_(k,j) = this->inv_cov_mask_(j,k);
        if (j!=k) {
          // apply mask to off-diagonal covariance elements
          this->inv_cov_mask_(j,k) *= this->get_mask(j);
          this->inv_cov_mask_(j,k) *= this->get_mask(k);
          // m(i,j) = m(j,i)
          this->inv_cov_mask_(k,j) = this->inv_cov_mask_(j,k);
        }
      }
      break;
    }
    default:
      spdlog::critical("{}: data format for covariance file = {} is invalid"
        "set_inv_cov", COV);
      exit(1);
  }
  this->inv_cov_mask_ = arma::inv(this->inv_cov_mask_);
  // apply mask again, to make sure numerical errors in matrix
  // inversion don't cause problems...
  // also, set diagonal elements corresponding to datavector elements
  // outside mask to zero, so that these elements don't contribute to chi2
  for (int i = 0; i < like.Ndata; i++) {
    this->inv_cov_mask_(i,i) *= this->get_mask(i)*this->get_mask(i);
    for (int j=0; j<i; j++) {
      this->inv_cov_mask_(i,j) *= this->get_mask(i)*this->get_mask(j);
      this->inv_cov_mask_(j,i) = this->inv_cov_mask_(i,j);
    }
  };
  this->cov_filename_ = COV;
  this->is_inv_cov_set_ = true;
}

arma::Col<double> ima::RealData::get_mask() const {
  return this->mask_;
}

int ima::RealData::get_mask(const int ci) const {
#ifdef DEBUG
  if (ci > like.Ndata || ci < 0) {
    spdlog::critical("\x1b[90m{}\x1b[0m: index i = {} is not valid (min = {}, max = {})",
      "get_mask", ci, 0.0, like.Ndata);
    exit(1);
  }
#endif
  return static_cast<int>(this->mask_(ci)+1e-13);
}

double ima::RealData::get_data(const int ci) const {
#ifdef DEBUG
  if (ci > like.Ndata || ci < 0) {
    spdlog::critical("\x1b[90m{}\x1b[0m: index i = {} is not valid (min = {}, max = {})",
      "get_data", ci, 0.0, like.Ndata);
    exit(1);
  }
#endif
    return this->data_(ci);
}

arma::Col<double> ima::RealData::get_data() const {
  return this->data_;
}

double ima::RealData::get_inv_cov(const int ci, const int cj) const {
#ifdef DEBUG
  if (ci > like.Ndata || ci < 0) {
    spdlog::critical("\x1b[90m{}\x1b[0m: index i = {} is not valid (min = {}, max = {})",
      "get_inv_cov", ci, 0.0, like.Ndata);
    exit(1);
  }
  if (cj > like.Ndata || cj < 0) {
    spdlog::critical("\x1b[90m{}\x1b[0m: index j = {} is not valid (min = {}, max = {})",
      "get_inv_cov", cj, 0.0, like.Ndata);
    exit(1);
  }
#endif
  return this->inv_cov_mask_(ci,cj);
}

double ima::RealData::get_chi2(std::vector<double> datavector) const {
#ifdef DEBUG
  if (!(this->is_data_set_)) {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} not set prior to this function call", "get_chi2",
      "data_vector");
    exit(1);
  }
  if (!(this->is_mask_set_)) {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} not set prior to this function call", "get_chi2",
      "mask");
    exit(1);
  }
  if (!(this->is_inv_cov_set_)) {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} not set prior to this function call", "get_chi2",
      "inv_cov");
    exit(1);
  }
#endif
  double chisqr = 0.0;
  for (int i=0; i<like.Ndata; i++) {
    if (this->get_mask(i)) {
      const double x = datavector[i]-this->get_data(i);
      for (int j=0; j<like.Ndata; j++) {
        if (this->get_mask(j)) {
          const double y = datavector[j]-this->get_data(j);
          chisqr += x*this->get_inv_cov(i,j)*y;
        }
      }
    }
  }
  if (chisqr<0.0) {
    spdlog::critical("\x1b[90m{}\x1b[0m: chi2 = {} (invalid)", "get_chi2", chisqr);
    exit(1);
  }
  return chisqr;
}

bool ima::RealData::is_mask_set() const {
  return this->is_mask_set_;
}

bool ima::RealData::is_data_set() const {
  return this->is_data_set_;
}

bool ima::RealData::is_inv_cov_set() const {
  return this->is_inv_cov_set_;
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Routines to test functions that were modified by Cocoa
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

arma::Mat<double> cpp_print_a() {
  const int Nchi = 1000;
  arma::Mat<double> r(Nchi, 2, arma::fill::zeros);
  const double chi_min = 10.;
  const double chi_max = 7000.;
  const double dlnchi = std::log(chi_max/chi_min)/(Nchi-1.);
  for(int i=0; i<Nchi; i++) {
    r(i,0) = chi_min*exp(dlnchi*i)/cosmology.coverH0;
    r(i,1) = a_chi(r(i,0));
  }
  return r;
}

arma::Mat<double> cpp_print_chi() {
  arma::Mat<double> r(Ntable.N_a, 2, arma::fill::zeros);
  const double da = ((1.0-limits.a_min)/(Ntable.N_a-1.0));
  for (int i=0; i<Ntable.N_a; i++) {
    r(i,0) = limits.a_min + i * da;
    r(i,1) = chi(r(i,0));
  }
  return r;
}

arma::Mat<double> cpp_print_dchi_da() {
  arma::Mat<double> r(Ntable.N_a, 2, arma::fill::zeros);
  const double da = ((1.0-limits.a_min)/(Ntable.N_a-1.0));
  for (int i=0; i<Ntable.N_a; i++) {
    r(i,0) = limits.a_min + i * da;
    r(i,1) = dchi_da(r(i,0));
  }
  return r;
}

arma::Mat<double> cpp_print_growth() {
  arma::Mat<double> r(Ntable.N_a, 2, arma::fill::zeros);
  const double da = ((1.0-limits.a_min)/(Ntable.N_a-1.0));
  for (int i=0; i<Ntable.N_a; i++) {
    r(i,0) = limits.a_min + i * da;
    r(i,1) = growfac(r(i,0));
  }
  return r;
}

arma::Mat<double> cpp_print_fgrowth() {
  arma::Mat<double> r(Ntable.N_a, 2, arma::fill::zeros);
  const double da = ((1.0-limits.a_min)/(Ntable.N_a-1.0));
  for (int i=0; i<Ntable.N_a; i++) {
    r(i,0) = limits.a_min + i * da;
    r(i,1) = f_growth(1.0/r(i,0)-1.0);
  }
  return r;
}

arma::Mat<double> cpp_print_hoverh0() {
  arma::Mat<double> r(Ntable.N_a, 2, arma::fill::zeros);
  const double da = ((1.0-limits.a_min)/(Ntable.N_a-1.));
  for (int i=0; i<Ntable.N_a; i++) {
    r(i,0) = limits.a_min + i * da;
    r(i,1) = hoverh0(r(i,0));
  }
  return r;
}

arma::Mat<double> cpp_print_ps() {
  const double da = ((1.0-limits.a_min)/(Ntable.N_a-1.));
  const double logkmin = std::log(limits.k_min_mpc*cosmology.coverH0);
  const double logkmax = std::log(limits.k_max_mpc_class*cosmology.coverH0);
  const double dk = (logkmax-logkmin)/(Ntable.N_k_nlin-1.0);
  arma::Mat<double> r(Ntable.N_a*Ntable.N_k_nlin, 4, arma::fill::zeros);
  for (int i=0; i<Ntable.N_a; i++) {
    const double a = limits.a_min + i * da;
    for (int j=0; j<Ntable.N_k_nlin; j++) {
      const double k = std::exp(logkmin + j * dk);
      r(i*Ntable.N_k_nlin+j,0) = a;
      r(i*Ntable.N_k_nlin+j,1) = k/cosmology.coverH0;
      r(i*Ntable.N_k_nlin+j,2) = p_lin(k, a);
      r(i*Ntable.N_k_nlin+j,3) = p_nonlin(k, a);
    }
  }
  return r;
}

arma::Mat<double> cpp_print_datavector() {
  std::vector<double> tmp = cpp_compute_data_vector();
  arma::Mat<double> r(tmp.size(), 2, arma::fill::zeros);
  for (int i=0; i<static_cast<int>(tmp.size()); i++) {
    r(i,0) = i;
    r(i,1) = tmp[i];
  }
  return r;
}

#ifdef PYBIND11
PYBIND11_MODULE(cosmolike_lsst_foutier_interface, m) {
    m.doc() = "CosmoLike Interface for LSST-FOURIER 3x2 Module";

    m.def("initial_setup", &cpp_initial_setup, "Def Setup");

    m.def("init_probes", &cpp_init_probes, "Init Probes", py::arg("possible_probes"));

    m.def("init_survey_parameters", &cpp_init_survey, "Init Survey", py::arg("surveyname"), py::arg("area"), py::arg("sigma_e"));

    m.def("init_cosmo_runmode", &cpp_init_cosmo_runmode,"Init Run Mode", py::arg("is_linear"));

    m.def("init_IA", &cpp_init_IA, "Init IA", py::arg("ia_model"));

    m.def("init_binning", &cpp_init_binning,"Init Binning", py::arg("Ncl"), py::arg("l_min"), py::arg("l_max"));

    m.def("init_lens_sample", &cpp_init_lens_sample,"Init Lens Sample", py::arg("multihisto_file"), py::arg("Ntomo"), py::arg("ggl_cut"));

    m.def("init_source_sample", &cpp_init_source_sample,"Init Source Sample", py::arg("multihisto_file"), py::arg("Ntomo"));

    m.def("init_size_data_vector", &cpp_init_size_data_vector, "Init Size Data Vector");

    m.def("init_linear_power_spectrum", &cpp_init_linear_power_spectrum, "Transfer Linear Matter Power Spectrum from Cobaya to Cosmolike",py::arg("log10k"), py::arg("z"), py::arg("lnP"));

    m.def("init_non_linear_power_spectrum", &cpp_init_non_linear_power_spectrum, "Transfer Matter Power Spectrum from Cobaya to Cosmolike", py::arg("log10k"), py::arg("z"), py::arg("lnP"));

    m.def("init_growth", &cpp_init_growth, "Transfer Growth Factor from Cobaya to Cosmolike", py::arg("z"), py::arg("G"));

    m.def("init_distances", &cpp_init_distances, "Transfer chi(z) from Cobaya to Cosmolike", py::arg("z"), py::arg("chi"));

    m.def("compute_chi2", &cpp_compute_chi2,"Get chi^2", py::arg("datavector"));

    m.def("compute_data_vector", &cpp_compute_data_vector,"Get data vector");

    m.def("set_nuisance_ia", &cpp_set_nuisance_ia_mpp,"Set Nuisance IA Parameters", py::arg("A1"), py::arg("A2"), py::arg("B_TA"));

    m.def("set_nuisance_bias", &cpp_set_nuisance_bias,"Set Nuisance Bias Parameters", py::arg("B1"), py::arg("B2"), py::arg("B_MAG"));

    m.def("set_nuisance_shear_calib", &cpp_set_nuisance_shear_calib,"Set Shear Calibration Parameters", py::arg("M"));

    m.def("set_nuisance_clustering_photoz", &cpp_set_nuisance_clustering_photoz,"Set Clustering Shear Photo-Z Parameters", py::arg("bias"));

    m.def("set_nuisance_shear_photoz", &cpp_set_nuisance_shear_photoz,"Set Shear Photo-Z Parameters", py::arg("bias"));

    m.def("set_cosmological_parameters", &cpp_set_cosmological_parameters,  "Set Cosmological Parameters", py::arg("omega_matter"), py::arg("hubble"), py::arg("is_cached"));

    m.def("init_data_real", &cpp_init_data_real,"Init cov, mask and data", py::arg("COV"), py::arg("MASK"), py::arg("DATA"));

    m.def("print_chi", &cpp_print_chi,"Print Distance");

    m.def("print_dchi_da", &cpp_print_dchi_da, "Print dchi/da");

    m.def("print_growth", &cpp_print_growth, "Print Growth");

    m.def("print_fgrowth", &cpp_print_fgrowth, "Print f Growth");

    m.def("print_hoverh0", &cpp_print_hoverh0, "Print hoverh0");

    m.def("print_ps", &cpp_print_ps, "Print Linear and Non-Linear Power Spectrum");

    m.def("print_a", &cpp_print_a,"Print Scale Factor");

    m.def("print_datavector", &cpp_print_datavector, "Print data vector");
}
#endif // PYBIND11

int main() {
  cpp_initial_setup();
  std::cout << "GOODBYE" << std::endl;
  exit(1);
}