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
#include "cosmolike/basics.h"
#include "cosmolike/bias.h"
#include "cosmolike/baryons.h"
#include "cosmolike/cosmo2D.h"
#include "cosmolike/cosmo3D.h"
#include "cosmolike/halo.h"
#include "cosmolike/radial_weights.h"
#include "cosmolike/recompute.h"
#include "cosmolike/pt_cfastpt.h"
#include "cosmolike/redshift_spline.h"
#include "cosmolike/structs.h"

using Vector = arma::Col<double>;
using STLVector = std::vector<double>;
using Matrix = arma::Mat<double>;
using Cube = arma::Cube<double>;

namespace cosmolike_interface
{

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// INIT FUNCTIONS
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_initial_setup()
{
  spdlog::cfg::load_env_levels();
  
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "initial_setup");

  // restart variables to 0 so error check can flag bad initialization
  tomo.shear_Nbin = 0;
  tomo.clustering_Nbin = 0;

  like.shear_shear = 0;
  like.shear_pos = 0;
  like.pos_pos = 0;

  like.lmin_kappacmb = 0;
  like.lmax_kappacmb = 0;

  like.is_cmb_bandpower = 0;
  like.lmin_bp = 0; 
  like.lmax_bp = 0;
  like.Nbp = 0;

  like.Ncl = 0;
  like.lmin = 0;
  like.lmax = 0;

  like.gk = 0;
  like.kk = 0;
  like.ks = 0;
  
  // bias
  gbias.b1_function = &b1_per_bin;

  // no priors
  like.clusterN = 0;
  like.clusterWL = 0;
  like.clusterCG = 0;
  like.clusterCC = 0;

  // reset bias
  for (int i = 0; i < MAX_SIZE_ARRAYS; i++)
  {
    gbias.b[i] = 0.0;
    gbias.b2[i] = 0.0;
    gbias.b_mag[i] = 0.0;
  }

  // reset IA
  for (int i = 0; i < MAX_SIZE_ARRAYS; i++)
  {
    nuisance.A_z[i] = 0.0;
    nuisance.A2_z[i] = 0.0;
    nuisance.b_ta_z[i] = 0.0;
  }

  // reset clustering photo-z stretch parameters
  for (int i = 0; i < MAX_SIZE_ARRAYS; i++)
  {
    nuisance.stretch_zphot_clustering[i] = 1.0;
  }

  // nonlimber ?
  // nonLimber: 0; Limber: 1;
  // Doesn't support IA=6 AND nonLimber together
  like.adopt_limber_gg = 0;
  like.adopt_limber_gammat = 1;

  like.high_def_integration = 1;

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "initial_setup");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_init_accuracy_boost(
    const double accuracy_boost, 
    const double sampling_boost,
    const int integration_accuracy
  )
{
  Ntable.N_a = static_cast<int>(ceil(Ntable.N_a*sampling_boost));
  
  Ntable.N_ell_TATT = static_cast<int>(ceil(Ntable.N_ell_TATT*sampling_boost));
  
  Ntable.N_k_lin = static_cast<int>(ceil(Ntable.N_k_lin*sampling_boost));
  
  Ntable.N_k_nlin = static_cast<int>(ceil(Ntable.N_k_nlin*sampling_boost));
  
  Ntable.N_ell = static_cast<int>(ceil(Ntable.N_ell*sampling_boost));
  
  Ntable.N_theta  = static_cast<int>(ceil(Ntable.N_theta*sampling_boost));

  Ntable.N_M = static_cast<int>(ceil(Ntable.N_M*sampling_boost));

  precision.low /= accuracy_boost;
  
  precision.medium /= accuracy_boost;
  
  precision.high /= accuracy_boost;
  
  precision.insane /= accuracy_boost; 

  like.high_def_integration = integration_accuracy;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_init_probes(std::string possible_probes)
{
  spdlog::info("\x1b[90m{}\x1b[0m: Begins", "init_probes");

  if (possible_probes.compare("xi") == 0)
  {
    like.shear_shear = 1;

    spdlog::info(
      "\x1b[90m{}\x1b[0m: {} = {} selected", 
      "init_probes", 
      "possible_probes", 
      "xi"
    );
  }
  else if (possible_probes.compare("wtheta") == 0)
  {
    like.pos_pos = 1;

    spdlog::info(
      "\x1b[90m{}\x1b[0m: {} = {} selected", 
      "init_probes", 
      "possible_probes", 
      "wtheta"
    );
  }
  else if (possible_probes.compare("gammat") == 0)
  {
    like.shear_pos = 1;

    spdlog::info(
      "\x1b[90m{}\x1b[0m: {} = {} selected", 
      "init_probes", 
      "possible_probes", 
      "gammat"
    );
  }
  else if (possible_probes.compare("2x2pt") == 0)
  {
    like.shear_pos = 1;
    like.pos_pos = 1;

    spdlog::info(
      "\x1b[90m{}\x1b[0m: {} = {} selected", 
      "init_probes", 
      "possible_probes", 
      "2x2pt"
    );
  }
  else if (possible_probes.compare("3x2pt") == 0)
  {
    like.shear_shear = 1;
    like.shear_pos = 1;
    like.pos_pos = 1;

    spdlog::info(
      "\x1b[90m{}\x1b[0m: {} = {} selected", 
      "init_probes", 
      "possible_probes", 
      "3x2pt"
    );
  }
  else if (possible_probes.compare("xi_ggl") == 0)
  {
    like.shear_shear = 1;
    like.shear_pos = 1;

    spdlog::info(
      "\x1b[90m{}\x1b[0m: {} = {} selected", 
      "init_probes", 
      "possible_probes", 
      "xi + ggl (2x2pt)"
    );
  }
  else if (possible_probes.compare("5x2pt") == 0)
  {
    like.shear_shear = 1;
    like.shear_pos = 1;
    like.pos_pos = 1;
    like.gk = 1;
    like.ks = 1;
    spdlog::info(
      "\x1b[90m{}\x1b[0m: {} = {} selected", 
      "init_probes", 
      "possible_probes", 
      "5x2pt"
    );
  }
  else if (possible_probes.compare("6x2pt") == 0)
  {
    like.shear_shear = 1;
    like.shear_pos = 1;
    like.pos_pos = 1;
    like.gk = 1;
    like.ks = 1;
    like.kk = 1;
    spdlog::info(
      "\x1b[90m{}\x1b[0m: {} = {} selected", 
      "init_probes", 
      "possible_probes", 
      "6x2pt"
    );
  }
  else if (possible_probes.compare("c3x2pt") == 0)
  {
    like.gk = 1;
    like.ks = 1;
    like.kk = 1;
    spdlog::info(
      "\x1b[90m{}\x1b[0m: {} = {} selected", 
      "init_probes", 
      "possible_probes", 
      "c3x2pt (gk + sk + kk)"
    );
  }
  else
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} = {} probe not supported",
      "init_probes", 
      "possible_probes", 
      possible_probes
    );
    exit(1);
  }

  spdlog::info("\x1b[90m{}\x1b[0m: Ends", "init_probes");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_init_survey(
    std::string surveyname, 
    double area, 
    double sigma_e
  )
{
  spdlog::info("\x1b[90m{}\x1b[0m: Begins", "init_survey");

  if (surveyname.size() > CHAR_MAX_SIZE - 1)
  {
    exit(1);
  }
  if (!(surveyname.size()>0))
  {
    spdlog::critical("{}: incompatible input", "init_survey");
    exit(1);
  }

  memcpy(survey.name, surveyname.c_str(), surveyname.size() + 1);
  
  survey.area = area;
  
  survey.sigma_e = sigma_e;

  spdlog::info("\x1b[90m{}\x1b[0m: Ends", "init_survey");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_init_cosmo_runmode(const bool is_linear)
{
  spdlog::info("\x1b[90m{}\x1b[0m: Begins", "init_cosmo_runmode");

  std::string mode = is_linear ? "linear" : "Halofit";
  
  const size_t size = mode.size();
  
  memcpy(pdeltaparams.runmode, mode.c_str(), size + 1);

  spdlog::info(
    "\x1b[90m{}\x1b[0m: {} = {} selected", 
    "init_cosmo_runmode", 
    "runmode", mode
  );
  
  spdlog::info("\x1b[90m{}\x1b[0m: Ends", "init_cosmo_runmode");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_init_cmb(
    const double lmin_kappa_cmb, 
    const double lmax_kappa_cmb, 
    const double fwhm, 
    std::string pixwin_file
  ) 
{
  spdlog::info("\x1b[90m{}\x1b[0m: Begins", "init_cmb");

  like.lmin_kappacmb = lmin_kappa_cmb;
  
  like.lmax_kappacmb = lmax_kappa_cmb;
  
  // fwhm = beam size in arcmin - cmb.fwhm = beam size in rad
  cmb.fwhm = fwhm*2.90888208665721580e-4;
  
  memcpy(cmb.pathHealpixWinFunc, pixwin_file.c_str(), pixwin_file.size()+1);
  
  spdlog::info("\x1b[90m{}\x1b[0m: Ends", "init_cmb");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_init_binning_cmb_bandpower(
    const int Nbp, 
    const double lmin, 
    const double lmax
  )
{
  spdlog::info("\x1b[90m{}\x1b[0m: Begins", "init_binning_cmb_bandpower");

  if (!(Nbp > 0))
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} = {} not supported", 
      "init_binning_cmb_bandpower", 
      "like.Nbp", 
      Nbp
    );
    exit(1);
  }

  spdlog::info(
    "\x1b[90m{}\x1b[0m: {} = {} selected.", 
    "init_binning_cmb_bandpower", 
    "Nbp", 
    Nbp
  );
  
  spdlog::info(
    "\x1b[90m{}\x1b[0m: {} = {} selected.", 
    "init_binning_cmb_bandpower", 
    "l_min", 
    lmin
  );

  spdlog::info(
    "\x1b[90m{}\x1b[0m: {} = {} selected.", 
    "init_binning_cmb_bandpower", 
    "l_max", 
    lmax
  );

  like.Nbp = Nbp;
  
  like.lmin_bp = lmin;
  
  like.lmax_bp = lmax;
  
  spdlog::info("\x1b[90m{}\x1b[0m: Ends", "init_binning_cmb_bandpower");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_init_binning_fourier(
    const int Ncl, 
    const double lmin, 
    const double lmax
  )
{
  spdlog::info("\x1b[90m{}\x1b[0m: Begins", "init_binning_fourier");

  if (!(Ncl > 0))
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} = {} not supported", 
      "init_binning_fourier", 
      "like.Ncl", Ncl
    );
    exit(1);
  }

  spdlog::info(
    "\x1b[90m{}\x1b[0m: {} = {} selected.", 
    "init_binning_fourier", 
    "Ncl", Ncl
  );

  spdlog::info(
    "\x1b[90m{}\x1b[0m: {} = {} selected.", 
    "init_binning_fourier", 
    "l_min", lmin
  );
  
  spdlog::info(
    "\x1b[90m{}\x1b[0m: {} = {} selected.", 
    "init_binning_fourier", 
    "l_max", lmax
  );

  like.Ncl = Ncl;
  
  like.lmin = lmin;
  
  like.lmax = lmax;
  
  const double logdl = (std::log(like.lmax) - std::log(like.lmin))/like.Ncl;
  
  like.ell = (double*) malloc(sizeof(double)*like.Ncl);
  
  for (int i = 0; i < like.Ncl; i++)
  {
    like.ell[i] = std::exp(std::log(like.lmin)+(i+0.5)*logdl);
  
    spdlog::debug(
      "\x1b[90m{}\x1b[0m: Bin {:d} - {} = {:.4e}, {} = {:.4e} and {} = {:.4e}",
      "init_binning_fourier", i, 
      "lmin", lmin, 
      "ell", like.ell[i], 
      "lmax", lmax
    );
  }

  spdlog::info("\x1b[90m{}\x1b[0m: Ends", "init_binning_fourier");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_init_binning(
    const int Ntheta, 
    const double theta_min_arcmin, 
    const double theta_max_arcmin
  )
{
  spdlog::info("\x1b[90m{}\x1b[0m: Begins", "init_binning");

  if (!(Ntheta > 0))
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} = {} not supported", 
      "init_binning",
      "like.Ntheta", Ntheta
    );
    exit(1);
  }

  spdlog::info(
    "\x1b[90m{}\x1b[0m: {} = {} selected.", 
    "init_binning", 
    "Ntheta", Ntheta
  );

  spdlog::info(
    "\x1b[90m{}\x1b[0m: {} = {} selected.", 
    "init_binning", 
    "theta_min_arcmin", theta_min_arcmin
  );

  spdlog::info(
    "\x1b[90m{}\x1b[0m: {} = {} selected.", 
    "init_binning", 
    "theta_max_arcmin",  theta_max_arcmin
  );

  like.Ntheta = Ntheta;
  
  like.vtmin = theta_min_arcmin * 2.90888208665721580e-4; // arcmin to rad conversion
  
  like.vtmax = theta_max_arcmin * 2.90888208665721580e-4; // arcmin to rad conversion
  
  const double logdt = (std::log(like.vtmax)-std::log(like.vtmin))/like.Ntheta;
  
  like.theta = (double*) calloc(like.Ntheta, sizeof(double));

  constexpr double x = 2./ 3.;

  for (int i = 0; i < like.Ntheta; i++)
  {
    const double thetamin = std::exp(log(like.vtmin) + (i + 0.0) * logdt);
    
    const double thetamax = std::exp(log(like.vtmin) + (i + 1.0) * logdt);
    
    like.theta[i] = x * (std::pow(thetamax, 3) - std::pow(thetamin, 3)) /
      (thetamax*thetamax - thetamin*thetamin);

    spdlog::debug(
      "\x1b[90m{}\x1b[0m: Bin {:d} - {} = {:.4e}, {} = {:.4e} and {} = {:.4e}",
      "init_binning", i, 
      "theta_min [rad]", thetamin, 
      "theta [rad]", like.theta[i], 
      "theta_max [rad]", thetamax
    );
  }

  spdlog::info("\x1b[90m{}\x1b[0m: Ends", "init_binning");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_init_lens_sample(
    std::string multihisto_file, 
    const int Ntomo
  )
{
  spdlog::info("\x1b[90m{}\x1b[0m: Begins", "init_lens_sample");

  if (tomo.shear_Nbin == 0)
  {
    spdlog::critical(
      "{}: {} not set prior to this function call", 
      "init_lens_sample", 
      "tomo.shear_Nbin"
    );
    exit(1);
  }
  if (multihisto_file.size()>CHAR_MAX_SIZE-1)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: insufficient pre-allocated char memory (max = {}) for"
      "the string: {}", 
      "init_lens_sample", 
      CHAR_MAX_SIZE-1, 
      multihisto_file
    );
    exit(1);
  }
  if (!(multihisto_file.size() > 0))
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: empty {} string not supported", 
      "init_lens_sample", 
      "multihisto_file"
    );
    exit(1);
  }
  if (!(Ntomo > 0) || Ntomo > MAX_SIZE_ARRAYS)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} = {} not supported (max = {})", 
      "init_lens_sample", 
      "Ntomo", Ntomo, 
      MAX_SIZE_ARRAYS
    );
    exit(1);
  }

  memcpy(
    redshift.clustering_REDSHIFT_FILE, 
    multihisto_file.c_str(), 
    multihisto_file.size()+1
  );

  redshift.clustering_photoz = 4;
  
  tomo.clustering_Nbin = Ntomo;
  
  tomo.clustering_Npowerspectra = tomo.clustering_Nbin;

  spdlog::info(
    "\x1b[90m{}\x1b[0m: {} = {} selected.", "init_lens_sample",
    "clustering_REDSHIFT_FILE", multihisto_file
  );

  spdlog::info(
    "\x1b[90m{}\x1b[0m: {} = {} selected.", "init_lens_sample",
    "clustering_Nbin", Ntomo
  );

  survey.ggl_overlap_cut = 0.0;
  
  spdlog::info(
    "\x1b[90m{}\x1b[0m: {} = {} selected.", "init_lens_sample",
    "survey.ggl_overlap_cut", survey.ggl_overlap_cut
  );

  pf_photoz(0.1, 0);
  {
    int n = 0;
    for (int i = 0; i < tomo.clustering_Nbin; i++)
    {
      for (int j = 0; j < tomo.shear_Nbin; j++)
      {
        n += test_zoverlap(i, j);
      }
    }
    tomo.ggl_Npowerspectra = n;

    spdlog::info(
      "\x1b[90m{}\x1b[0m: tomo.ggl_Npowerspectra = {}",
      "init_lens_sample", 
      tomo.ggl_Npowerspectra
    );
  }

  spdlog::info("\x1b[90m{}\x1b[0m: Ends", "init_lens_sample");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_init_source_sample(
    std::string multihisto_file, 
    const int Ntomo
  )
{
  spdlog::info("\x1b[90m{}\x1b[0m: Begins", "init_source_sample");

  if (multihisto_file.size() > CHAR_MAX_SIZE - 1)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: insufficient pre-allocated char memory (max = {}) for"
      "the string: {}", 
      "init_source_sample", 
      CHAR_MAX_SIZE-1, 
      multihisto_file
    );
    exit(1);
  }
  if (!(multihisto_file.size() > 0))
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: empty {} string not supported",
      "init_source_sample", 
      "multihisto_file"
    );
    exit(1);
  }
  if (!(Ntomo > 0) || Ntomo > MAX_SIZE_ARRAYS)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} = {} not supported (max = {})",
      "init_source_sample", 
      "Ntomo", 
      Ntomo, 
      MAX_SIZE_ARRAYS
    );
    exit(1);
  }

  // convert std::string to char*
  memcpy(
    redshift.shear_REDSHIFT_FILE, 
    multihisto_file.c_str(), 
    multihisto_file.size() + 1
  );

  redshift.shear_photoz = 4;
  tomo.shear_Nbin = Ntomo;
  tomo.shear_Npowerspectra = tomo.shear_Nbin * (tomo.shear_Nbin + 1) / 2;

  spdlog::info(
    "\x1b[90m{}\x1b[0m: tomo.shear_Npowerspectra = {}", 
    "init_source_sample", 
    tomo.shear_Npowerspectra
  );

  for (int i=0; i<tomo.shear_Nbin; i++)
  {
    nuisance.bias_zphot_shear[i] = 0.0;

    spdlog::info(
      "\x1b[90m{}\x1b[0m: bin {} - {} = {}.",
      "init_source_sample", 
      i, 
      "<z_s>", 
      zmean_source(i)
    );
  }

  spdlog::info(
    "\x1b[90m{}\x1b[0m: {} = {} selected.", 
    "init_source_sample",
    "shear_REDSHIFT_FILE", 
    multihisto_file
  );

  spdlog::info(
    "\x1b[90m{}\x1b[0m: {} = {} selected.", 
    "init_source_sample",
    "shear_Nbin", 
    Ntomo
  );

  spdlog::info("\x1b[90m{}\x1b[0m: Ends", "init_source_sample");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_init_cmb_bandpower(
    const int is_cmb_bandpower, 
    const int is_cmb_kkkk_cov_from_sim, 
    const double alpha
  )
{
  spdlog::info("\x1b[90m{}\x1b[0m: Begins", "init_cmb_bandpower");

  like.is_cmb_bandpower = is_cmb_bandpower;
  
  like.is_cmb_kkkk_cov_from_sim = is_cmb_kkkk_cov_from_sim;
  
  like.alpha_Hartlap_kkkk = alpha;

  spdlog::info("\x1b[90m{}\x1b[0m: Ends", "init_cmb_bandpower");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_init_cmb_bandpower_data(std::string BINMAT, std::string OFFSET)
{
  spdlog::info("\x1b[90m{}\x1b[0m: Begins", "init_cmb_bandpower_data");

  IPCMB::get_instance().set_cmb_binning_mat(BINMAT);
  
  IPCMB::get_instance().set_cmb_theory_offset(OFFSET);
  
  spdlog::info("\x1b[90m{}\x1b[0m: Ends", "init_cmb_bandpower_data");

  return;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_init_baryon_pca_scenarios(std::string scenarios)
{
  spdlog::info("\x1b[90m{}\x1b[0m: Begins", "init_baryon_pca_scenarios");

  BaryonScenario::get_instance().set_scenarios(scenarios);

  spdlog::info("\x1b[90m{}\x1b[0m: Ends", "init_baryon_pca_scenarios");

  return;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_init_baryons_contamination(
    const bool use_baryonic_sims_contamination,
    const std::string which_baryonic_sim
  )
{
  spdlog::info("\x1b[90m{}\x1b[0m: Begins", "init_baryons_contamination");

  spdlog::info(
    "\x1b[90m{}\x1b[0m: {} = {} selected",
    "init_baryons_contamination", 
    "use_baryonic_simulations", use_baryonic_sims_contamination
  );

  if (use_baryonic_sims_contamination)
  {
    // Determine simulation name and ID
    // simulation names only have _ as deliminator, e.g. owls_AGN, antilles
    // simulation IDs are indicated by -, e.g. antilles-1, antilles-386, BAHAMAS-3
    
    // To accommodate default behavior, 
    // owls_AGN_T80 = owls_AGN-1 // owls_AGN_T85 = owls_AGN-2
    // owls_AGN_T87 = owls_AGN-3 // BAHAMAS_T78 = BAHAMAS-1
    // BAHAMAS_T76 = BAHAMAS-2   // BAHAMAS_T80 = BAHAMAS-3
    
    int sim_id;
    std::string sim_name;
    
    if (which_baryonic_sim.rfind("owls_AGN") != std::string::npos)
    {
      sim_name = "owls_AGN";
      const std::string  sim_tag = which_baryonic_sim.substr(9);
      
      if ((sim_tag == "T80") || (sim_tag == "1")) 
      {
        sim_id = 1;
      }
      else if ((sim_tag == "T85") || (sim_tag == "2")) 
      {
        sim_id = 2;
      }
      else if ((sim_tag == "T87") || (sim_tag == "3")) 
      {
        sim_id = 3;
      }
      else
      {
        spdlog::critical(
          "\x1b[90m{}\x1b[0m: {} = {} not supported", 
          "cpp_init_baryons_contamination", 
          "which_baryonic_simulations_contamination", which_baryonic_sim
        );
        exit(1);
      }
    } 
    else if (which_baryonic_sim.rfind("BAHAMAS") != std::string::npos)
    {
      sim_name = "BAHAMAS";
      const std::string sim_tag = which_baryonic_sim.substr(8);
      
      if ((sim_tag == "T78") || (sim_tag == "1")) 
      {
        sim_id = 1;
      }
      else if ((sim_tag == "T76") || (sim_tag == "2")) 
      {
        sim_id = 2;
      }
      else if ((sim_tag == "T80") || (sim_tag == "3")) 
      {
        sim_id = 3;
      }
      else
      {
        spdlog::critical(
          "\x1b[90m{}\x1b[0m: {} = {} not supported", 
          "cpp_init_baryons_contamination", 
          "which_baryonic_simulations_contamination", which_baryonic_sim
        );
        exit(1);
      }
    } 
    else
    { // Second, check whether there's deliminator - for simulation ID
      if (which_baryonic_sim.rfind('-') != std::string::npos)
      { // Found - deliminator, parse simulation name and sim ID
        sim_name = which_baryonic_sim.substr(0, found1);
        sim_id = std::stoi(which_baryonic_sim.substr(found1 + 1));
      } 
      else
      { // Not found - deliminator, only parse simulation name, sim ID = 1
        
        sim_name = which_baryonic_sim;
        
        sim_id = 1;
      }
    }

    // Initialize baryon feedback scenario
    init_baryons(sim_name.c_str(), sim_id);

    spdlog::info(
      "\x1b[90m{}\x1b[0m: {} = {} selected (Sim {} ID.{})",
      "init_baryons_contamination", 
      "which_baryonic_simulations_contamination",
      which_baryonic_sim, 
      sim_name, 
      sim_id
    );
  }
  else
  {
    reset_bary_struct();
  }

  spdlog::info("\x1b[90m{}\x1b[0m: Ends", "init_baryons_contamination");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_init_size_6x2pt_data_vector()
{
  spdlog::info("\x1b[90m{}\x1b[0m: Begins", "init_size_data_vector");

  if (tomo.shear_Nbin == 0)
  {
    spdlog::debug(
      "{}: {} not set prior to this function call",
      "init_size_data_vector", 
      "tomo.shear_Nbin"
    );
  }
  if (tomo.clustering_Nbin == 0)
  {
    spdlog::debug(
      "{}: {} not set prior to this function call",
      "init_size_data_vector", 
      "tomo.clustering_Nbin"
    );
  }
  if (like.Ntheta == 0) 
  {
    spdlog::critical(
      "{}: {} not set prior to this function call",
      "init_size_data_vector", 
      "like.Ntheta"
    );
    exit(1);
  }
  if (like.is_cmb_bandpower == 0)
  {
    if (like.Ncl == 0)
    {
      spdlog::critical(
        "{}: {} not set prior to this function call", 
        "init_size_data_vector", 
        "like.Ncl"
      );
      exit(1);
    }
  }
  else if (like.is_cmb_bandpower == 1)
  {
    if (like.Nbp == 0)
    {
      spdlog::critical(
        "{}: {} not set prior to this function call", 
        "init_size_data_vector", "like.Nbp"
      );
      exit(1);
    }    
  }
  else
  {
    spdlog::critical(
        "{}: {} not set prior to this function call", 
        "init_size_data_vector", 
        "is_cmb_bandpower"
    );
  }

  like.Ndata = like.Ntheta*(2*tomo.shear_Npowerspectra + 
                              tomo.ggl_Npowerspectra + 
                              tomo.clustering_Npowerspectra + 
                              tomo.shear_Nbin + 
                              tomo.clustering_Nbin
                            );

  if (like.is_cmb_bandpower == 0)
  {
    like.Ndata += like.Ncl;
  }
  else
  {
    like.Ndata += like.Nbp;
  }
  
  spdlog::info("\x1b[90m{}\x1b[0m: {} = {} selected.", 
               "init_size_data_vector", 
               "Ndata", 
               like.Ndata
              );
  
  spdlog::info("\x1b[90m{}\x1b[0m: Ends", "init_size_data_vector");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_init_linear_power_spectrum(
    STLVector io_log10k,
    STLVector io_z, 
    STLVector io_lnP
  )
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_linear_power_spectrum");

  {
    bool debug_fail = false;
    if (io_z.size()*io_log10k.size() != io_lnP.size())
    {
      debug_fail = true;
    }
    else
    {
      if (io_z.size() == 0 || io_log10k.size() == 0)
      {
        debug_fail = true;
      }
    }
    if (debug_fail)
    {
      spdlog::critical(
        "\x1b[90m{}\x1b[0m: incompatible input w/ k.size = {}, z.size = {}, "
        "and lnP.size = {}", "init_linear_power_spectrum", io_log10k.size(),
        io_z.size(), io_lnP.size());
      exit(1);
    }

    if(io_z.size() < 5 || io_log10k.size() < 5)
    {
      spdlog::critical(
        "\x1b[90m{}\x1b[0m: bad input w/ k.size = {}, z.size = {}, "
        "and lnP.size = {}", "init_linear_power_spectrum", io_log10k.size(),
        io_z.size(), io_lnP.size());
      exit(1);
    }
  }

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

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_linear_power_spectrum");

  return;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_init_non_linear_power_spectrum(
    STLVector io_log10k,
    STLVector io_z, 
    STLVector io_lnP
  )
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_non_linear_power_spectrum");

  {
    bool debug_fail = false;
    
    if (io_z.size()*io_log10k.size() != io_lnP.size())
    {
      debug_fail = true;
    }
    else
    {
      if (io_z.size() == 0)
      {
        debug_fail = true;
      }
    }

    if (debug_fail)
    {
      spdlog::critical(
        "\x1b[90m{}\x1b[0m: incompatible input w/ k.size = {}, z.size = {}, "
        "and lnP.size = {}", 
        "init_non_linear_power_spectrum", 
        io_log10k.size(),
        io_z.size(), 
        io_lnP.size()
      );
      exit(1);
    }

    if (io_z.size() < 5 || io_log10k.size() < 5)
    {
      spdlog::critical(
        "\x1b[90m{}\x1b[0m: bad input w/ k.size = {}, z.size = {}, "
        "and lnP.size = {}", 
        "init_non_linear_power_spectrum", 
        io_log10k.size(),
        io_z.size(), 
        io_lnP.size()
      );
      exit(1);
    }
  }

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

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_non_linear_power_spectrum");

  return;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// Growth: D = G * a
void cpp_init_growth(STLVector io_z, STLVector io_G)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_growth");

  {
    bool debug_fail = false;
    
    if (io_z.size() != io_G.size())
    {
      debug_fail = true;
    }
    else
    {
      if (io_z.size() == 0)
      {
        debug_fail = true;
      }
    }
    
    if (debug_fail)
    {
      spdlog::critical(
        "\x1b[90m{}\x1b[0m: incompatible input w/ z.size = {} and G.size = {}",
        "init_growth", 
        io_z.size(), 
        io_G.size()
      );
      exit(1);
    }
  }

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

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_growth");

  return;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_init_distances(STLVector io_z, STLVector io_chi)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "init_distances");

  {
    bool debug_fail = false;
    
    if (io_z.size() != io_chi.size())
    {
      debug_fail = true;
    }
    else
    {
      if (io_z.size() == 0)
      {
        debug_fail = true;
      }
    }
    
    if (debug_fail)
    {
      spdlog::critical(
        "\x1b[90m{}\x1b[0m: incompatible input w/ z.size = {} and G.size = {}",
        "init_distances",
        io_z.size(),
        io_chi.size()
      );
      exit(1);
    }
  }

  int nz = static_cast<int>(io_z.size());
  double* vz = io_z.data();
  double* vchi = io_chi.data();
  setup_chi(&nz, &vz, &vchi, 1);

  // force initialization - imp to avoid seg fault when openmp is on
  const double io_a = 1.0;
  chi(io_a);

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "init_distances");

  return;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// SET FUNCTIONS
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_set_nuisance_shear_calib(STLVector M)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "set_nuisance_shear_calib");

  if (tomo.shear_Nbin == 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} = 0 is invalid", 
      "set_nuisance_shear_calib",
      "shear_Nbin"
    );
    exit(1);
  }

  if (tomo.shear_Nbin != static_cast<int>(M.size()))
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: incompatible input w/ size = {} (!= {})",
      "set_nuisance_shear_calib", 
      M.size(), 
      tomo.shear_Nbin
    );
    exit(1);
  }

  for (int i=0; i<tomo.shear_Nbin; i++)
  {
    nuisance.shear_calibration_m[i] = M[i];
  }

   spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "set_nuisance_shear_calib");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_set_nuisance_shear_photoz(STLVector SP)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "set_nuisance_shear_photoz");

  if (tomo.shear_Nbin == 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} = 0 is invalid",
      "set_nuisance_shear_photoz",
      "shear_Nbin"
    );
    exit(1);
  }

  if (tomo.shear_Nbin != static_cast<int>(SP.size()))
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: incompatible input w/ size = {} (!= {})",
      "set_nuisance_shear_photoz",
      SP.size(),
      tomo.shear_Nbin
    );
    exit(1);
  }

  for (int i=0; i<tomo.shear_Nbin; i++)
  {
    nuisance.bias_zphot_shear[i] = SP[i];
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "set_nuisance_shear_photoz");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_set_nuisance_clustering_photoz(STLVector CP)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "set_nuisance_clustering_photoz");

  if (tomo.clustering_Nbin == 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} = 0 is invalid",
      "set_nuisance_clustering_photoz",
      "clustering_Nbin"
    );
    exit(1);
  }

  if (tomo.clustering_Nbin != static_cast<int>(CP.size()))
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: incompatible input w/ size = {} (!= {})",
      "set_nuisance_clustering_photoz",
      CP.size(),
      tomo.clustering_Nbin
    );
    exit(1);
  }

  for (int i=0; i<tomo.clustering_Nbin; i++)
  {
    nuisance.bias_zphot_clustering[i] = CP[i];
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "set_nuisance_clustering_photoz");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_set_nuisance_clustering_photoz_stretch(STLVector CP)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "set_nuisance_clustering_photoz_stretch");

  if (tomo.clustering_Nbin == 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} = 0 is invalid",
      "set_nuisance_clustering_photoz_stretch",
      "clustering_Nbin"
    );
    exit(1);
  }

  if (tomo.clustering_Nbin != static_cast<int>(CP.size()))
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: incompatible input w/ size = {} (!= {})",
      "set_nuisance_clustering_photoz_stretch",
      CP.size(),
      tomo.clustering_Nbin
    );
    exit(1);
  }

  for (int i=0; i<tomo.clustering_Nbin; i++)
  {
    nuisance.stretch_zphot_clustering[i] = CP[i];
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "set_nuisance_clustering_photoz_stretch");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_set_nuisance_linear_bias(STLVector B1)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "set_nuisance_linear_bias");

  if (tomo.clustering_Nbin == 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} = 0 is invalid",
      "set_nuisance_linear_bias", "clustering_Nbin"
    );
    exit(1);
  }

  if (tomo.clustering_Nbin != static_cast<int>(B1.size()))
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: incompatible input w/ size = {} (!= {})",
      "set_nuisance_linear_bias", 
      B1.size(), 
      tomo.clustering_Nbin
    );
    exit(1);
  }

  for (int i=0; i<tomo.clustering_Nbin; i++)
  {
    gbias.b[i] = B1[i];
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "set_nuisance_linear_bias");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_set_nuisance_nonlinear_bias(
    STLVector B1,
    STLVector B2
  )
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "set_nuisance_nonlinear_bias");

  if (tomo.clustering_Nbin == 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} = 0 is invalid",
      "set_nuisance_nonlinear_bias", 
      "clustering_Nbin"
    );
    exit(1);
  }

  if ( tomo.clustering_Nbin != static_cast<int>(B1.size()) ||
       tomo.clustering_Nbin != static_cast<int>(B2.size()))
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: incompatible input w/ sizes = {} and {} (!= {})",
      "set_nuisance_nonlinear_bias", 
      B1.size(), 
      B2.size(), 
      tomo.clustering_Nbin
    );
    exit(1);
  }

  constexpr double tmp = -4./7.;
  
  for (int i=0; i<tomo.clustering_Nbin; i++)
  {
    gbias.b2[i] = B2[i];
    gbias.bs2[i] = almost_equal(B2[i], 0.) ? 0 : tmp*(B1[i]-1.0);
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "set_nuisance_nonlinear_bias");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_set_nuisance_magnification_bias(STLVector B_MAG)
{
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "set_nuisance_magnification_bias");

  if (tomo.clustering_Nbin == 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} = 0 is invalid",
      "set_nuisance_magnification_bias",
      "clustering_Nbin"
    );
    exit(1);
  }

  if (tomo.clustering_Nbin != static_cast<int>(B_MAG.size()))
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: incompatible input w/ size = {} (!= {})",
      "set_nuisance_magnification_bias", 
      B_MAG.size(), 
      tomo.clustering_Nbin
    );
    exit(1);
  }

  for (int i=0; i<tomo.clustering_Nbin; i++)
  {
    gbias.b_mag[i] = B_MAG[i];
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "set_nuisance_magnification_bias");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_set_nuisance_bias(STLVector B1, STLVector B2, STLVector B_MAG)
{
  cpp_set_nuisance_linear_bias(B1);
  
  cpp_set_nuisance_nonlinear_bias(B1, B2);
  
  cpp_set_nuisance_magnification_bias(B_MAG);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_set_pm(STLVector pm)
{  
  PointMass::get_instance().set_pm_vector(pm); 
  return;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// GET FUNCTIONS
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

Matrix cpp_get_cov_masked()
{
  return IP::get_instance().get_cov_masked();
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

Matrix cpp_get_cov_masked_sqzd()
{
  return IP::get_instance().get_cov_masked_sqzd();
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int cpp_get_mask(const int i)
{
  return IP::get_instance().get_mask(i);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int cpp_get_ndata()
{
  return IP::get_instance().get_ndata();
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int cpp_get_ndata_sqzd()
{
  return IP::get_instance().get_ndata_sqzd();
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int cpp_get_index_sqzd(const int i)
{
  return IP::get_instance().get_index_sqzd(i);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

STLVector expand_ndata_from_masked_sqzd(STLVector sqzd)
{
  Vector tmp = IP::get_instance().expand_ndata_from_masked_sqzd(Vector(sqzd));

  // The conversion between STL vector and python np array is cleaner
  // arma:Col is cast to 2D np array with 1 column (not as nice!)
  
  STLVector res(tmp.n_elem, 0.0);
  
  for(int i=0; i<static_cast<int>(tmp.n_elem); i++)
  {
    res[i] = tmp(i);
  }

  return res;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int cpp_get_baryon_pca_nscenarios()
{
  return BaryonScenario::get_instance().nscenarios();
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

std::string cpp_get_baryon_pca_scenario_name(const int i)
{
  return BaryonScenario::get_instance().get_scenario(i);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// COMPUTE FUNCTIONS
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double cpp_compute_chi2(STLVector datavector)
{
  return IP::get_instance().get_chi2(datavector);
}

double cpp_compute_pm(const int zl, const int zs, const double theta)
{
  return PointMass::get_instance().get_pm(zl,zs,theta);
}

double cpp_compute_baryon_ratio(const double log10k, const double a)
{
  const double KNL = pow(10.0,log10k)*cosmology.coverH0;
  return PkRatio_baryons(KNL, a);
}

// ordering = (1,2,3,4,5,6) => Cosmic Shear, ggl, gg, gk, sk, kk
// ordering = (2,1,3,4,5,6) => ggl, Cosmic Shear, gg, gk, sk, kk ...
STLVector cpp_compute_data_vector_6x2pt_masked(arma::Col<uword>::fixed<6> ordering)
{ 
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "compute_data_vector_masked");

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

  if (like.is_cmb_bandpower == 1)
  {
    if (!(like.Nbp > 0))
    {
      spdlog::critical(
        "\x1b[90m{}\x1b[0m: {} not set prior to this function call",
        "compute_data_vector_masked", 
        "like.Ncl"
      );
      exit(1);
    }

    if (like.lmin_bp < 0 or like.lmax_bp <= 0)
    { // check ell range
      spdlog::critical(
        "\x1b[90m{}\x1b[0m: {} and {} are invalid",
        "compute_data_vector_masked", 
        "like.lmin_bp", 
        "like.lmax_bp"
      );
      exit(1);
    }

    if (!IP::get_instance().is_cmb_binmat_set())
    { // check binning matrix and CMB lensing band power offset
      spdlog::critical(
        "\x1b[90m{}\x1b[0m: {} not set prior to this function call",
        "compute_data_vector_masked", 
        "cmb_binning_matrix_with_correction"
      );
      exit(1);
    }

    if (!IP::get_instance().is_cmb_offset_set())
    {
      spdlog::critical(
        "\x1b[90m{}\x1b[0m: {} not set prior to this function call",
        "compute_data_vector_masked", 
        "cmb_theory_offset"
      );
      exit(1);
    }
  }
  else if (like.is_cmb_bandpower == 0) 
  {
    if (!(like.Ncl > 0))
    {
      spdlog::critical(
        "\x1b[90m{}\x1b[0m: {} not set prior to this function call",
        "compute_data_vector_masked", 
        "like.Ncl"
      );
      exit(1);
    }

    if (like.lmin < 0 or like.lmax <= 0)
    {
      spdlog::critical(
        "\x1b[90m{}\x1b[0m: {} and {} are invalid",
        "compute_data_vector_masked", 
        "like.lmin", 
        "like.lmax"
      );
      exit(1);
    }
  }
  else
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} not set prior to this function call",
      "compute_data_vector_masked", 
      "like.is_cmb_bandpower"
    );
    exit(1);
  }

  if (!IP::get_instance().is_mask_set())
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} not set prior to this function call",
      "compute_data_vector_masked", 
      "mask"
    );
    exit(1);
  }

  if (like.lmin_kappacmb <= 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} not set prior to this function call",
      "compute_data_vector_masked", 
      "like.lmin_kappacmb"
    );
    exit(1);
  }
  
  if (like.lmax_kappacmb <= 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} not set prior to this function call",
      "compute_data_vector_masked", 
      "like.lmax_kappacmb"
    );
    exit(1);
  }
  
  if (cmb.fwhm <= 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} not set prior to this function call",
      "compute_data_vector_masked", 
      "cmb.fwhm"
    );
    exit(1);
  }

  arma::Col<uword>::fixed<6> indices = arma::sort_index(ordering);

  arma::Col<uword>::fixed<6> sizes =
    {
      2*like.Ntheta*tomo.shear_Npowerspectra,
      like.Ntheta*tomo.ggl_Npowerspectra,
      like.Ntheta*tomo.clustering_Npowerspectra,
      like.Ntheta*tomo.clustering_Nbin,
      like.Ntheta*tomo.shear_Nbin,
      like.is_cmb_bandpower  == 1 ? like.Nbp : like.Ncl 
    }

  if(like.Ndata != static_cast<int>(arma::accu(sizes)))
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: logical error, internal inconsistency on the size of data vector",
      "cpp_compute_data_vector_6x2pt_masked"
    );
    exit(1);
  }

  arma::Col<uword>::fixed<6> start = {0,0,0,0,0,0};

  for(int i=0; i<6; i++)
  {
    for(int j=0; j<indices(i); j++)
    {
      start(i) += sizes(indices(j));
    }
  }
  
  std::vector<double> data_vector(like.Ndata, 0.0);
  
  if (like.shear_shear == 1)
  {
    for (int nz=0; nz<tomo.shear_Npowerspectra; nz++)
    {
      const int z1 = Z1(nz);
      const int z2 = Z2(nz);
  
      for (int i=0; i<like.Ntheta; i++)
      {
        int index = start(0) + like.Ntheta*nz + i;
        
        if (cpp_get_mask(index))
        {
          data_vector[index] = xi_pm_tomo(1, i, z1, z2, 1)*
            (1.0 + nuisance.shear_calibration_m[z1])*
            (1.0 + nuisance.shear_calibration_m[z2]);
        }
        
        index += like.Ntheta*tomo.shear_Npowerspectra;
        
        if (cpp_get_mask(index))
        {
          data_vector[index] = xi_pm_tomo(-1, i, z1, z2, 1)*
            (1.0 + nuisance.shear_calibration_m[z1])*
            (1.0 + nuisance.shear_calibration_m[z2]);
        }
      }
    }
  }

  if (like.shear_pos == 1)
  {
    for (int nz=0; nz<tomo.ggl_Npowerspectra; nz++)
    {
      const int zl = ZL(nz);
      const int zs = ZS(nz);
      
      for (int i=0; i<like.Ntheta; i++)
      {
        const int index = start(1) + like.Ntheta*nz + i;
        
        if (cpp_get_mask(index))
        {
          data_vector[index] = 
            ( w_gammat_tomo(i, zl, zs, like.adopt_limber_gammat) + 
              cpp_compute_pm(zl, zs, like.theta[i])
            )*(1.0+nuisance.shear_calibration_m[zs]);
        }
      }
    }
  }

  if (like.pos_pos == 1)
  {
    for (int nz=0; nz<tomo.clustering_Npowerspectra; nz++)
    {
      for (int i=0; i<like.Ntheta; i++)
      {
        const int index = start(2) + like.Ntheta*nz + i;
        if (cpp_get_mask(index))
        {
          data_vector[index] = w_gg_tomo(i, nz, nz, like.adopt_limber_gg);
        }
      }
    }
  }

  if (like.gk == 1)
  {
    for (int nz=0; nz<tomo.clustering_Nbin; nz++)
    {
      for (int i=0; i<like.Ntheta; i++) 
      {
        const int index = start(3) + like.Ntheta*nz + i;
        if (cpp_get_mask(index))
        {
          data_vector[index] = w_gk_tomo(i, nz, 1);
        }
      }
    }
  }

  if (like.ks == 1) 
  {
    for (int nz=0; nz<tomo.shear_Nbin; nz++)
    {
      for (int i=0; i<like.Ntheta; i++)
      {
        const int index = start(4) + like.Ntheta*nz + i; 
        if (cpp_get_mask(index))
        {
          data_vector[index] = 
            w_ks_tomo(i, nz, 1)*(1.0 + nuisance.shear_calibration_m[nz]);
        }
      }
    }
  }
  
  if (like.kk == 1)
  {
    if (like.is_cmb_bandpower == 0)
    {
      for (int i=0; i<like.Ncl; i++)  
      {
        const int index = start(5) + i; 
        if (cpp_get_mask(index))
        {
          const double l = like.ell[i];
          if (l <= limits.LMIN_tab)
          {
            data_vector[index] = C_kk_limber_nointerp(l, 0, 0);
          }
          else
          {
            data_vector[index] = C_kk_limber(l);
          }
        }
      }
    }
    else
    {
      for (int L=like.lmin_bp; L<like.lmax_bp + 1; L++)
      {
        const double Ckk = (L <= limits.LMIN_tab) ? 
          C_kk_limber_nointerp((double) L, 0, 0) : 
          C_kk_limber((double) L);

        const int i = L - like.lmin_bp;

        for (int j=0; j<like.Nbp; j++) // Loop through bandpower bins
        { 
          const int index = start(5) + j; 
          
          if (cpp_get_mask(index))
          {
            double x = 
              IP::get_instance().get_binning_matrix_with_correction(j, i);
        
            data_vector[index] += (Ckk * x);
          }
        }
      }

      for (int j=0; j<like.Nbp; j++) // offset due to marginalizing over primary CMB
      { 
        const int index = start(5) + j;
        if (cpp_get_mask(index))
        {
          data_vector[index] -= IP::get_instance().get_cmb_theory_offset(j);
        }
      }
    }
  }
  
  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "compute_data_vector_masked");

  return data_vector;
}

// ordering = (1,2,3) => Cosmic Shear, ggl, gg
// ordering = (2,1,3) => ggl, Cosmic Shear, gg  ...
STLVector cpp_compute_data_vector_3x2pt_masked(arma::Col<uword>::fixed<3> ordering)
{ 
  spdlog::debug("\x1b[90m{}\x1b[0m: Begins", "compute_data_vector_masked");

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

  if (!IP::get_instance().is_mask_set())
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} not set prior to this function call",
      "compute_data_vector_masked", 
      "mask"
    );
    exit(1);
  }

  arma::Col<uword>::fixed<3> indices = arma::sort_index(ordering);

  arma::Col<uword>::fixed<3> sizes =
    {
      2*like.Ntheta*tomo.shear_Npowerspectra,
      like.Ntheta*tomo.ggl_Npowerspectra,
      like.Ntheta*tomo.clustering_Npowerspectra
    }

  if(like.Ndata != static_cast<int>(arma::accu(sizes)))
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: logical error, internal inconsistency on the size of data vector",
      "compute_data_vector_3x2pt_masked"
    );
    exit(1);
  }

  arma::Col<uword>::fixed<3> start = {0,0,0};

  for(int i=0; i<3; i++)
  {
    for(int j=0; j<indices(i); j++)
    {
      start(i) += sizes(indices(j));
    }
  }
  
  std::vector<double> data_vector(like.Ndata, 0.0);
  
  if (like.shear_shear == 1)
  {
    for (int nz=0; nz<tomo.shear_Npowerspectra; nz++)
    {
      const int z1 = Z1(nz);
      const int z2 = Z2(nz);
  
      for (int i=0; i<like.Ntheta; i++)
      {
        int index = start(0) + like.Ntheta*nz + i;
        
        if (cpp_get_mask(index))
        {
          data_vector[index] = xi_pm_tomo(1, i, z1, z2, 1)*
            (1.0 + nuisance.shear_calibration_m[z1])*
            (1.0 + nuisance.shear_calibration_m[z2]);
        }
        
        index += like.Ntheta*tomo.shear_Npowerspectra;
        
        if (cpp_get_mask(index))
        {
          data_vector[index] = xi_pm_tomo(-1, i, z1, z2, 1)*
            (1.0 + nuisance.shear_calibration_m[z1])*
            (1.0 + nuisance.shear_calibration_m[z2]);
        }
      }
    }
  }

  if (like.shear_pos == 1)
  {
    for (int nz=0; nz<tomo.ggl_Npowerspectra; nz++)
    {
      const int zl = ZL(nz);
      const int zs = ZS(nz);
      
      for (int i=0; i<like.Ntheta; i++)
      {
        const int index = start(1) + like.Ntheta*nz + i;
        
        if (cpp_get_mask(index))
        {
          data_vector[index] = 
            ( w_gammat_tomo(i, zl, zs, like.adopt_limber_gammat) + 
              cpp_compute_pm(zl, zs, like.theta[i])
            )*(1.0+nuisance.shear_calibration_m[zs]);
        }
      }
    }
  }

  if (like.pos_pos == 1)
  {
    for (int nz=0; nz<tomo.clustering_Npowerspectra; nz++)
    {
      for (int i=0; i<like.Ntheta; i++)
      {
        const int index = start(2) + like.Ntheta*nz + i;
        if (cpp_get_mask(index))
        {
          data_vector[index] = w_gg_tomo(i, nz, nz, like.adopt_limber_gg);
        }
      }
    }
  }

  spdlog::debug("\x1b[90m{}\x1b[0m: Ends", "compute_data_vector_masked");

  return data_vector;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// RESET FUNCTIONS
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void cpp_reset_baryionic_struct()
{
  reset_bary_struct();
  return;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Class IP MEMBER FUNCTIONS
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

bool IP::is_mask_set() const
{
  return this->is_mask_set_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

bool IP::is_data_set() const
{
  return this->is_data_set_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

bool IP::is_inv_cov_set() const
{
  return this->is_inv_cov_set_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void IP::set_data(std::string DATA)
{
  if (!(this->is_mask_set_))
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} not set prior to this function call", 
      "set_data",
      "mask"
    );
    exit(1);
  }

  this->data_masked_.set_size(this->ndata_);
  
  this->data_filename_ = DATA;
  
  this->data_masked_sqzd_.set_size(this->ndata_masked_);

  Matrix table = read_table(DATA);

  for(int i=0; i<like.Ndata; i++)
  {
    this->data_masked_(i) = table(i,1);
    
    this->data_masked_(i) *= this->get_mask(i);

    if(this->get_mask(i) == 1)
    {
      if(this->get_index_sqzd(i) < 0)
      {
        spdlog::critical(
          "\x1b[90m{}\x1b[0m: logical error, internal"
          " inconsistent mask operation", 
          "set_data"
        );
        exit(1);
      }

      this->data_masked_sqzd_(this->get_index_sqzd(i)) = this->data_masked_(i);
    }
  }

  this->is_data_set_ = true;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void IP::set_mask(std::string MASK, arma::Col<uword>::fixed<6> ordering)
{
  if (!(like.Ndata>0))
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} not set prior to this function call",
      "set_mask", 
      "like.Ndata"
    );
    exit(1);
  }

  this->ndata_ = like.Ndata;

  this->mask_.set_size(this->ndata_);

  Matrix table = read_table(MASK);
  
  for (int i=0; i<this->ndata_; i++)
  {
    this->mask_(i) = static_cast<int>(table(i,1)+1e-13);
    
    if (!(this->mask_(i) == 0 || this->mask_(i) == 1))
    {
      spdlog::critical(
        "\x1b[90m{}\x1b[0m: inconsistent mask", 
        "set_mask"
      );
      exit(1);
    }
  }

  arma::Col<uword>::fixed<6> indices = arma::sort_index(ordering);

  arma::Col<uword>::fixed<6> sizes =
    {
      2*like.Ntheta*tomo.shear_Npowerspectra,
      like.Ntheta*tomo.ggl_Npowerspectra,
      like.Ntheta*tomo.clustering_Npowerspectra,
      like.Ntheta*tomo.clustering_Nbin,
      like.Ntheta*tomo.shear_Nbin,
      like.is_cmb_bandpower  == 1 ? like.Nbp : like.Ncl 
    }

  arma::Col<uword>::fixed<6> start = {0,0,0,0,0,0};

  for(int i=0; i<6; i++)
  {
    for(int j=0; j<indices(i); j++)
    {
      start(i) += sizes(indices(j));
    }
  }

  if (like.shear_shear == 0)
  {
    const int N = start(0);
    const int M = N + sizes(0);
    for (int i=N; i<M; i++)
    {
      this->mask_(i) = 0;
    }
  }

  if (like.shear_pos == 0)
  {
    const int N = start(1);
    const int M = N + sizes(1);
    for (int i=N; i<M; i++)
    {
      this->mask_(i) = 0;
    }
  }

  if (like.pos_pos == 0)
  {
    const int N = start(2);
    const int M = N + sizes(2);
    for (int i=N; i<M; i++)
    {
      this->mask_(i) = 0;
    }
  }

  if (like.gk == 0) 
  {
    const int N = start(3);
    const int M = N + sizes(3);;
    for (int i=N; i<M; i++) 
    {
      this->mask_(i) = 0.0;
    }
  }

  if (like.ks == 0) 
  {
    const int N = start(4);
    const int M = N + sizes(4);
    for (int i=N; i<M; i++) 
    {
      this->mask_(i) = 0.0;
    }
  }

  if (like.kk == 0) 
  {
    const int N = start(5);
    const int M = N + sizes(5);
    
    for (int i=N; i<M; i++) 
    {
      this->mask_(i) = 0.0;
    }
  }

  this->mask_filename_ = MASK;
  
  this->ndata_masked_ = arma::accu(this->mask_);
  
  if(!(this->ndata_masked_>0))
  {
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
    this->ndata_masked_
  );

  this->index_sqzd_.set_size(this->ndata_);
  {
    double j=0;
    
    for(int i=0; i<this->ndata_; i++)
    {
      if(this->get_mask(i) > 0)
      {
        this->index_sqzd_(i) = j;
        j++;
      }
      else
      {
        this->index_sqzd_(i) = -1;
      }
    }

    if(j != this->ndata_masked_)
    {
      spdlog::critical(
        "\x1b[90m{}\x1b[0m: logical error, internal inconsistent mask operation",
        "set_mask"
      );
      exit(1);
    }
  }

  this->is_mask_set_ = true;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void IP::set_inv_cov(std::string COV)
{
  if (!(this->is_mask_set_))
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} not set prior to this function call",
      "set_inv_cov",
      "mask"
    );
    exit(1);
  }

  Matrix table = read_table(COV); // this reads cov!

  this->cov_masked_.set_size(this->ndata_, this->ndata_);
  this->cov_masked_.zeros();

  this->inv_cov_masked_.set_size(this->ndata_, this->ndata_);
  this->inv_cov_masked_.zeros();

  switch (table.n_cols)
  {
    case 3:
    {
      for (int i=0; i<static_cast<int>(table.n_rows); i++)
      {
        const int j = static_cast<int>(table(i,0));
        const int k = static_cast<int>(table(i,1));

        this->cov_masked_(j,k) = table(i,2);
        this->inv_cov_masked_(j,k) = table(i,2);

        if (j!=k)
        {
          // apply mask to off-diagonal covariance elements
          this->cov_masked_(j,k) *= this->get_mask(j);
          this->cov_masked_(j,k) *= this->get_mask(k);

          this->inv_cov_masked_(j,k) *= this->get_mask(j);
          this->inv_cov_masked_(j,k) *= this->get_mask(k);

          // m(i,j) = m(j,i)
          this->cov_masked_(k,j) = this->cov_masked_(j,k);
          this->inv_cov_masked_(k,j) = this->inv_cov_masked_(j,k);
        }
      };
      break;
    }
    case 4:
    {
      for (int i=0; i<static_cast<int>(table.n_rows); i++)
      {
        const int j = static_cast<int>(table(i,0));
        const int k = static_cast<int>(table(i,1));

        this->cov_masked_(j,k) = table(i,2) + table(i,3);
        this->inv_cov_masked_(j,k) = table(i,2) + table(i,3);

        if (j!=k)
        {
          // apply mask to off-diagonal covariance elements
          this->cov_masked_(j,k) *= this->get_mask(j);
          this->cov_masked_(j,k) *= this->get_mask(k);

          this->inv_cov_masked_(j,k) *= this->get_mask(j);
          this->inv_cov_masked_(j,k) *= this->get_mask(k);

          // m(i,j) = m(j,i)
          this->cov_masked_(k,j) = this->cov_masked_(j,k);
          this->inv_cov_masked_(k,j) = this->inv_cov_masked_(j,k);
        }
      };
      break;
    }
    case 10:
    {
      for (int i=0; i<static_cast<int>(table.n_rows); i++)
      {
        const int j = static_cast<int>(table(i,0));
        const int k = static_cast<int>(table(i,1));

        this->cov_masked_(j,k) = table(i,8) + table(i,9);
        this->inv_cov_masked_(j,k) = table(i,8) + table(i,9);

        if (j!=k)
        {
          // apply mask to off-diagonal covariance elements
          this->cov_masked_(j,k) *= this->get_mask(j);
          this->cov_masked_(j,k) *= this->get_mask(k);

          this->inv_cov_masked_(j,k) *= this->get_mask(j);
          this->inv_cov_masked_(j,k) *= this->get_mask(k);

          // m(i,j) = m(j,i)
          this->cov_masked_(k,j) = this->cov_masked_(j,k);
          this->inv_cov_masked_(k,j) = this->inv_cov_masked_(j,k);
        }
      }
      break;
    }
    default:
    {
      spdlog::critical(
        "{}: data format for covariance file = {} is invalid",
        "set_inv_cov", 
        COV
      );
      exit(1);
    }
  }

  // Apply the Hartlap factor before inversion is equivalent to apply it 
  // after inversion and use block-wise inversion to the whole matrix
  if ((like.is_cmb_bandpower == 1) && (like.is_cmb_kkkk_cov_from_sim == 1))
  {
    int N5x2pt = this->ndata_ - like.Nbp;
    for (int i=N5x2pt; i<this->ndata_; i++)
    {
      for (int j=N5x2pt; j<this->ndata_; j++)
      {
        this->inv_cov_masked_(i,j) /= like.alpha_Hartlap_kkkk;
      }
    }
  }

  Vector eigvals = arma::eig_sym(this->inv_cov_masked_);
  for(int i=0; i<this->ndata_; i++)
  {
    if(eigvals(i) < 0)
    {
      spdlog::critical("{}: covmat not positive definite!", "set_inv_cov");
      exit(-1);
    }
  }

  this->inv_cov_masked_ = arma::inv(this->inv_cov_masked_);

  // apply mask again to make sure numerical errors in matrix inversion don't 
  // cause problems. Also, set diagonal elements corresponding to datavector
  // elements outside mask to 0, so that they don't contribute to chi2
  for (int i=0; i<this->ndata_; i++)
  {
    this->inv_cov_masked_(i,i) *= this->get_mask(i)*this->get_mask(i);
    for (int j=0; j<i; j++)
    {
      this->inv_cov_masked_(i,j) *= this->get_mask(i)*this->get_mask(j);
      this->inv_cov_masked_(j,i) = this->inv_cov_masked_(i,j);
    }
  };
  
  this->cov_filename_ = COV;
  
  this->is_inv_cov_set_ = true;

  this->cov_masked_sqzd_.set_size(this->ndata_masked_, this->ndata_masked_);

  this->inv_cov_masked_sqzd_.set_size(this->ndata_masked_, this->ndata_masked_);

  for(int i=0; i<this->ndata_; i++)
  {
    for(int j=0; j<this->ndata_; j++)
    {
      if((this->mask_(i)>0.99) && (this->mask_(j)>0.99))
      {
        if(this->get_index_sqzd(i) < 0)
        {
          spdlog::critical(
            "\x1b[90m{}\x1b[0m: logical error, internal"
            " inconsistent mask operation", 
            "set_inv_cov"
          );
          exit(1);
        }
        if(this->get_index_sqzd(j) < 0)
        {
          spdlog::critical(
            "\x1b[90m{}\x1b[0m: logical error, internal"
            " inconsistent mask operation", 
            "set_inv_cov"
          );
          exit(1);
        }

        this->cov_masked_sqzd_(this->get_index_sqzd(i),this->get_index_sqzd(j)) 
          = this->cov_masked_(i,j);

        this->inv_cov_masked_sqzd_(this->get_index_sqzd(i),this->get_index_sqzd(j)) 
          = this->inv_cov_masked_(i,j);
      }
    }
  }
}



// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int IP::get_mask(const int ci) const
{
  if (ci > like.Ndata || ci < 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: index i = {} is not valid (min = {}, max = {})",
      "get_mask", 
      ci, 
      0, 
      like.Ndata
    );
    exit(1);
  }

  return this->mask_(ci);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int IP::get_index_sqzd(const int ci) const
{
  if (ci > like.Ndata || ci < 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: index i = {} is not valid (min = {}, max = {})", 
      "get_index_sqzd", 
      ci, 
      0, 
      like.Ndata
    );
    exit(1);
  }

  return this->index_sqzd_(ci);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double IP::get_dv_masked(const int ci) const
{
  if (ci > like.Ndata || ci < 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: index i = {} is not valid (min = {}, max = {})",
      "get_dv_masked", ci, 0, like.Ndata);
    exit(1);
  }

  return this->data_masked_(ci);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double IP::get_dv_masked_sqzd(const int ci) const
{
  if (ci > like.Ndata || ci < 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: index i = {} is not valid (min = {}, max = {})",
      "get_dv_masked_sqzd", ci, 0, like.Ndata);
    exit(1);
  }

  return this->data_masked_sqzd_(ci);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double IP::get_inv_cov_masked(
    const int ci, 
    const int cj
  ) const
{
  if (ci > like.Ndata || ci < 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: index i = {} is not valid (min = {}, max = {})",
      "get_inv_cov_masked", 
      ci, 
      0, 
      like.Ndata
    );
    exit(1);
  }

  if (cj > like.Ndata || cj < 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: index j = {} is not valid (min = {}, max = {})",
      "get_inv_cov_masked", 
      cj, 
      0, 
      like.Ndata
    );
    exit(1);
  }

  return this->inv_cov_masked_(ci, cj);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double IP::get_inv_cov_masked_sqzd(
    const int ci, 
    const int cj
  ) const
{
  if (ci > like.Ndata || ci < 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: index i = {} is not valid (min = {}, max = {})",
      "get_inv_cov_masked_sqzd", 
      ci, 
      0, 
      like.Ndata
    );
    exit(1);
  }
  
  if (cj > like.Ndata || cj < 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: index j = {} is not valid (min = {}, max = {})",
      "get_inv_cov_masked_sqzd", 
      cj, 
      0, 
      like.Ndata
    );
    exit(1);
  }

  return this->inv_cov_masked_sqzd_(ci, cj);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int IP::get_ndata() const
{
  return this->ndata_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int IP::get_ndata_sqzd() const
{
  return this->ndata_masked_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

arma::Col<uword> IP::get_mask() const
{
  return this->mask_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

Vector IP::get_dv_masked() const
{
  return this->data_masked_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

Vector IP::get_dv_masked_sqzd() const
{
  return this->data_masked_sqzd_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

Matrix IP::get_inv_cov_masked() const
{
  return this->inv_cov_masked_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

Matrix IP::get_inv_cov_masked_sqzd() const
{
  return this->inv_cov_masked_sqzd_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

Matrix IP::get_cov_masked() const
{
  return this->cov_masked_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double IP::get_chi2(STLVector datavector) const
{
  if (!(this->is_data_set_))
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} not set prior to this function call",
      "get_chi2",
      "data_vector"
    );
    exit(1);
  }
  if (!(this->is_mask_set_))
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} not set prior to this function call",
      "get_chi2",
      "mask"
    );
    exit(1);
  }
  if (!(this->is_inv_cov_set_))
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} not set prior to this function call",
      "get_chi2",
      "inv_cov"
    );
    exit(1);
  }

  double chi2 = 0.0;
  
  for (int i=0; i<like.Ndata; i++)
  {
    if (this->get_mask(i))
    {
      const double x = datavector.at(i) - this->get_dv_masked(i);
      
      for (int j=0; j<like.Ndata; j++)
      {
        if (this->get_mask(j))
        {
          const double y = datavector.at(j) - this->get_dv_masked(j);
          chi2 += x*this->get_inv_cov_masked(i,j)*y;
        }
      }
    }
  }

  if (chi2 < 0.0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: chi2 = {} (invalid)", 
      "get_chi2", 
      chi2
    );
    exit(1);
  }

  return chi2;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

Matrix IP::get_cov_masked_sqzd() const
{
  return this->cov_masked_sqzd_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

Vector IP::expand_dim_from_masked_sqzd(Vector input) const
{
  if (this->ndata_masked_ != static_cast<int>(input.n_elem))
  {
    spdlog::critical("\x1b[90m{}\x1b[0m: {} invalid input vector",
      "expand_dim_from_masked_sqzd"
    );
    exit(1);
  }

  Vector result(this->ndata_, arma::fill::zeros);

  for(int i=0; i<this->ndata_; i++)
  {
    if(this->mask_(i) > 0.99)
    {
      if(this->get_index_sqzd(i) < 0)
      {
        spdlog::critical(
          "\x1b[90m{}\x1b[0m: logical error, internal "
          "inconsistent mask operation",
          "expand_dim_from_masked_sqzd"
        );
        exit(1);
      }

      result(i) = input(this->get_index_sqzd(i));
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
// Class IPCMB MEMBER FUNCTIONS
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

bool IPCMB::is_cmb_binmat_set() const
{
  return this->is_cmb_binmat_set_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

bool IPCMB::is_cmb_offset_set() const
{
  return this->is_cmb_offset_set_;
}

void IPCMB::set_cmb_binning_mat(std::string BINMAT)
{
  if(like.is_cmb_bandpower == 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} == 0, incompatible choice", 
      "set_cmb_binning_mat", 
      "like.is_cmb_bandpower"
    );
    exit(1);
  }

  if (!(like.Nbp>0))
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} not set prior to this call",
      "set_cmb_binning_mat", 
      "like.Nbp"
    );
    exit(1);
  } 

  if (!((like.lmin_bp>=0) && (like.lmax_bp>0)))
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} and {} not set prior to this call",
      "set_cmb_binning_mat", 
      "like.lmin_bp", 
      "like.lmax_bp"
    );
    exit(1);
  }

  this->nbp_ = like.Nbp;
  
  this->ncl_ = like.lmax_bp - like.lmin_bp + 1;
  
  this->cmb_binning_matrix_with_correction_.set_size(this->nbp_, this->ncl_);
  
  this->binmat_filename_ = BINMAT;

  Matrix table = read_table(BINMAT);
  for (int i=0; i<this->nbp_; i++)
  {
    for (int j=0; j<this->ncl_; j++)
    {
      this->cmb_binning_matrix_with_correction_(i,j) = table(i,j);
    }
  }

  spdlog::info(
    "\x1b[90m{}\x1b[0m: binning matrix file {} has {} x {} elements",
    "set_cmb_binning_mat", 
    BINMAT, 
    this->nbp_, 
    this->ncl_
  );

  this->is_cmb_binmat_set_ = true;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void IPCMB::set_cmb_theory_offset(std::string OFFSET)
{
  if (!(like.Nbp>0))
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: {} not set prior to this call", 
      "set_cmb_theory_offset", 
      "like.Nbp"
    );
    exit(1);
  }

  this->nbp_ = like.Nbp;
  
  this->cmb_theory_offset_.set_size(this->nbp_);
  
  this->offset_filename_ = OFFSET;

  Matrix table = read_table(OFFSET);
  for (int i=0; i<this->nbp_; i++)
  {
    this->cmb_theory_offset_(i) = static_cast<double>(table(i,0));
  }

  spdlog::info(
    "\x1b[90m{}\x1b[0m: {} has {} elements", 
    "set_cmb_theory_offset", 
    OFFSET, 
    this->nbp_
  );

  this->is_cmb_offset_set_ = true;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double IPCMB::get_binning_matrix_with_correction(
    const int ci, 
    const int cj
  ) const
{
  if (ci > this->nbp_ || ci < 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: index i = {} is not valid (min = {}, max = {})",
      "get_binning_matrix_with_correction", 
      ci, 
      0, 
      this->nbp_
    );
    exit(1);
  }
  
  if (cj > this->ncl_ || cj < 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: index j = {} is not valid (min = {}, max = {})",
      "get_binning_matrix_with_correction", 
      cj, 
      0, 
      this->ncl_
    );
    exit(1);
  }

  return this->cmb_binning_matrix_with_correction_(ci, cj);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double IPCMB::get_cmb_theory_offset(const int ci) const
{
  if (ci > this->nbp_ || ci < 0)
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: index i = {} is not valid (min = {}, max = {})",
      "get_cmb_theory_offset", ci, 0.0, this->nbp_);
    exit(1);
  }
  return this->cmb_theory_offset_(ci);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

Matrix IPCMB::get_binning_matrix_with_correction() const
{
  return this->cmb_binning_matrix_with_correction_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

Matrix IPCMB::get_cmb_theory_offset() const
{
  return this->cmb_theory_offset_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Class PointMass MEMBER FUNCTIONS
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void PointMass::set_pm_vector(STLVector pm)
{
  this->pm_ = pm;
  return;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

STLVector PointMass::get_pm_vector() const
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

  return 4*G_over_c2*this->pm_[zl]*1.e+13*g_tomo(a_lens, zs)/(theta*theta)/(chi_lens*a_lens);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// cia GLOBAL FUNCTIONS
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

Matrix read_table(const std::string file_name)
{
  std::ifstream input_file(file_name);

  if (!input_file.is_open())
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: file {} cannot be opened",
      "read_table",
      file_name
    );
    exit(1);
  }

  // --------------------------------------------------------
  // Read the entire file into memory
  // --------------------------------------------------------

  std::string tmp;
  
  input_file.seekg(0,std::ios::end);
  
  tmp.resize(static_cast<size_t>(input_file.tellg()));
  
  input_file.seekg(0,std::ios::beg);
  
  input_file.read(&tmp[0],tmp.size());
  
  input_file.close();
  
  if(tmp.empty())
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: file {} is empty",
      "read_table",
      file_name
    );
    exit(1);
  }
  
  // --------------------------------------------------------
  // Second: Split file into lines
  // --------------------------------------------------------
  
  std::vector<std::string> lines;
  lines.reserve(50000);

  boost::trim_if(tmp,boost::is_any_of("\t "));
  
  boost::trim_if(tmp,boost::is_any_of("\n"));
  
  boost::split(lines, tmp,boost::is_any_of("\n"), boost::token_compress_on);
  
  // Erase comment/blank lines
  auto check = [](std::string mystr) -> bool
  {
    return boost::starts_with(mystr, "#");
  };
  lines.erase(std::remove_if(lines.begin(), lines.end(), check), lines.end());
  
  // --------------------------------------------------------
  // Third: Split line into words
  // --------------------------------------------------------

  Matrix result;
  size_t ncols = 0;
  
  { // first line
    std::vector<std::string> words;
    words.reserve(100);
    
    boost::split(
      words,lines[0], 
      boost::is_any_of(" \t"),
      boost::token_compress_on
    );
    
    ncols = words.size();
    
    result.set_size(lines.size(), ncols);
    
    for (size_t j=0; j<ncols; j++)
    {
      result(0,j) = std::stod(words[j]);
    }
  }

  #pragma omp parallel for
  for (size_t i=1; i<lines.size(); i++)
  {
    std::vector<std::string> words;
    
    boost::split(
      words, 
      lines[i], 
      boost::is_any_of(" \t"),
      boost::token_compress_on
    );
    
    if (words.size() != ncols)
    {
      spdlog::critical(
        "\x1b[90m{}\x1b[0m: file {} is not well formatted"
        " (regular table required)", 
        "read_table", 
        file_name
      );
      exit(1);
    }
    
    for (size_t j=0; j<ncols; j++)
    {
      result(i,j) = std::stod(words[j]);
    }
  };
  
  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

STLVector convert_arma_col_to_stl_vector(Vector in)
{
  STLVector out(in.n_elem, 0.0);

  for(int i=0; i<static_cast<int>(in.n_elem); i++)
  {
    out[i] = in(i);
  }

  return out;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// BaryonScenario MEMBER FUNCTIONS
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int BaryonScenario::nscenarios() const
{
  return this->nscenarios_;
}

std::string BaryonScenario::get_scenario(const int i) const
{
  return this->scenarios_.at(i);
}

void BaryonScenario::set_scenarios(std::string scenarios)
{
  std::vector<std::string> lines;
  lines.reserve(50);

  // First: trim tab or newline
  boost::trim_if(scenarios, boost::is_any_of("\t "));
  
  boost::trim_if(scenarios, boost::is_any_of("\n"));

  if (scenarios.empty())
  {
    spdlog::critical(
      "\x1b[90m{}\x1b[0m: invalid string input (empty)",
      "init_baryon_pca_scenarios"
    );
    exit(1);
  }
  
  spdlog::info(
    "\x1b[90m{}\x1b[0m: Registering baryon scenarios for PCA", 
    "set_scenarios"
  );

  // --------------------------------------------------------
  // Second: Split scenarios into lines
  // --------------------------------------------------------

  boost::split(
    lines, 
    scenarios, 
    boost::is_any_of("/"), 
    boost::token_compress_on
  );
  
  // --------------------------------------------------------
  // Third: Expand abbreviated scenarios
  // --------------------------------------------------------

  int nscenarios = 0;
  
  for (auto it=lines.begin(); it != lines.end(); ++it)
  {
    // --------------------------------------------------------
    // Count occurrences of - (dashes)
    // --------------------------------------------------------
    size_t pos = 0;
    size_t count = 0; 
    while ((pos = (*it).find("-", pos)) != std::string::npos) 
    {
      ++count;
      ++pos;
    }

    if (count == 2)
    {
      // --------------------------------------------------------
      // Found an abbreviated scenario (e.g. antilles-3-12)
      // --------------------------------------------------------

      std::string sim_name = (*it).substr(0, (*it).find("-", 0));
      
      size_t pstart = (*it).find("-", 0) + 1; 
      
      size_t pend = (*it).find("-", pstart) + 1;
      
      int start = std::stoi((*it).substr(pstart, pend-pstart-1));
      
      int end = std::stoi((*it).substr(pend));
      
      for (int i=start; i<=end; i++)
      {
        this->scenarios_[nscenarios++] = sim_name + "-" + std::to_string(i);
        
        spdlog::info(
          "\x1b[90m{}\x1b[0m: Scenario {} = {}", 
          "set_scenarios", 
          nscenarios, sim_name + "-" + std::to_string(i)
        );
      }
    } 
    else if (count > 2)
    {
      // --------------------------------------------------------
      // Sanity check, there can't be more than two "-"
      // --------------------------------------------------------

      spdlog::critical(
        "\x1b[90m{}\x1b[0m: {} = {} probe not supported",
        "set_scenarios", 
        "scenario", 
        *it
      );
      exit(1);
    } 
    else
    {
      this->scenarios_[nscenarios++] = *it;
      
      spdlog::info(
        "\x1b[90m{}\x1b[0m: Scenario {} = {}", 
        "set_scenarios", 
        nscenarios, 
        *it
      );
    }
  }

  this->nscenarios_ = nscenarios;

  spdlog::info(
    "\x1b[90m{}\x1b[0m: {} scenarios are registered", 
    "set_scenarios", 
    this->nscenarios_
  );

  spdlog::info(
    "\x1b[90m{}\x1b[0m: Registering baryon scenarios for PCA done!", 
    "set_scenarios"
  );
  
  return;
}

} // end namespace cosmolike_interface

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------