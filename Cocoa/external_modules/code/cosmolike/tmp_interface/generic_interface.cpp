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

// boost library
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/replace.hpp>

// cosmolike
#include "cosmolike/basics.h"
#include "cosmolike/bias.h"
#include "cosmolike/baryons.h"
#include "cosmolike/cosmo2D.h"
#include "cosmolike/cosmo3D.h"
#include "cosmolike/IA.h"
#include "cosmolike/halo.h"
#include "cosmolike/radial_weights.h"
#include "cosmolike/pt_cfastpt.h"
#include "cosmolike/redshift_spline.h"
#include "cosmolike/structs.h"



#include "cosmolike/generic_interface.hpp"

using vector = arma::Col<double>;
using matrix = arma::Mat<double>;
using cube = arma::Cube<double>;

// Why the cpp functions accept and return STL vectors (instead of arma:Col)?
// Answer: the conversion between STL vector and python np array is cleaner
// Answer: arma:Col is cast to 2D np array with 1 column (not as nice!)

namespace cosmolike_interface
{

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

void initial_setup()
{
  spdlog::cfg::load_env_levels();
  
  spdlog::debug("{}: Begins", "initial_setup");

  like.shear_shear = 0;
  like.shear_pos = 0;
  like.pos_pos = 0;

  like.lmin_kappacmb = 0;
  like.lmax_kappacmb = 0;

  like.is_bandpower = 0;
  like.lmin_bp = 0; 
  like.lmax_bp = 0;
  like.Nbp = 0;

  like.Ncl = 0;
  like.lmin = 0;
  like.lmax = 0;

  like.gk = 0;
  like.kk = 0;
  like.ks = 0;
  
    // no priors
  like.clusterN = 0;
  like.clusterWL = 0;
  like.clusterCG = 0;
  like.clusterCC = 0;

  // reset bias
  reset_redshift_struct();
  reset_nuisance_struct();
  reset_cosmology_struct();
  reset_tomo_struct();
  reset_Ntable_struct();
  reset_like_struct();

  like.adopt_limber_gg = 0;

  std::string mode = "Halofit";
  memcpy(pdeltaparams.runmode, mode.c_str(), mode.size() + 1);

  spdlog::debug("{}: Ends", "initial_setup");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void init_accuracy_boost(
    const double accuracy_boost, 
    const double sampling_boost,
    const int integration_accuracy
  )
{
   // (imp on notebooks where users can change sampling_boost interactively
  static int N_a = 0;
  static int N_ell = 0;

  if (N_a == 0) N_a = Ntable.N_a;
  Ntable.N_a = static_cast<int>(ceil(N_a*sampling_boost));
  
  if (N_ell == 0) N_ell = Ntable.N_ell;
  Ntable.N_ell = static_cast<int>(ceil(N_ell*sampling_boost));

  /*  
  Ntable.N_k_lin = 
    static_cast<int>(ceil(Ntable.N_k_lin*sampling_boost));
  
  Ntable.N_k_nlin = 
    static_cast<int>(ceil(Ntable.N_k_nlin*sampling_boost));

  Ntable.N_M = 
    static_cast<int>(ceil(Ntable.N_M*sampling_boost));
  */

  Ntable.high_def_integration = integration_accuracy;

  // update cache
  Ntable.random = RandomNumber::get_instance().get();
}





// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void init_cosmo_runmode(const bool is_linear)
{
  spdlog::debug("{}: Begins", "init_cosmo_runmode");

  std::string mode = is_linear ? "linear" : "Halofit";
  
  const size_t size = mode.size();
  
  memcpy(pdeltaparams.runmode, mode.c_str(), size + 1);

  spdlog::debug("{}: {} = {} selected", 
    "init_cosmo_runmode", "runmode", mode);
  
  spdlog::debug("{}: Ends", "init_cosmo_runmode");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void init_data_vector_size_real_space(arma::Col<int>::fixed<6> exclude)
{ // data vector size on Cosmolike: like.Ndata
  arma::Col<int>::fixed<6> ndv = 
    {
      Ntable.Ntheta*2*tomo.shear_Npowerspectra,
      Ntable.Ntheta*tomo.ggl_Npowerspectra,
      Ntable.Ntheta*tomo.clustering_Npowerspectra,
      Ntable.Ntheta*redshift.shear_nbin,
      Ntable.Ntheta*redshift.clustering_nbin,
      0
    };

  if (IPCMB::get_instance().is_bandpower() == 0)
    ndv(5) = like.Ncl;
  else
    ndv(5) = IPCMB::get_instance().get_nbandpower();

  like.Ndata = 0.0;

  for(int i=0; i<exclude.n_elem; i++)
    if (exclude(i) > 0)
      like.Ndata += ndv(i);
}



// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void init_data_vector_size_6x2pt_real_space()
{
  spdlog::debug("{}: Begins", "init_data_vector_size_6x2pt_real_space");

  if (tomo.shear_Npowerspectra == 0)
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "init_data_vector_size_6x2pt_real_space", "tomo.shear_Npowerspectra");
    exit(1);
  }
  if (tomo.ggl_Npowerspectra == 0)
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "init_data_vector_size_6x2pt_real_space", "tomo.ggl_Npowerspectra");
    exit(1);
  }
  if (tomo.clustering_Npowerspectra == 0)
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "init_data_vector_size_6x2pt_real_space", "tomo.clustering_Npowerspectra");
    exit(1);
  }
  if (redshift.shear_nbin == 0)
  {
    spdlog::critical( "{}: {} not set prior to this function call",
      "init_data_vector_size_6x2pt_real_space", "redshift.shear_nbin");
    exit(1);
  }
  if (redshift.clustering_nbin == 0)
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "init_data_vector_size_6x2pt_real_space", "redshift.clustering_nbin");
    exit(1);
  }
  if (Ntable.Ntheta == 0) 
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "init_data_vector_size_6x2pt_real_space", "Ntable.Ntheta");
    exit(1);
  }
  if (IPCMB::get_instance().is_bandpower() == 0)
  {
    if (like.Ncl == 0)
    {
      spdlog::critical("{}: {} not set prior to this function call", 
        "init_data_vector_size_6x2pt_real_space", "like.Ncl");
      exit(1);
    }
  }
  else if (IPCMB::get_instance().is_bandpower() == 1)
  {
    if (IPCMB::get_instance().get_nbandpower() == 0)
    {
      spdlog::critical("{}: {} not set prior to this function call", 
        "init_data_vector_size_6x2pt_real_space", 
        "IPCMB::get_instance().get_nbandpower()");
      exit(1);
    }    
  }
  else
  {
    spdlog::critical("{}: {} not set prior to this function call", 
      "init_data_vector_size_6x2pt_real_space", "is_bandpower");
  }

  arma::Col<int>::fixed<6> exclude = {1, 1, 1, 1, 1, 1};

  init_data_vector_size_real_space(exclude);

  spdlog::debug("{}: {} = {} selected.", 
    "init_data_vector_size_6x2pt_real_space", "Ndata", like.Ndata);
  
  spdlog::debug("{}: Ends", "init_data_vector_size_6x2pt_real_space");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void init_data_vector_size_fourier_space(
    arma::Col<int>::fixed<6> exclude
  )
{ // data vector size on Cosmolike: like.Ndata
  arma::Col<int>::fixed<6> ndv = 
    {
      like.Ncl*tomo.shear_Npowerspectra,
      like.Ncl*tomo.ggl_Npowerspectra,
      like.Ncl*tomo.clustering_Npowerspectra,
      like.Ncl*redshift.shear_nbin,
      like.Ncl*redshift.clustering_nbin,
      0
    };

  if (IPCMB::get_instance().is_bandpower() == 0)
    ndv(5) = like.Ncl;
  else
    ndv(5) = IPCMB::get_instance().get_nbandpower();

  like.Ndata = 0.0;

  for(int i=0; i<exclude.n_elem; i++)
    if (exclude(i) > 0)
      like.Ndata += ndv(i);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void init_data_vector_size_6x2pt_fourier_space()
{
  spdlog::debug("{}: Begins", 
      "init_data_vector_size_6x2pt_fourier_space");

  if (tomo.shear_Npowerspectra == 0)
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "init_data_vector_size_6x2pt_fourier_space", "tomo.shear_Npowerspectra");
    exit(1);
  }
  if (tomo.ggl_Npowerspectra == 0)
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "init_data_vector_size_6x2pt_fourier_space", "tomo.ggl_Npowerspectra");
    exit(1);
  }
  if (tomo.clustering_Npowerspectra == 0)
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "init_data_vector_size_6x2pt_fourier_space", "tomo.clustering_Npowerspectra");
    exit(1);
  }
  if (redshift.shear_nbin == 0)
  {
    spdlog::critical( "{}: {} not set prior to this function call",
      "init_data_vector_size_6x2pt_fourier_space", "redshift.shear_nbin");
    exit(1);
  }
  if (redshift.clustering_nbin == 0)
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "init_data_vector_size_6x2pt_fourier_space", "redshift.clustering_nbin");
    exit(1);
  }
  if (like.Ncl == 0) 
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "init_data_vector_size_6x2pt_real_space", "like.Ncl");
    exit(1);
  }
  if (IPCMB::get_instance().is_bandpower() == 0)
  {
    if (like.Ncl == 0)
    {
      spdlog::critical("{}: {} not set prior to this function call", 
        "init_data_vector_size_6x2pt_real_space", "like.Ncl");
      exit(1);
    }
  }
  else if (IPCMB::get_instance().is_bandpower() == 1)
  {
    if (IPCMB::get_instance().get_nbandpower() == 0)
    {
      spdlog::critical("{}: {} not set prior to this function call", 
        "init_data_vector_size_6x2pt_real_space", 
        "IPCMB::get_instance().get_nbandpower()");
      exit(1);
    }    
  }
  else
  {
    spdlog::critical("{}: {} not set prior to this function call", 
      "init_data_vector_size_6x2pt_real_space", "is_bandpower");
  }

  arma::Col<int>::fixed<6> exclude = {1, 1, 1, 1, 1, 1};

  init_data_vector_size_fourier_space(exclude);
  
  spdlog::debug("{}: {} = {} selected.", 
    "init_data_vector_size_6x2pt_fourier_space", "Ndata", like.Ndata);
  
  spdlog::debug("{}: Ends", 
    "init_data_vector_size_6x2pt_fourier_space");
}





// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void init_survey(
    std::string surveyname, 
    double area, 
    double sigma_e
  )
{
  spdlog::debug("{}: Begins", "init_survey");

  boost::trim_if(surveyname, boost::is_any_of("\t "));
  
  surveyname = boost::algorithm::to_lower_copy(surveyname);
    
  if (surveyname.size() > CHAR_MAX_SIZE - 1)
  {
    spdlog::critical(
      "{}: survey name too large for Cosmolike "
      "(C char memory overflow)", "init_survey");
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

  spdlog::debug("{}: Ends", "init_survey");
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

double get_baryon_power_spectrum_ratio(const double log10k, const double a)
{
  const double KNL = pow(10.0, log10k)*cosmology.coverH0;
  return PkRatio_baryons(KNL, a);
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


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void compute_gk_real_masked(vector& data_vector, const int start)
{
  if (like.gk == 1)
  {
    for (int nz=0; nz<redshift.clustering_nbin; nz++)
    {
      for (int i=0; i<Ntable.Ntheta; i++) 
      {
        const int index = start + Ntable.Ntheta*nz + i;
        
        if (IP::get_instance().get_mask(index))
        {
          data_vector(index) = w_gk_tomo(i, nz, 1);
        }
      }
    }
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void compute_ks_real_masked(vector& data_vector, const int start)
{
  if (like.ks == 1) 
  {
    for (int nz=0; nz<redshift.shear_nbin; nz++)
    {
      for (int i=0; i<Ntable.Ntheta; i++)
      {
        const int index = start + Ntable.Ntheta*nz + i; 
        
        if (IP::get_instance().get_mask(index))
        {
          data_vector(index) = 
            w_ks_tomo(i, nz, 1)*(1.0 + nuisance.shear_calibration_m[nz]);
        }
      }
    }
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void compute_kk_fourier_masked(vector& data_vector, const int start)
{
  if (like.kk == 1)
  {
    if (!IPCMB::get_instance().is_bandpower())
    {
      for (int i=0; i<IPCMB::get_instance().get_ncl(); i++)  
      {        
        if (IP::get_instance().get_mask(start+i))
        {
          const double l = like.ell[i];
          data_vector(start+i) = 
            (l<=limits.LMIN_tab) ? C_kk_limber_nointerp(l,0) : C_kk_limber(l);
        }
      }
    }
    else
    {
      const int Nbp = IPCMB::get_instance().get_nbandpower();
      std::array<double,2> lrange = IPCMB::get_instance().lrange_bandpower();
      
      for (int L=lrange[0]; L<lrange[1]+1; L++)
      {
        const double Ckk = (L <= limits.LMIN_tab) ? 
          C_kk_limber_nointerp((double) L, 0) : C_kk_limber((double) L);

        const int i = L-lrange[0];

        for (int j=0; j<Nbp; j++) // Loop through bandpower bins
          if (IP::get_instance().get_mask(start+j))
            data_vector(start+j) += (
              Ckk*IPCMB::get_instance().get_binned_matrix_with_correction(j,i));
      }

      for (int j=0; j<Nbp; j++) // offset due to marginalizing over primary CMB
        if (IP::get_instance().get_mask(start+j))
          data_vector(start+j) -= IPCMB::get_instance().get_theory_offset(j);
    }
  }
}   


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

vector compute_data_vector_6x2pt_real_masked_any_order(
    arma::Col<int>::fixed<6> order
  )
{ // order = (1,2,3,4,5,6) => Cosmic Shear, ggl, gg, gk, sk, kk
  spdlog::debug("{}: Begins", "compute_data_vector_6x2pt_real_masked_any_order");

  if (!(0 > redshift.shear_nbin))
  {
    spdlog::critical("{}: {} = 0 is invalid",
        "compute_data_vector_6x2pt_real_masked_any_order","redshift.shear_Nbin");
    exit(1);
  }
  if (!(0 > Ntable.Ntheta))
  {
    spdlog::critical("{}: {} = 0 is invalid",
      "compute_data_vector_6x2pt_real_masked_any_order","Ntable.Ntheta");
    exit(1);
  }
  if (IPCMB::get_instance().is_bandpower())
  {
    if (!(IPCMB::get_instance().get_nbandpower() > 0))
    {
      spdlog::critical("{}: {} not set prior to this function call",
        "compute_data_vector_6x2pt_real_masked_any_order", "IPCMB::nbandpower");
      exit(1);
    } 
    if (!(0 > IPCMB::get_instance().get_lrange_bandpower()[0]) or 
        !(0 > IPCMB::get_instance().get_lrange_bandpower()[1]))
    {
      spdlog::critical("{}: {} is invalid",
          "compute_data_vector_6x2pt_real_masked_any_order", 
          "IPCMB::lrange_bandpower");
      exit(1);
    }
    if (!IPCMB::get_instance().is_binmat_set())
    { // check binning matrix and CMB lensing band power offset
      spdlog::critical("{}: {} not set prior to this function call",
        "compute_data_vector_6x2pt_real_masked_any_order", 
        "IPCMB::cmb_binning_matrix_with_correction");
      exit(1);
    }
    if (!IPCMB::get_instance().is_offset_set())
    {
      spdlog::critical("{}: {} not set prior to this function call",
        "compute_data_vector_6x2pt_real_masked_any_order", "cmb_theory_offset");
      exit(1);
    }
  }
  else
  {
    if (!(IPCMB::get_instance().get_ncl() > 0))
    {
      spdlog::critical("{}: {} not set prior to this function call",
        "compute_data_vector_6x2pt_real_masked_any_order", "IPCMB::ncl");
      exit(1);
    }
    if (!(0 > IPCMB::get_instance().get_lrange_kappacmb()[0]) or 
        !(0 > IPCMB::get_instance().get_lrange_kappacmb()[1]))
    {
      spdlog::critical("{}: {} is invalid",
        "compute_data_vector_6x2pt_real_masked_any_order","IPCMB::lrange_kappacmb");
      exit(1);
    }
  }
  if (!IP::get_instance().is_mask_set())
  {
    spdlog::critical(
        "{}: {} not set prior to this function call",
        "compute_data_vector_6x2pt_real_masked_any_order", 
        "mask"
      );
    exit(1);
  }
  if (!(IPCMB::get_instance()get_lrange_kappacmb()[0]>0) || 
      !(IPCMB::get_instance()get_lrange_kappacmb()[1]>0))
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "compute_data_vector_6x2pt_real_masked_any_order", "IPCMB::lrange_kappacmb");
    exit(1);
  }
  if (cmb.fwhm <= 0)
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "compute_data_vector_6x2pt_real_masked_any_order", "cmb.fwhm");
    exit(1);
  }

  constexpr int sz = 6;

  auto indices = arma::conv_to<arma::Col<int>>::from(
      arma::stable_sort_index(order, "ascend")
    );

  arma::Col<int>::fixed<sz> sizes =
    {
      2*Ntable.Ntheta*tomo.shear_Npowerspectra,
      Ntable.Ntheta*tomo.ggl_Npowerspectra,
      Ntable.Ntheta*tomo.clustering_Npowerspectra,
      Ntable.Ntheta*redshift.clustering_nbin,
      Ntable.Ntheta*redshift.shear_nbin,
        IPCMB::get_instance().is_bandpower() == 1 ? 
        IPCMB::get_instance().get_nbandpower() : 
        IPCMB::get_instance().get_ncl();
    };

  arma::Col<int>::fixed<sz> start = {0,0,0,0,0,0};

  for(int i=0; i<sz; i++)
    for(int j=0; j<indices(i); j++)
      start(i) += sizes(indices(j));
  
  vector data_vector(like.Ndata, arma::fill::zeros);
  
  compute_ss_real_masked(data_vector, start(0));

  compute_gs_real_masked(data_vector, start(1));

  compute_gg_real_masked(data_vector, start(2));

  compute_gk_real_masked(data_vector, start(3));

  compute_ks_real_masked(data_vector, start(4));
  
  compute_kk_fourier_masked(data_vector, start(5));
  
  spdlog::debug("{}: Ends", "compute_data_vector_6x2pt_real_masked_any_order");

  return data_vector;
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Class IPCMB MEMBER FUNCTIONS
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

bool IPCMB::is_binmat_set() const
{
  return this->is_binmat_set_;
}

bool IPCMB::is_offset_set() const
{
  return this->is_offset_set_;
}

bool IPCMB::is_bandpower() const
{
  return this->is_bandpower_;
}

bool IPCMB::is_cov_kkkk_term_from_sim() const
{
  return this->is_cov_kkkk_term_from_sim_;
}

void IPCMB::set_binning_mat(std::string binned_matrix_filename)
{
  if(this->is_bandpower_ == 0)
  {
    spdlog::critical("{}: {} == 0, incompatible choice",
      "IPCMB::set_binning_mat", "IPCMB::is_bandpower");
    exit(1);
  }
  if (!(this->nbandpower_>0))
  {
    spdlog::critical("{}: {} not set prior to this call",
      "IPCMB::set_binning_mat", "IPCMB::nbandpower"); 
    exit(1);
  } 
  if (!(0 > this->lrange_bandpower_[0]) or !(0 > this->lrange_bandpower_[1]))
  {
    spdlog::critical("{}: {} not set prior to this call",
        "IPCMB::set_binning_mat", "IPCMB::lrange_bandpower");
    exit(1);
  }
  
  this->ncl_ = this->lrange_bandpower_[1] - this->lrange_bandpower_[0] + 1;
  
  this->binmat_with_correction_.set_size(this->nbandpower_,this->ncl_);
  
  this->binmat_filename_ = binned_matrix_filename;

  matrix table = read_table(binned_matrix_filename);
  #pragma omp parallel for
  for (int i=0; i<this->nbandpower_; i++)
    for (int j=0; j<this->ncl_; j++)
      this->binmat_with_correction_(i,j) = table(i,j);

  spdlog::debug("{}: CMB binning matrix from file {} has {} x {} elements",
    "IPCMB::set_binning_mat",binned_matrix_filename,this->nbandpower_,this->ncl_);

  this->is_binmat_set_ = true;
}

void IPCMB::set_theory_offset(std::string theory_offset_filename)
{
  if (!(like.Nbp>0))
  {
    spdlog::critical("{}: {} not set prior to this call", 
      "IPCMB::set_theory_offset", "like.Nbp");
    exit(1);
  }

  this->nbandpower_ = like.Nbp;
  
  this->theory_offset_.set_size(this->nbandpower_);
  
  this->offset_filename_ = theory_offset_filename;

  matrix table = read_table(theory_offset_filename);
  #pragma omp parallel for
  for (int i=0; i<this->nbandpower_; i++)
    this->theory_offset_(i) = static_cast<double>(table(i,0));

  spdlog::debug("{}: CMB theory offset from file {} has {} elements", 
    "IPCMB::set_theory_offset", theory_offset_filename, this->nbandpower_);

  this->is_offset_set_ = true;
}

double IPCMB::get_binned_matrix_with_correction(const int ci, const int cj) const
{
  if (ci > this->nbandpower_ || ci < 0)
  {
    spdlog::critical("{}: index i = {} is not valid (min = {}, max = {})",
      "IPCMB::get_binned_matrix_with_correction", ci, 0, this->nbandpower_);
    exit(1);
  }
  if (cj > this->ncl_ || cj < 0)
  {
    spdlog::critical("{}: index j = {} is not valid (min = {}, max = {})",
      "IPCMB::get_binned_matrix_with_correction", cj, 0, this->ncl_);
    exit(1);
  }
  return this->binmat_with_correction_(ci, cj);
}

double IPCMB::get_theory_offset(const int ci) const
{
  if (ci > this->nbandpower_ || ci < 0)
  {
    spdlog::critical("{}: index i = {} is not valid (min = {}, max = {})",
      "IPCMB::get_theory_offset", ci, 0.0, this->nbandpower_);
    exit(1);
  }
  return this->theory_offset_(ci);
}

matrix IPCMB::get_binned_matrix_with_correction() const
{
  return this->binmat_with_correction_;
}

matrix IPCMB::get_theory_offset() const
{
  return this->theory_offset_;
}

int IPCMB::get_nbandpower() const
{
  return this->nbandpower_;
}

double IPCMB::get_alpha_hartlap_cov_term_kkkk() const
{
  return this->alpha_hartlap_cov_kkkk_;
}

std::array<double,2> IPCMB::get_lrange_kappa() const
{
  return this->lrange_kappacmb_;
}

std::array<double,2> IPCMB::get_lrange_bandpower() const
{
  return this->lrange_bandpower_;
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

} // end namespace cosmolike_interface

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------