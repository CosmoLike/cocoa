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
#include "cosmolike/interface_masked_data.hpp"
#include "cosmolike/structs.h"

namespace cosmolike_interface
{

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void init_probes(std::string possible_probes)
{
  boost::trim_if(possible_probes, boost::is_any_of("\t "));
  possible_probes = boost::algorithm::to_lower_copy(possible_probes);

  spdlog::debug("{}: Begins", "init_probes");

  if (possible_probes.compare("xi") == 0)
  {
    like.shear_shear = 1;

    spdlog::debug("{}: {} = {} selected", 
      "init_probes", "possible_probes", "xi");
  }
  else if (possible_probes.compare("wtheta") == 0)
  {
    like.pos_pos = 1;

    spdlog::debug("{}: {} = {} selected", 
      "init_probes", "possible_probes", "wtheta");
  }
  else if (possible_probes.compare("gammat") == 0)
  {
    like.shear_pos = 1;

    spdlog::debug("{}: {} = {} selected", 
      "init_probes", "possible_probes", "gammat");
  }
  else if (possible_probes.compare("2x2pt") == 0)
  {
    like.shear_pos = 1;
    like.pos_pos = 1;

    spdlog::debug("{}: {} = {} selected", 
      "init_probes", "possible_probes", "2x2pt");
  }
  else if (possible_probes.compare("3x2pt") == 0)
  {
    like.shear_shear = 1;
    like.shear_pos = 1;
    like.pos_pos = 1;

    spdlog::debug("{}: {} = {} selected", 
      "init_probes", "possible_probes", "3x2pt");
  }
  else if (possible_probes.compare("xi_ggl") == 0)
  {
    like.shear_shear = 1;
    like.shear_pos = 1;

    spdlog::debug("{}: {} = {} selected", 
      "init_probes", "possible_probes", "xi + ggl (2x2pt)");
  }
  else if (possible_probes.compare("xi_gg") == 0)
  {
    like.shear_shear = 1;
    like.pos_pos = 1;

    spdlog::debug("{}: {} = {} selected", 
      "init_probes", "possible_probes", "xi + ggl (2x2pt)");
  }
  else if (possible_probes.compare("5x2pt") == 0)
  {
    like.shear_shear = 1;
    like.shear_pos = 1;
    like.pos_pos = 1;
    like.gk = 1;
    like.ks = 1;
    
    spdlog::debug("{}: {} = {} selected", 
      "init_probes", "possible_probes", "5x2pt");
  }
  else if (possible_probes.compare("6x2pt") == 0)
  {
    like.shear_shear = 1;
    like.shear_pos = 1;
    like.pos_pos = 1;
    like.gk = 1;
    like.ks = 1;
    like.kk = 1;
    
    spdlog::debug("{}: {} = {} selected", 
      "init_probes", "possible_probes", "6x2pt");
  }
  else if (possible_probes.compare("c3x2pt") == 0)
  {
    like.gk = 1;
    like.ks = 1;
    like.kk = 1;
    
    spdlog::debug("{}: {} = {} selected", 
      "init_probes", "possible_probes", "c3x2pt (gk + sk + kk)");
  }
  else
  {
    spdlog::critical("{}: {} = {} probe not supported",
      "init_probes", "possible_probes", possible_probes);
    exit(1);
  }

  spdlog::debug("{}: Ends", "init_probes");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void init_binning_fourier(const int Nells, const int lmin, const int lmax)
{
  spdlog::debug("{}: Begins", "init_binning_fourier");

  if (!(Nells > 0))
  {
    spdlog::critical("{}: {} = {} not supported", "init_binning_fourier", 
      "Number of l modes (Nells)", Nells);
    exit(1);
  }

  spdlog::debug("{}: {} = {} selected.", "init_binning_fourier", "Nells", Nells);
  spdlog::debug("{}: {} = {} selected.", "init_binning_fourier", "l_min", lmin);
  spdlog::debug("{}: {} = {} selected.", "init_binning_fourier", "l_max", lmax);

  like.Ncl  = Nells;
  like.lmin = lmin;
  like.lmax = lmax;
  
  const double logdl = (std::log(like.lmax) - std::log(like.lmin))/like.Ncl;
  
  like.ell = (double*) malloc(sizeof(double)*like.Ncl);
  
  for (int i=0; i<like.Ncl; i++)
  {
    like.ell[i] = std::exp(std::log(like.lmin) + (i + 0.5)*logdl);
  
    spdlog::debug("{}: Bin {:d} - {} = {:.4e}, {} = {:.4e} and {} = {:.4e}",
      "init_binning_fourier", i, "lmin", lmin, "ell", like.ell[i], "lmax", lmax);
  }

  spdlog::debug("{}: Ends", "init_binning_fourier");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void init_binning_real_space(const int Ntheta, 
    const double theta_min_arcmin, const double theta_max_arcmin)
{
  spdlog::debug("{}: Begins", "init_binning_real_space");

  if (!(Ntheta > 0))
  {
    spdlog::critical("{}: {} = {} not supported", "init_binning_real_space",
      "NthetaNtheta", Ntheta);
    exit(1);
  }

  spdlog::debug("{}: {} = {} selected.", "init_binning_real_space", 
    "Ntheta", Ntheta);
  spdlog::debug("{}: {} = {} selected.", "init_binning_real_space", 
    "theta_min_arcmin", theta_min_arcmin);
  spdlog::debug("{}: {} = {} selected.", "init_binning_real_space", 
    "theta_max_arcmin", theta_max_arcmin);

  Ntable.Ntheta = Ntheta;
  Ntable.vtmin  = theta_min_arcmin * 2.90888208665721580e-4; // arcmin to rad conv
  Ntable.vtmax  = theta_max_arcmin * 2.90888208665721580e-4; // arcmin to rad conv
  
  const double logdt=(std::log(Ntable.vtmax)-std::log(Ntable.vtmin))/Ntable.Ntheta;
  
  for (int i=0; i<Ntable.Ntheta; i++)
  {
    const double thetamin = std::exp(log(Ntable.vtmin) + (i + 0.) * logdt);
    const double thetamax = std::exp(log(Ntable.vtmin) + (i + 1.) * logdt);
    
    const double theta = (2./3.) * (std::pow(thetamax,3) - std::pow(thetamin,3)) /
                                   (thetamax*thetamax    - thetamin*thetamin);

    spdlog::debug("{}: Bin {:d} - {} = {:.4e}, {} = {:.4e} and {} = {:.4e}",
      "init_binning_real_space", i, "theta_min [rad]", thetamin, "theta [rad]", 
      theta, "theta_max [rad]", thetamax);
  }

  spdlog::debug("{}: Ends", "init_binning_real_space");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void init_data_3x2pt_real_space(std::string cov, std::string mask, 
  std::string data, arma::Col<int>::fixed<3> order)
{
  spdlog::debug("{}: Begins", "init_data_3x2pt_real_space");

  init_data_vector_size_3x2pt_real_space();
  
  IP::get_instance().set_mask(mask, order);   // set_mask must be called first
  
  IP::get_instance().set_data(data);
  
  IP::get_instance().set_inv_cov(cov);

  spdlog::debug("{}: Ends", "init_data_3x2pt_real_space");

  return;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void init_data_6x2pt_real_space(std::string cov, std::string mask, 
  std::string data, arma::Col<int>::fixed<6> order)
{
  spdlog::debug("{}: Begins", "init_data_6x2pt_real_space");

  init_data_vector_size_6x2pt_fourier_space();
  
  IP::get_instance().set_mask(mask, order);   // set_mask must be called first
  
  IP::get_instance().set_data(data);
  
  IP::get_instance().set_inv_cov(cov);

  spdlog::debug("{}: Ends", "init_data_6x2pt_real_space");

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

void IP::set_data(std::string datavector_filename)
{
  if (!(this->is_mask_set_))
  {
    spdlog::critical(
      "{}: {} not set prior to this function call","set_data","mask");
    exit(1);
  }

  this->data_masked_.set_size(this->ndata_);
    
  this->data_masked_sqzd_.set_size(this->ndata_sqzd_);

  this->data_filename_ = datavector_filename;

  matrix table = read_table(datavector_filename);
  
  if (table.n_rows != this->ndata_)
  {
    spdlog::critical("{}: inconsistent data vector", "IP::set_data");
    exit(1);
  }

  for(int i=0; i<like.Ndata; i++)
  {
    this->data_masked_(i) = table(i,1);
    
    this->data_masked_(i) *= this->get_mask(i);

    if(this->get_mask(i) == 1)
    {
      if(this->get_index_sqzd(i) < 0)
      {
        spdlog::critical("{}: logical error, internal"
          " inconsistent mask operation", "IP::set_data");
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

void IP::set_mask(std::string mask_filename, arma::Col<int>::fixed<3> order)
{
  if (!(like.Ndata>0))
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "IP::set_mask", "like.Ndata");
    exit(1);
  }

  constexpr int sz = 3;

  this->ndata_ = like.Ndata;
  
  this->mask_.set_size(this->ndata_);

  this->mask_filename_ = mask_filename;

  matrix table = read_table(mask_filename);

  if (table.n_rows != this->ndata_)
  {
    spdlog::critical("{}: inconsistent mask", "IP::set_mask");
    exit(1);
  }
  
  for (int i=0; i<this->ndata_; i++)
  {
    this->mask_(i) = static_cast<int>(table(i,1) + 1e-13);
    
    if (!(0 == this->mask_(i) || 1 == this->mask_(i)))
    {
      spdlog::critical("{}: inconsistent mask", "IP::set_mask");
      exit(1);
    }
  }

  auto indices = arma::conv_to<arma::Col<int>>::from(
      arma::stable_sort_index(order, "ascend")
    );

  arma::Col<int>::fixed<sz> sizes =
    {
      2*Ntable.Ntheta*tomo.shear_Npowerspectra,
      Ntable.Ntheta*tomo.ggl_Npowerspectra,
      Ntable.Ntheta*tomo.clustering_Npowerspectra,
    };

  arma::Col<int>::fixed<sz> start = {0,0,0};

  for(int i=0; i<sz; i++)
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

  this->ndata_sqzd_ = arma::accu(this->mask_);
  
  if(!(this->ndata_sqzd_>0))
  {
    spdlog::critical("{}: mask file {} left no data points after masking",
      "IP::set_mask", mask_filename);
    exit(1);
  }
  
  spdlog::debug("{}: mask file {} left {} non-masked elements after masking",
    "IP::set_mask", mask_filename, this->ndata_sqzd_);

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

    if(j != this->ndata_sqzd_)
    {
      spdlog::critical("{}: logical error, internal inconsistent mask operation",
          "IP::set_mask");
      exit(1);
    }
  }

  this->is_mask_set_ = true;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void IP::set_mask(std::string mask_filename, arma::Col<int>::fixed<6> order)
{
  if (!(like.Ndata>0))
  {
    spdlog::critical(
        "{}: {} not set prior to this function call",
        "IP::set_mask", 
        "like.Ndata"
      );
    exit(1);
  }

  constexpr int sz = 6;

  this->ndata_ = like.Ndata;

  this->mask_.set_size(this->ndata_);

  this->mask_filename_ = mask_filename;

  matrix table = read_table(mask_filename);
  
  if (table.n_rows != this->ndata_)
  {
    spdlog::critical("{}: inconsistent mask", "IP::set_mask");
    exit(1);
  }

  for (int i=0; i<this->ndata_; i++)
  {
    this->mask_(i) = static_cast<int>(table(i,1)+1e-13);
    
    if (!(this->mask_(i) == 0 || this->mask_(i) == 1))
    {
      spdlog::critical("{}: inconsistent mask", "IP::set_mask");
      exit(1);
    }
  }

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
      IPCMB::get_instance().is_bandpower()  == 1 ? 
        IPCMB::get_instance().get_nbandpower() : 
        IPCMB::get_instance().get_ncl()
    };

  arma::Col<int>::fixed<sz> start = {0,0,0,0,0,0};

  for(int i=0; i<sz; i++)
    for(int j=0; j<indices(i); j++)
      start(i) += sizes(indices(j));

  if (like.shear_shear == 0)
  {
    const int N = start(0);
    const int M = N + sizes(0);
    for (int i=N; i<M; i++)
      this->mask_(i) = 0;
  }

  if (like.shear_pos == 0)
  {
    const int N = start(1);
    const int M = N + sizes(1);
    for (int i=N; i<M; i++)
      this->mask_(i) = 0;
  }

  if (like.pos_pos == 0)
  {
    const int N = start(2);
    const int M = N + sizes(2);
    for (int i=N; i<M; i++)
      this->mask_(i) = 0;
  }

  if (like.gk == 0) 
  {
    const int N = start(3);
    const int M = N + sizes(3);;
    for (int i=N; i<M; i++) 
      this->mask_(i) = 0.0;
  }

  if (like.ks == 0) 
  {
    const int N = start(4);
    const int M = N + sizes(4);
    for (int i=N; i<M; i++) 
      this->mask_(i) = 0.0;
  }

  if (like.kk == 0) 
  {
    const int N = start(5);
    const int M = N + sizes(5);
    for (int i=N; i<M; i++) 
      this->mask_(i) = 0.0;
  }
  
  this->ndata_sqzd_ = arma::accu(this->mask_);
  
  if(!(this->ndata_sqzd_>0))
  {
    spdlog::critical("{}: mask file {} left no data points after masking",
      "IP::set_mask", mask_filename);
    exit(1);
  }
  
  spdlog::debug("{}: mask file {} left {} non-masked elements "
    "after masking", "IP::set_mask",  mask_filename, this->ndata_sqzd_);

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
        this->index_sqzd_(i) = -1;
    }
    if(j != this->ndata_sqzd_)
    {
      spdlog::critical("{}: logical error, internal "
        "inconsistent mask operation", "IP::set_mask");
      exit(1);
    }
  }

  this->is_mask_set_ = true;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void IP::set_inv_cov(std::string covariance_filename)
{
  if (!(this->is_mask_set_))
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "IP::set_inv_cov", "mask");
    exit(1);
  }

  this->cov_filename_ = covariance_filename;

  matrix table = read_table(covariance_filename); 
  
  this->cov_masked_.set_size(this->ndata_, this->ndata_);
  this->cov_masked_.zeros();

  switch (table.n_cols)
  {
    case 3:
    {
      for (int i=0; i<static_cast<int>(table.n_rows); i++)
      {
        const int j = static_cast<int>(table(i,0));
        const int k = static_cast<int>(table(i,1));

        this->cov_masked_(j,k) = table(i,2);
        if (j!=k)
        {
          // apply mask to off-diagonal covariance elements
          this->cov_masked_(j,k) *= this->get_mask(j);
          this->cov_masked_(j,k) *= this->get_mask(k);
          // m(i,j) = m(j,i)
          this->cov_masked_(k,j) = this->cov_masked_(j,k);
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
        if (j!=k)
        {
          // apply mask to off-diagonal covariance elements
          this->cov_masked_(j,k) *= this->get_mask(j);
          this->cov_masked_(j,k) *= this->get_mask(k);
          // m(i,j) = m(j,i)
          this->cov_masked_(k,j) = this->cov_masked_(j,k);
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
        if (j!=k)
        {
          // apply mask to off-diagonal covariance elements
          this->cov_masked_(j,k) *= this->get_mask(j);
          this->cov_masked_(j,k) *= this->get_mask(k);
          // m(i,j) = m(j,i)
          this->cov_masked_(k,j) = this->cov_masked_(j,k);
        }
      }
      break;
    }
    default:
    {
      spdlog::critical("{}: data format for covariance file = {} is invalid",
        "IP::set_inv_cov", covariance_filename);
      exit(1);
    }
  }

  // Apply the Hartlap factor before inversion is equivalent to apply it 
  // after inversion and use block-wise inversion to the whole matrix
  if (IPCMB::get_instance().is_bandpower() && 
      IPCMB::get_instance().is_cov_kkkk_term_from_sim())
  {
    const int N5x2pt = this->ndata_ - IPCMB::get_instance().get_nbandpower();
    for (int i=N5x2pt; i<this->ndata_; i++)
      for (int j=N5x2pt; j<this->ndata_; j++)
        this->cov_masked_(i,j) /= IPCMB::get_instance().get_alpha_hartlap_cov_term_kkkk();
  }

  arma::Col<double> eigvals = arma::eig_sym(this->cov_masked_);
  for(int i=0; i<this->ndata_; i++)
  {
    if(eigvals(i) < 0)
    {
      spdlog::critical("{}: covmat (masked) not positive definite!","IP::set_inv_cov");
      exit(-1);
    }
  }

  this->inv_cov_masked_ = arma::inv(this->cov_masked_);

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
  
  this->cov_masked_sqzd_.set_size(this->ndata_sqzd_, this->ndata_sqzd_);

  this->inv_cov_masked_sqzd_.set_size(this->ndata_sqzd_, this->ndata_sqzd_);

  for(int i=0; i<this->ndata_; i++)
  {
    for(int j=0; j<this->ndata_; j++)
    {
      if((this->mask_(i)>0.99) && (this->mask_(j)>0.99))
      {
        if(this->get_index_sqzd(i) < 0)
        {
          spdlog::critical("{}: logical error, internal"
            " inconsistent mask operation", "IP::set_inv_cov");
          exit(1);
        }
        if(this->get_index_sqzd(j) < 0)
        {
          spdlog::critical("{}: logical error, internal"
            " inconsistent mask operation", "IP::set_inv_cov");
          exit(1);
        }

        this->cov_masked_sqzd_(this->get_index_sqzd(i),
            this->get_index_sqzd(j)) = this->cov_masked_(i,j);

        this->inv_cov_masked_sqzd_(this->get_index_sqzd(i),
            this->get_index_sqzd(j)) = this->inv_cov_masked_(i,j);
      }
    }
  }

  this->is_inv_cov_set_ = true;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int IP::get_mask(const int ci) const
{
  if (ci > like.Ndata || ci < 0)
  {
    spdlog::critical("{}: index i = {} is not valid (min = {}, max = {})",
      "IP::get_mask", ci, 0, like.Ndata);
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
    spdlog::critical("{}: index i = {} is not valid (min = {}, max = {})", 
      "IP::get_index_sqzd", ci, 0, like.Ndata);
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
    spdlog::critical("{}: index i = {} is not valid (min = {}, max = {})",
      "IP::get_dv_masked", ci, 0, like.Ndata);
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
    spdlog::critical("{}: index i = {} is not valid (min = {}, max = {})",
      "IP::get_dv_masked_sqzd", ci, 0, like.Ndata);
    exit(1);
  }
  return this->data_masked_sqzd_(ci);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double IP::get_inv_cov_masked(const int ci, const int cj) const
{
  if (ci > like.Ndata || ci < 0)
  {
    spdlog::critical("{}: index i = {} is not valid (min = {}, max = {})",
      "IP::get_inv_cov_masked", ci, 0, like.Ndata);
    exit(1);
  }
  if (cj > like.Ndata || cj < 0)
  {
    spdlog::critical("{}: index j = {} is not valid (min = {}, max = {})",
      "IP::get_inv_cov_masked", cj, 0, like.Ndata);
    exit(1);
  }
  return this->inv_cov_masked_(ci, cj);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double IP::get_inv_cov_masked_sqzd(const int ci, const int cj) const
{
  if (ci > like.Ndata || ci < 0)
  {
    spdlog::critical("{}: index i = {} is not valid (min = {}, max = {})",
      "IP::get_inv_cov_masked_sqzd", ci, 0, like.Ndata);
    exit(1);
  }
  if (cj > like.Ndata || cj < 0)
  {
    spdlog::critical("{}: index j = {} is not valid (min = {}, max = {})",
      "IP::get_inv_cov_masked_sqzd", cj, 0, like.Ndata);
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
  return this->ndata_sqzd_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

arma::Col<int> IP::get_mask() const
{
  return this->mask_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

arma::Col<double> IP::get_dv_masked() const
{
  return this->data_masked_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

arma::Col<double> IP::get_dv_masked_sqzd() const
{
  return this->data_masked_sqzd_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

matrix IP::get_inv_cov_masked() const
{
  return this->inv_cov_masked_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

matrix IP::get_inv_cov_masked_sqzd() const
{
  return this->inv_cov_masked_sqzd_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

matrix IP::get_cov_masked() const
{
  return this->cov_masked_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double IP::get_chi2(arma::Col<double> datavector) const
{
  if (!(this->is_data_set_))
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "IP::get_chi2", "data_vector");
    exit(1);
  }
  if (!(this->is_mask_set_))
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "IP::get_chi2", "mask");
    exit(1);
  }
  if (!(this->is_inv_cov_set_))
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "IP::get_chi2", "inv_cov");
    exit(1);
  }
  if (datavector.n_elem != like.Ndata)
  {
    spdlog::critical("{}: "
      "incompatible data vector (theory size = {}, data size = {})",
      "IP::get_chi2", datavector.n_elem, like.Ndata);
    exit(1);
  }

  double chi2 = 0.0;
  
  for (int i=0; i<like.Ndata; i++)
  {
    if (this->get_mask(i))
    {
      const double x = datavector(i) - this->get_dv_masked(i);
      
      for (int j=0; j<like.Ndata; j++)
      {
        if (this->get_mask(j))
        {
          const double y = datavector(j) - this->get_dv_masked(j);
          chi2 += x*this->get_inv_cov_masked(i,j)*y;
        }
      }
    }
  }

  if (chi2 < 0.0)
  {
    spdlog::critical("{}: chi2 = {} (invalid)", "IP::get_chi2", chi2);
    exit(1);
  }
  return chi2;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

matrix IP::get_cov_masked_sqzd() const
{
  return this->cov_masked_sqzd_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

arma::Col<double> IP::expand_theory_data_vector_from_sqzd(arma::Col<double> input) const
{
  if (this->ndata_sqzd_ != static_cast<int>(input.n_elem))
  {
    spdlog::critical("{}: invalid input data vector",
      "IP::expand_theory_data_vector_from_sqzd");
    exit(1);
  }

  arma::Col<double> result(this->ndata_,arma::fill::zeros);
  for(int i=0; i<this->ndata_; i++)
  {
    if(this->mask_(i) > 0.99)
    {
      if(this->get_index_sqzd(i) < 0)
      {
        spdlog::critical("{}: logical error, inconsistent mask operation",
          "IP::expand_theory_data_vector_from_sqzd");
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

arma::Col<double> IP::sqzd_theory_data_vector(arma::Col<double> input) const
{
  if (this->ndata_ != static_cast<int>(input.n_elem))
  {
    spdlog::critical("{}: invalid input data vector",
      "IP::sqzd_theory_data_vector");
    exit(1);
  }

  arma::Col<double> result(this->ndata_sqzd_, arma::fill::zeros);
  for (int i=0; i<this->ndata_; i++)
    if (this->get_mask(i) > 0.99)
      result(this->get_index_sqzd(i)) = input(i);
  return result;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
} // namespace cosmolike_interface