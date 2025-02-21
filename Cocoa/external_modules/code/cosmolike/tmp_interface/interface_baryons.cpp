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
#include "cosmolike/structs.h"
#include "cosmolike/interface_aux.hpp"
#include "cosmolike/interface_baryons.hpp"
#include "cosmolike/interface_datavector.hpp"

using vector = arma::Col<double>;
using matrix = arma::Mat<double>;
using cube   = arma::Cube<double>;

namespace cosmolike_interface
{

std::tuple<std::string,int> get_baryon_sim_name_and_tag(std::string sim)
{
  // Desired Convention:
  // (1) Python input: not be case sensitive
  // (2) simulation names only have "_" as deliminator, e.g., owls_AGN.
  // (3) simulation IDs are indicated by "-", e.g., antilles-1.
 
  boost::trim_if(sim, boost::is_any_of("\t "));
  sim = boost::algorithm::to_lower_copy(sim);
  
  { // Count occurrences of - (dashes)
    size_t pos = 0; 
    size_t count = 0; 
    std::string tmp = sim;

    while ((pos = tmp.rfind("-")) != std::string::npos) 
    {
      tmp = tmp.substr(pos+1);
      ++count;
    }

    if (count > 1)
    {
      const std::string etmp = "get_baryon_sim_name_and_tag";
      const std::string emsg = "{}: Scenario {} not supported (too many dashes)";
      spdlog::critical(emsg, etmp, sim);
      exit(1);
    }
  }

  if (sim.rfind("owls_agn") != std::string::npos)
  {
    boost::replace_all(sim, "owls_agn", "owls_AGN");
    boost::replace_all(sim, "_t80", "-1");
    boost::replace_all(sim, "_t85", "-2");
    boost::replace_all(sim, "_t87", "-3");
  } 
  else if (sim.rfind("bahamas") != std::string::npos)
  {
    boost::replace_all(sim, "bahamas", "BAHAMAS");
    boost::replace_all(sim, "_t78", "-1");
    boost::replace_all(sim, "_t76", "-2");
    boost::replace_all(sim, "_t80", "-3");
  } 
  else if (sim.rfind("hzagn") != std::string::npos)
  {
    boost::replace_all(sim, "hzagn", "HzAGN");
  }
  else if (sim.rfind("tng") != std::string::npos)
  {
    boost::replace_all(sim, "tng", "TNG");
  }
  
  std::string name;
  int tag;

  if (sim.rfind('-') != std::string::npos)
  {
    const size_t pos = sim.rfind('-');
    name = sim.substr(0, pos);
    tag = std::stoi(sim.substr(pos + 1));
  } 
  else
  { 
    name = sim;
    tag = 1; 
  }

  return std::make_tuple(name, tag);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void init_baryons_contamination(std::string sim)
{ // OLD API
  spdlog::debug("{}: Begins", "init_baryons_contamination");

  auto [name, tag] = get_baryon_sim_name_and_tag(sim);
  
  spdlog::debug("{}: Baryon simulation w/ Name = {} & Tag = {} selected",
    "init_baryons_contamination", name, tag);

  std::string tmp = name + "-" + std::to_string(tag);

  init_baryons(tmp.c_str());
  
  spdlog::debug("{}: Ends", "init_baryons_contamination");
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

#ifdef HDF5LIB
void init_baryons_contamination(std::string sim, std::string all_sims_hdf5_file)
{ // NEW API
  spdlog::debug("{}: Begins", "init_baryons_contamination");

  auto [name, tag] = get_baryon_sim_name_and_tag(sim);
       
  spdlog::debug("{}: Baryon simulation w/ Name = {} & Tag = {} selected",
    "init_baryons_contamination", name, tag);

  init_baryons_from_hdf5_file(name.c_str(), tag, all_sims_hdf5_file.c_str());

  spdlog::debug("{}: Ends", "init_baryons_contamination");
}
#endif

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

arma::Mat<double> compute_baryon_pcas_3x2pt_real(arma::Col<int>::fixed<3> order)
{
  const int ndata = IP::get_instance().get_ndata();

  const int ndata_sqzd = IP::get_instance().get_ndata_sqzd();
  
  const int nscenarios = BaryonScenario::get_instance().nscenarios();

  // Compute Cholesky Decomposition of the Covariance Matrix --------------
  spdlog::debug("{}: Computing Cholesky Decomposition of"
    " the Covariance Matrix begins", "compute_baryon_pcas_3x2pt_real");

  arma::Mat<double> L = arma::chol(IP::get_instance().get_cov_masked_sqzd(), "lower");

  arma::Mat<double> inv_L = arma::inv(L);

  spdlog::debug("{}: Computing Cholesky Decomposition of"
    " the Covariance Matrix ends", "compute_baryon_pcas_3x2pt_real");

  // Compute Dark Matter data vector --------------------------------------
  spdlog::debug("{}: Computing DM only data vector begins", 
    "compute_baryon_pcas_3x2pt_real");
  
  cosmology.random = RandomNumber::get_instance().get();
  
  reset_bary_struct(); // make sure there is no baryon contamination

  arma::Col<double> dv_dm = IP::get_instance().sqzd_theory_data_vector(
      compute_data_vector_3x2pt_real_masked_any_order(order)
    );

  spdlog::debug("{}: Computing DM only data vector ends", 
    "compute_baryon_pcas_3x2pt_real");

  // Compute data vector for all Baryon scenarios -------------------------
  arma::Mat<double> D = arma::Mat<double>(ndata_sqzd, nscenarios);

  for (int i=0; i<nscenarios; i++)
  {
    spdlog::debug("{}: Computing contaminated data vector"
      " with baryon scenario {} begins", "compute_baryon_pcas_3x2pt_real",
      BaryonScenario::get_instance().get_scenario(i));

    // line below to force clear cosmolike cosmology cache
    cosmology.random = RandomNumber::get_instance().get();

    init_baryons_contamination(BaryonScenario::get_instance().get_scenario(i));

    arma::Col<double> dv = IP::get_instance().sqzd_theory_data_vector(
        compute_data_vector_3x2pt_real_masked_any_order(order)
      );

    D.col(i) = dv - dv_dm;

    spdlog::debug("{}: Computing contaminated data vector"
      " with baryon scenario {} ends", "compute_baryon_pcas_3x2pt_real",
      BaryonScenario::get_instance().get_scenario(i));
  }

  reset_bary_struct();
  
  // line below to force clear cosmolike cosmology cache
  cosmology.random = RandomNumber::get_instance().get();

  // weight the diff matrix by inv_L; then SVD ----------------------------  
  arma::Mat<double> U, V;
  arma::Col<double> s;
  arma::svd(U, s, V, inv_L * D);

  // compute PCs ----------------------------------------------------------
  arma::Mat<double> PC = L * U; 

  // Expand the number of dims --------------------------------------------
  arma::Mat<double> R = arma::Mat<double>(ndata, nscenarios); 

  for (int i=0; i<nscenarios; i++)
    R.col(i)=IP::get_instance().expand_theory_data_vector_from_sqzd(PC.col(i));

  return R;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

arma::Mat<double> compute_baryon_pcas_6x2pt(arma::Col<int>::fixed<6> order)
{
  const int ndata = IP::get_instance().get_ndata();

  const int ndata_sqzd = IP::get_instance().get_ndata_sqzd();
  
  const int nscenarios = BaryonScenario::get_instance().nscenarios();

  // Compute Cholesky Decomposition of the Covariance Matrix --------------
  
  spdlog::debug("{}: Computing Cholesky Decomposition of"
    " the Covariance Matrix begins",  "compute_baryon_pcas_3x2pt");

  arma::Mat<double> L = arma::chol(IP::get_instance().get_cov_masked_sqzd(), "lower");

  arma::Mat<double> inv_L = arma::inv(L);

  spdlog::debug("{}: Computing Cholesky Decomposition of"
    " the Covariance Matrix ends", "compute_baryon_pcas_3x2pt");

  // Compute Dark Matter data vector --------------------------------------
  
  spdlog::debug("{}: Computing DM only data vector begins", 
    "compute_baryon_pcas_3x2pt");
  
  // line below to force clear cosmolike cosmology cache
  cosmology.random = RandomNumber::get_instance().get();
  
  reset_bary_struct(); // make sure there is no baryon contamination

  vector dv_dm = IP::get_instance().sqzd_theory_data_vector(
      compute_data_vector_6x2pt_real_masked_any_order(order)
    );

  spdlog::debug("{}: Computing DM only data vector ends", 
    "compute_baryon_pcas_3x2pt");

  // Compute data vector for all Baryon scenarios -------------------------
  
  arma::Mat<double> D = arma::Mat<double>(ndata_sqzd, nscenarios);

  for (int i=0; i<nscenarios; i++)
  {
    spdlog::debug("{}: Computing contaminated data vector"
      " with baryon scenario {} begins", "compute_baryon_pcas_3x2pt",
      BaryonScenario::get_instance().get_scenario(i));

    // line below to force clear cosmolike cosmology cache
    cosmology.random = RandomNumber::get_instance().get();

    init_baryons_contamination(BaryonScenario::get_instance().get_scenario(i));

    vector dv = IP::get_instance().sqzd_theory_data_vector(
        compute_data_vector_6x2pt_real_masked_any_order(order)
      );

    D.col(i) = dv - dv_dm;

    spdlog::debug("{}: Computing contaminated data vector"
      " with baryon scenario {} ends", "compute_baryon_pcas_3x2pt",
      BaryonScenario::get_instance().get_scenario(i));
  }

  reset_bary_struct();
  
  // line below to force clear cosmolike cosmology cache ------------------
  cosmology.random = RandomNumber::get_instance().get();
  
  // weight the diff matrix by inv_L; then SVD ----------------------------  
  arma::Mat<double> U, V;
  vector s;
  arma::svd(U, s, V, inv_L * D);
  // compute PCs ----------------------------------------------------------
  arma::Mat<double> PC = L * U; 
  // Expand the number of dims --------------------------------------------
  arma::Mat<double> R = arma::Mat<double>(ndata, nscenarios); 

  for (int i=0; i<nscenarios; i++)
    R.col(i) = IP::get_instance().expand_theory_data_vector_from_sqzd(PC.col(i));

  return R;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

bool BaryonScenario::is_pcs_set() const
{
  return this->is_pcs_set_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

bool BaryonScenario::is_scenarios_set() const
{
  return this->is_scenarios_set_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int BaryonScenario::nscenarios() const
{
  if (!this->is_scenarios_set_)
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "BaryonScenario::nscenarios", "Baryon Scenarios");
    exit(1);
  }
  return this->nscenarios_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

std::string BaryonScenario::get_scenario(const int i) const
{
  if (!this->is_scenarios_set_)
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "BaryonScenario::get_scenario", "Baryon Scenarios");
    exit(1);
  }
  return this->scenarios_.at(i);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void BaryonScenario::set_scenarios(std::string scenarios)
{
  std::vector<std::string> lines;
  lines.reserve(50);

  boost::trim_if(scenarios, boost::is_any_of("\t "));
  boost::trim_if(scenarios, boost::is_any_of("\n"));

  if (scenarios.empty())
  {
    spdlog::critical("{}: invalid string input (empty)",
      "BaryonScenario::set_scenarios");
    exit(1);
  }
  
  spdlog::debug("{}: Selecting baryon scenarios for PCA", 
    "BaryonScenario::set_scenarios");

  boost::split(
      lines, 
      scenarios, 
      boost::is_any_of("/ \t"), 
      boost::token_compress_on
    );
  
  int nscenarios = 0;
  
  for (auto it=lines.begin(); it != lines.end(); ++it)
  {
    auto [name, tag] = get_baryon_sim_name_and_tag(*it);

    this->scenarios_[nscenarios++] = name + "-" + std::to_string(tag);
  }

  this->nscenarios_ = nscenarios;

  spdlog::debug("{}: {} scenarios are registered", 
    "BaryonScenario::set_scenarios", this->nscenarios_);

  spdlog::debug("{}: Registering baryon scenarios for PCA done!", 
    "BaryonScenario::set_scenarios");
  
  this->is_scenarios_set_ = true;
  return;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void BaryonScenario::set_pcs(arma::Mat<double> eigenvectors)
{
  this->eigenvectors_ = eigenvectors;
  this->is_pcs_set_ = true;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

arma::Mat<double> BaryonScenario::get_pcs() const
{
  if (!this->is_pcs_set_)
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "BaryonScenario::get_pcs", "PC eigenvectors");
    exit(1);
  }
  return this->eigenvectors_;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double BaryonScenario::get_pcs(const int ci, const int cj) const
{
  if (!this->is_pcs_set_)
  {
    spdlog::critical("{}: {} not set prior to this function call",
      "BaryonScenario::get_pcs", "PC eigenvectors");
    exit(1);
  }
  return this->eigenvectors_(ci, cj); 
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
} // namespace cosmolike_interface
