#include <carma.h>
#include <armadillo>

#ifndef __COSMOLIKE_GENERIC_INTERFACE_COSMOLOGY_HPP
#define __COSMOLIKE_GENERIC_INTERFACE_COSMOLOGY_HPP

namespace cosmolike_interface
{

void set_cosmological_parameters(
    const double omega_matter,
    const double hubble
  );

void set_distances(
    arma::Col<double> io_z, 
    arma::Col<double> io_chi
  );

void set_growth(
    arma::Col<double> io_z, 
    arma::Col<double> io_G
  );

void set_linear_power_spectrum(
    arma::Col<double> io_log10k,
    arma::Col<double> io_z, 
    arma::Col<double> io_lnP
  );

void set_non_linear_power_spectrum(
    arma::Col<double> io_log10k,
    arma::Col<double> io_z, 
    arma::Col<double> io_lnP
  );

} // namespace cosmolike_interface