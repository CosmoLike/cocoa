#include <carma.h>
#include <armadillo>

#ifndef __COSMOLIKE_GENERIC_INTERFACE_REDSHIFT_HPP
#define __COSMOLIKE_GENERIC_INTERFACE_REDSHIFT_HPP

namespace cosmolike_interface
{

void init_redshift_distributions_from_files(
    std::string lens_multihisto_file, 
    const int lens_ntomo,
    std::string source_multihisto_file, 
    const int source_ntomo
  );

void set_nuisance_clustering_photoz(arma::Col<double> CP);

void set_nuisance_shear_photoz(arma::Col<double> SP);

} // namespace cosmolike_interface