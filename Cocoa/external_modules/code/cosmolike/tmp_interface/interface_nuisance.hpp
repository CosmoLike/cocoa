#include <carma.h>
#include <armadillo>

#ifndef __COSMOLIKE_GENERIC_INTERFACE_NUISANCE_HPP
#define __COSMOLIKE_GENERIC_INTERFACE_NUISANCE_HPP

void init_bias(arma::Col<double> bias_z_evol_model);

void init_IA(const int IA_MODEL, const int IA_REDSHIFT_EVOL);


void set_nuisance_bias(
    arma::Col<double> B1, 
    arma::Col<double> B2, 
    arma::Col<double> B_MAG
  );


void set_nuisance_IA(
    arma::Col<double> A1, 
    arma::Col<double> A2,
    arma::Col<double> BTA
  );

void set_nuisance_magnification_bias(
    arma::Col<double> B_MAG
  );

void set_nuisance_nonlinear_bias(
    arma::Col<double> B1,
    arma::Col<double> B2
  );

void set_nuisance_shear_calib(
    arma::Col<double> M
  );

} // namespace cosmolike_interface