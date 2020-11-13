#include <assert.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "basics.h"
#include "bias.h"
#include "cosmo3D.h"
#include "redshift_spline.h"
#include "structs.h"


// ---------------------------------------------------------------------------
// galaxy bias & HOD routines
// ---------------------------------------------------------------------------

// fitting formula from Lazeyras et al. 2016 (Eq. 5.2)
// https://arxiv.org/abs/1511.01096
double b2_from_b1(double b1) {
  return 0.412 - 2.143 * b1 + 0.929 * b1 * b1 + 0.008 * b1 * b1 * b1;
}

// e.g. https://arxiv.org/pdf/1405.1447v4.pdf
double bs2_from_b1(double b1) {
  return -4. / 7. * (b1 - 1.0);
}

double b3nl_from_b1(double b1) {
  return (b1 - 1.0);
}

double b1_per_bin_evolv(double z, int ni) {
  return gbias.b[ni] * pow((1. + z) / (1. + zmean(ni)), 1.0);
}

double b1_per_bin(double z __attribute__((unused)), int ni) {
  return gbias.b[ni];
}

double b1_per_bin_pass_evolv(double z, int ni) {
  double z_evolv_passiv = growfac(1. / (z + 1.)) /
                          growfac(1. / (1. + 0.5 * (tomo.clustering_zmin[ni] +
                                                    tomo.clustering_zmax[ni])));
  return (gbias.b[ni] - 1.) / z_evolv_passiv + 1.;
}

double b1_growth_scaling(double z, int ni __attribute__((unused))) {
  return gbias.b[0] / (growfac(1. / (z + 1.)) / growfac(1.));
}

double b1_powerlaw(double z, int ni __attribute__((unused))) {
  return gbias.b[0] * pow(1 + z, gbias.b[1]);
}

double bgal_z(double z, int ni) {
  // bias evolution within redshift bin, used by
  // clustering/G-G-lensing routines without HOD modeling
  // change this into desired redshift evolution as function of (z, z_pivot =
  // gbias[ni][1]), with z_evolv(z_pivot) =1 e.g. z_evolv =
  // pow((1+z)/(1+gbias[ni][1]),0.5); passive evolution  following Fry 1996:
  // b(z) = (b(z_pivot)-1)*D(z_pivot)/D(z) + 1

  // no tomography ->  bias = 1
  if (ni == -1) {
    return 1.0;
  }

  // -------------------------------------------------------------
  // new default: use b1_function specified in init -
  // if different from bgal_z, to avoid circular call
  // -------------------------------------------------------------
  if ((gbias.b1_function) && (gbias.b1_function != &bgal_z)) {
    return gbias.b1_function(z, ni);
  }

  // -------------------------------------------------------------
  // old options, for backward compatability
  // -------------------------------------------------------------
  // very broad redshift distribution -> no bias evolution model
  if (redshift.clustering_photoz == 4) {
    return gbias.b[ni];
  }
  else { // old default: passive evolution
    double z_evolv_passiv =
        growfac(1. / (z + 1.)) /
        growfac(1. / (1. + 0.5 * (tomo.clustering_zmin[ni] +
                                  tomo.clustering_zmax[ni])));
    return (gbias.b[ni] - 1.) / z_evolv_passiv + 1.;
  }
}