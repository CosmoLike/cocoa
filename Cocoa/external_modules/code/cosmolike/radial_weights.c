#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "log.c/src/log.h"

#include "basics.h"
#include "bias.h"
#include "cosmo2D_fourier.h"
#include "cosmo3D.h"
#include "redshift_spline.h"
#include "structs.h"

double W_kappa(double a, double fK, double nz) {
  double wkappa = 1.5 * cosmology.Omega_m * fK / a * g_tomo(a, (int)nz);
  if (cosmology.MGSigma != 0.) {
    wkappa *= (1. + MG_Sigma(a));
  }
  return wkappa;
}

double W2_kappa(double a, double fK, double nz) {
  const double tmp = (1.5*cosmology.Omega_m*fK/a);
  double wkappa = tmp*tmp*g2_tomo(a, (int) nz);
  if(cosmology.MGSigma != 0.){
    wkappa *= pow((1.+MG_Sigma(a)),2.);
  }
  return wkappa;
}

double W_mag(double a, double fK, double nz) {
  double wmag = 1.5 * cosmology.Omega_m * fK / a * g_lens(a, (int) nz);
  if (cosmology.MGSigma != 0.) {
    wmag *= (1. + MG_Sigma(a));
  }
  return wmag;
}

double W_gal(double a, double nz, double chi, double hoverh0) {
  double wgal = gbias.b1_function(1. / a - 1., (int) nz) *
    pf_photoz(1. / a - 1., (int) nz) * hoverh0;
  double wmag = gbias.b_mag[(int) nz] * 1.5 * cosmology.Omega_m * f_K(chi) /
    a * g_lens(a, (int) nz);
  if (cosmology.MGSigma != 0.) {
    wmag *= (1. + MG_Sigma(a));
  }
  return wgal + wmag;
}

double W_source(double a, double nz, double hoverh0) {
  return zdistr_photoz(1. / a - 1., (int)nz) * hoverh0;
}

double f_rsd(double aa) {
#ifdef DEBUG
  if(!(aa>0)) {
    log_fatal("Error on f_rsd: null scale factor");
    exit(1);
  }
#endif
  const double z = 1.0/aa-1.0;
  return f_growth(z);
}

double W_RSD(double l, double a0, double a1, double nz) {
  double w = (1 + 8. * l) / ((2. * l + 1.) * (2. * l + 1.)) *
             pf_photoz(1. / a0 - 1., (int) nz) * hoverh0(a0) * f_rsd(a0);
  w -= 4. / (2 * l + 3.) * sqrt((2 * l + 1.) / (2 * l + 3.)) *
       pf_photoz(1. / a1 - 1., (int) nz) * hoverh0(a1) * f_rsd(a1);
  return w;
}

double W_HOD(double a, double nz, double hoverh0) {
  return pf_photoz(1. / a - 1., (int) nz) * hoverh0;
}