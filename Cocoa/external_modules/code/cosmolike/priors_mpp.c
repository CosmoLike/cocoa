#include "priors_mpp.h"
#include "structs.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void set_b1_priors_mpp(double b1_min __attribute__((unused)), double b1_max __attribute__((unused))) {
  printf("priors_mpp.c: set_b1_priors_mpp: flat bias priors should be applied "
         "in python script\n");
}

void set_b2_priors_mpp(double b2_min __attribute__((unused)), double b2_max __attribute__((unused))) {
  printf("priors_mpp.c: set_b2_priors_mp: flat bias priors should be applied "
         "in python script\n");
}

void set_ia_priors_mpp(double A_min __attribute__((unused)), double A_max __attribute__((unused))) {
  printf("priors_mpp.c: set_ia_priors_mpp: flat A priors should be applied in "
         "python script\n");
}

void set_shear_priors_mpp(double *mean_m, double *sigma_m) {
  printf("Setting Gaussian Priors on shear calibration\n");
  for (int i = 0; i < 10; i++) {
    prior.shear_calibration_m[i][0] = mean_m[i];
    prior.shear_calibration_m[i][1] = sigma_m[i];
    if (sigma_m[i]) {
      printf("zn = %d (mean,sigma) = (%e, %e)\n", i,
             prior.shear_calibration_m[i][0], prior.shear_calibration_m[i][1]);
      like.shearcalib = 1;
    }
  }
  printf("\n");
}

void set_wlphotoz_priors_mpp(double *bias_photoz_s, double *sigma_b_photoz_s) {
  printf("Setting Gaussian Priors on source n(z) bias\n");
  for (int i = 0; i < 10; i++) {
    prior.bias_zphot_shear[i][0] = bias_photoz_s[i];
    prior.bias_zphot_shear[i][1] = sigma_b_photoz_s[i];
    if (sigma_b_photoz_s[i]) {
      printf("zn = %d (mean,sigma) = (%e, %e)\n", i,
             prior.bias_zphot_shear[i][0], prior.bias_zphot_shear[i][1]);
      like.wlphotoz = 1;
    }
  }
  printf("\n");
}

void set_clphotoz_priors_mpp(double *bias_photoz_l, double *sigma_b_photoz_l) {
  printf("Setting Gaussian Priors on lens n(z) bias\n");
  for (int i = 0; i < 10; i++) {
    prior.bias_zphot_clustering[i][0] = bias_photoz_l[i];
    prior.bias_zphot_clustering[i][1] = sigma_b_photoz_l[i];
    if (sigma_b_photoz_l[i]) {
      printf("zn = %d (mean,sigma) = (%e, %e)\n", i,
             prior.bias_zphot_clustering[i][0],
             prior.bias_zphot_clustering[i][1]);
      like.clphotoz = 1;
    }
  }
  printf("\n");
}
