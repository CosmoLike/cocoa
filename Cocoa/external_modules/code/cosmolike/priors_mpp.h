#ifndef __COSMOLIKE_PRIORS_MPP_H
#define __COSMOLIKE_PRIORS_MPP_H
#ifdef __cplusplus
extern "C" {
#endif

void set_shear_priors_mpp(double *mean_m, double *sigma_m);

void set_wlphotoz_priors_mpp(double *bias_photoz_s, double *sigma_photoz_s);

void set_clphotoz_priors_mpp(double *bias_photoz_l, double *sigma_photoz_l);

void set_ia_priors_mpp(double A_min, double A_max);

void set_b1_priors_mpp(double b1_min, double b1_max);

void set_b2_priors_mpp(double b2_min, double b2_max);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD