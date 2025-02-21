#ifndef __COSMOLIKE_HALO_H
#define __COSMOLIKE_HALO_H
#ifdef __cplusplus
extern "C" {
#endif

// HALO BIAS OPTIONS ---------------------------
#define HALO_BIAS_TINKER_2010 0

// HMF OPTIONS ---------------------------------
#define HMF_TINKER_2010 0

// CONCENTRATION OPTIONS -----------------------
#define CONCENTRATION_BHATTACHARYA_2013 0

// HALO PROFILE OPTIONS OPTIONS -----------------------
#define HALO_PROFILE_NFW 0

double p_mm(const double k, const double a);

double p_gm(const double k, const double a, const int ni);

double p_gg(const double k, const double a, const int ni, const int nj);

double p_my(const double k, const double a);

double p_yy(const double k, const double a);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD