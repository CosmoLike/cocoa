#ifndef __COSMOLIKE_COSMO2D_FULLSKY_H
#define __COSMOLIKE_COSMO2D_FULLSKY_H
#ifdef __cplusplus
extern "C" {
#endif

// galaxy clustering tomography 2PCF galaxies in bins ni, nj,
// computed as sum P_l(cos(like.theta[nt])*C(l,ni,nj)
double w_tomo_fullsky(int nt, int ni, int nj);

// G-G lensing, lens bin ni, source bin nj, including IA contamination
// if like.IA = 3
double w_gamma_t_fullsky(int nt,int ni, int nj);

// shear tomography correlation functions, including IA contamination
// if like.IA = 3
double xi_pm_fullsky(int pm, int nt, int ni, int nj);

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// CMB x LSS Functions
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double w_gk_fullsky(int nt, int ni);

double w_ks_fullsky(int nt, int ni);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD