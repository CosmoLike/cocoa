#ifndef __COSMOLIKE_COSMO2D_H
#define __COSMOLIKE_COSMO2D_H
#ifdef __cplusplus
extern "C" {
#endif

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Limber Approximation
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Angular Power Spectra - no interpolation
// ----------------------------------------------------------------------------

//shear tomography power spectra of source galaxy bins ni, nj
double C_shear_tomo_nointerp(double l, int ni, int nj);

//G-G lensing power spectrum, lens bin ni, source bin nj
double C_gl_tomo_nointerp(double l, int ni, int nj);

// galaxy clustering power spectrum bins ni, nj (w/o look-up tables)
double C_cl_tomo_nointerp(double l, int ni, int nj);

// galaxy position x kappa CMB, lens z-bin nl
double C_gk_tomo_nointerp(double l, int nl);

// shear x kappa CMB, for source z-bin ns
double C_ks_tomo_nointerp(double l, int ns);

// kappa CMB x kappa CMB
double C_kk_nointerp(double l);

// ----------------------------------------------------------------------------
// Angular Power Spectra - with interpolation
// ----------------------------------------------------------------------------

double C_shear_tomo(double l, int ni, int nj);

double C_gl_tomo(double l, int ni, int nj);

double C_cl_tomo(double l, int ni, int nj);

double C_gk_tomo(double l, int ni);

double C_ks_tomo(double l, int ni);

double C_kk(double l);

// ----------------------------------------------------------------------------
// Correlation Functions (real space) - flat sky
// ----------------------------------------------------------------------------

double xi_pm_tomo(int pm, double theta, int ni, int nj);

double w_gamma_t_tomo(double theta,int ni, int nj);

double w_tomo(double theta, int ni, int nj);

// ----------------------------------------------------------------------------
// Correlation Functions (real space) - full sky
// ----------------------------------------------------------------------------

double w_tomo_fullsky(int nt, int ni, int nj);

double w_gamma_t_tomo_fullsky(int nt, int ni, int nj);

double xi_pm_tomo_fullsky(int pm, int nt, int ni, int nj);

double w_gk_tomo_fullsky(int nt, int ni);

double w_ks_tomo_fullsky(int nt, int ni);

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Non-Limber
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Correlation Functions (real space)
// ----------------------------------------------------------------------------

double w_gamma_t_fullsky_nonLimber(int nt, int ni, int nj);

// w(theta) including non-Limber+RSD
double w_tomo_fullsky_nonLimber(int nt, int ni, int nj);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD