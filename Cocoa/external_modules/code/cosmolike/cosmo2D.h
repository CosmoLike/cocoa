#ifndef __COSMOLIKE_COSMO2D_H
#define __COSMOLIKE_COSMO2D_H
#ifdef __cplusplus
extern "C" {
#endif

// Naming convention:
// g = galaxy positions ("g" as in "galaxy")
// k = kappa CMB ("k" as in "kappa")
// s = kappa from source galaxies ("s" as in "shear")

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

double C_ss_tomo_nointerp(double l, int ni, int nj);

double C_ss_tomo_withIA_nointerp(double s, int ni, int nj);

double C_gs_tomo_nointerp(double l, int ni, int nj);

double C_gs_tomo_withIA_nointerp(double s, int nl, int ns);

double C_gg_tomo_nointerp(double l, int ni, int nj);

double C_gk_tomo_nointerp(double l, int nl);

double C_ks_tomo_nointerp(double l, int ns);

double C_ks_tomo_withIA_nointerp(double s, int ni);

double C_kk_nointerp(double l);

// ----------------------------------------------------------------------------
// Angular Power Spectra - with interpolation
// ----------------------------------------------------------------------------

double C_ss_tomo(double l, int ni, int nj);

double C_ss_tomo_withIA(double l, int ni, int nj);

double C_gs_tomo(double l, int ni, int nj);

double C_gs_tomo_withIA(double l, int ni, int nj);

double C_gg_tomo(double l, int ni, int nj);

double C_gk_tomo(double l, int ni);

double C_ks_tomo(double l, int ni);

double C_kk(double l);

// ----------------------------------------------------------------------------
// Correlation Functions (real space) - flat sky
// ----------------------------------------------------------------------------

double xi_pm_tomo(int pm, double theta, int ni, int nj);

double w_gs_tomo(double theta, int ni, int nj);

double w_gg_tomo(double theta, int ni, int nj);

// ----------------------------------------------------------------------------
// Correlation Functions (real space) - full sky
// ----------------------------------------------------------------------------

double xi_pm_tomo_fullsky(int pm, int nt, int ni, int nj);

double w_gs_tomo_fullsky(int nt, int ni, int nj);

double w_gg_tomo_fullsky(int nt, int ni, int nj);

double w_gk_tomo_fullsky(int nt, int ni);

double w_ks_tomo_fullsky(int nt, int ni);

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Non-Limber (+RSD)
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Correlation Functions (real space)
// ----------------------------------------------------------------------------

double w_gs_fullsky_nonlimber(int nt, int ni, int nj);

double w_gg_tomo_fullsky_nonlimber(int nt, int ni, int nj);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD