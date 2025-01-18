#ifndef __COSMOLIKE_COSMO2D_H
#define __COSMOLIKE_COSMO2D_H
#ifdef __cplusplus
extern "C" {
#endif

// ----------------------------------------------------------------------------
// Naming convention:
// ----------------------------------------------------------------------------
// c = cluster position ("c" as in "cluster")
// g = galaxy positions ("g" as in "galaxy")
// k = kappa CMB ("k" as in "kappa")
// s = kappa from source galaxies ("s" as in "shear")

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Correlation Functions (real Space) - Full Sky - bin average
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// ss in real space
double xi_pm_tomo(
    const int pm, const int nt, const int ni, const int nj, const int limber);

// gs in real space
double w_gammat_tomo(const int nt, const int ni, const int nj, const int limber);

double w_gg_tomo(const int nt, const int ni, const int nj, const int limber);

double w_gk_tomo(const int nt, const int ni, const int limber);

double w_ks_tomo(const int nt, const int ni, const int limber);

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Limber Approximation (Angular Power Spectrum)
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double C_ss_tomo_limber(const double l, const int ni, const int nj, const int EE);

double C_gs_tomo_limber(const double l, const int ni, const int nj);

double C_gg_tomo_limber(const double l, const int ni, const int nj);

double C_ks_tomo_limber(const double l, const int ni);

double C_gk_tomo_limber(const double l, const int ni);

double C_kk_limber(const double l);

double C_gy_tomo_limber(const double l, const int ni);

double C_ys_tomo_limber(const double l, const int ni);

double C_ky_limber(const double l);

double C_yy_limber(const double l);

// ----------------------------------------------------------------------------
// Non-Interpolated Version (Will compute the Integral at every call)
// ----------------------------------------------------------------------------

double C_ss_tomo_limber_nointerp(double l, int ni, int nj, int EE, const int init);

double C_gs_tomo_limber_nointerp(const double l, const int ni, const int nj, const int init);

double C_gg_tomo_limber_nointerp(double l, int ni, int nj, const int init);

double C_gg_tomo_limber_nointerp(double l, int ni, int nj, const int init);

double C_gk_tomo_limber_nointerp(double l, int nl, const int init);

double C_ks_tomo_limber_nointerp(double l, int ns, const int init);

double C_kk_limber_nointerp(double l, const int init);

double C_gy_tomo_limber_nointerp(double l, int ni, const int init);

double C_ys_tomo_limber_nointerp(double l, int ni, const int init);

double C_ky_limber_nointerp(const double l, const int init);

double C_yy_limber_nointerp(const double l, const int init);

// ----------------------------------------------------------------------------
// Integrands 
// ----------------------------------------------------------------------------

double int_for_C_ss_tomo_limber(double a, void* params);

double int_for_C_gs_tomo_limber(double a, void* params);

double int_for_C_gg_tomo_limber(double a, void* params);

double int_for_C_gk_tomo_limber(double a, void* params);

double int_for_C_ks_tomo_limber(double a, void* params);

double int_for_C_kk_limber(double a, void* params);

double int_for_C_gy_tomo_limber(double a, void* params);

double int_for_C_ys_tomo_limber(double a, void* params);

double int_for_C_ky_limber(double a, void* params);

double int_for_C_yy_limber(double a, void *params);

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Non Limber (Angular Power Spectrum)
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

void C_cl_tomo(int L, int ni, int nj, double *Cl, double dev, double tolerance);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD
