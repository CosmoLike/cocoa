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

// ss in real space has a special name
double xi_pm_tomo(
    const int pm, 
    const int nt, 
    const int ni, 
    const int nj, 
    const int limber
  );

// gs in real space has a special name
double w_gammat_tomo(int nt, int ni, int nj, int limber);

double w_gg_tomo(int nt, int ni, int nj, int limber);

double w_gk_tomo(int nt, int ni, int limber);

double w_ks_tomo(int nt, int ni, int limber);

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Non Limber (Angular Power Spectrum)
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// Mixture of non-Limber and Limber of C_cl (G-G lensing)
void C_gl_tomo(int L, int nl, int ns, double *Cl, double dev, double tolerance);

// Mixture of non-Limber and Limber of C_cl (galaxy clustering)
void C_cl_tomo(int L, int ni, int nj, double *Cl, double dev, double tolerance);

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Limber Approximation (Angular Power Spectrum)
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double C_ss_TATT_EE_tomo_limber_nointerp(double l, int ni, int nj, 
  const int init_static_vars_only); // IA=(5||6)

double C_ss_TATT_EE_tomo_limber(double l, int ni, int nj, 
  const int force_no_recompute); // IA=(5||6)


double C_ss_TATT_BB_tomo_limber_nointerp(double l, int ni, int nj, 
  const int init_static_vars_only); // IA=(5||6)

double C_ss_TATT_BB_tomo_limber(double l, int ni, int nj, 
  const int force_no_recompute); // IA=(5||6)


double C_ss_tomo_limber_nointerp(double l, int ni, int nj, 
  int use_linear_ps, const int init_static_vars_only);

double C_ss_tomo_limber(double l, int ni, int nj, 
  const int force_no_recompute); // IA=(0||3||4)


// works with IA=(0||3||4||5||6)
double C_gs_tomo_limber_nointerp(double l, int ni, int nj, 
  int use_linear_ps, const int init_static_vars_only);

double C_gs_tomo_limber(double l, int ni, int nj, 
  const int force_no_recompute);


// ----------------------------------------------------------------------------
// All functions below can run w/ like.IA=0 || like.IA=3 || like.IA=4
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Interpolated Version (On first call, Func will create an interpolation table)
// ----------------------------------------------------------------------------

double C_gg_tomo_limber(double l, int ni, int nj, const int force_no_recompute);

double C_ks_tomo_limber(double l, int ni, const int force_no_recompute);

double C_kk_limber(double l, const int force_no_recompute);

double C_gy_tomo_limber(double l, int ni, const int force_no_recompute);

double C_ys_tomo_limber(double l, int ni, const int force_no_recompute);

double C_ky_limber(double l, const int force_no_recompute);

double C_yy_limber(double l, const int force_no_recompute);

double C_gk_tomo_limber(double l, int ni, const int force_no_recompute);

// ----------------------------------------------------------------------------
// Non-Interpolated Version (Will compute the Integral at every call)
// ----------------------------------------------------------------------------

// WARNING: C_gg beyond linear bias for cross-tomography bins not yet supported
double C_gg_tomo_limber_nointerp(double l, int ni, int nj, 
  int use_linear_ps, const int init_static_vars_only);

double C_gk_tomo_limber_nointerp(double l, int nl, 
  int use_linear_ps, const int init_static_vars_only);

double C_ks_tomo_limber_nointerp(double l, int ns, 
  int use_linear_ps, const int init_static_vars_only);

double C_kk_limber_nointerp(double l, int use_linear_ps, 
  const int init_static_vars_only);

double C_gy_tomo_limber_nointerp(double l, int ni, 
  int use_linear_ps,  const int init_static_vars_only);

double C_ys_tomo_limber_nointerp(double l, int ni, 
  int use_linear_ps, const int init_static_vars_only);

double C_ky_limber_nointerp(double l, 
  int use_linear_ps, const int init_static_vars_only);

double C_yy_limber_nointerp(double l, 
  int use_linear_ps, const int init_static_vars_only);

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Integrand for Angular Power Spectrum - Limber Approximation 
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double int_for_C_ss_TATT_EE_tomo_limber(double a, void* params);

double int_for_C_ss_TATT_BB_tomo_limber(double a, void* params);

double int_for_C_ss_tomo_limber(double a, void* params);

double int_for_C_gs_TATT_tomo_limber(double a, void* params);

double int_for_C_gs_tomo_limber(double a, void* params);

double int_for_C_gs_tomo_limber_withb2(double a, void* params);

double int_for_C_gg_tomo_limber(double a, void* params);

double int_for_C_gg_tomo_limber_withb2(double a, void* params);

double int_for_C_gk_tomo_limber(double a, void* params);

double int_for_C_gk_tomo_limber_withb2(double a, void* params);

double int_for_C_ks_tomo_limber(double a, void* params);

double int_for_C_kk_limber(double a, void* params);

double int_for_C_gy_tomo_limber(double a, void* params);

double int_for_C_ys_tomo_limber(double a, void* params);

double int_for_C_ky_limber(double a, void* params);

double int_for_C_yy_limber(double a, void *params);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD
