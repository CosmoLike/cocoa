#ifndef __COSMOLIKE_COSMO2D_CLUSTERH
#define __COSMOLIKE_COSMO2D_CLUSTERH
#ifdef __cplusplus
extern "C" {
#endif

// ----------------------------------------------------------------------------
// Naming convention: (same as cosmo2D.h)
// ----------------------------------------------------------------------------
// c = cluster position ("c" as in "cluster")
// g = galaxy positions ("g" as in "galaxy")
// k = kappa CMB ("k" as in "kappa")
// s = kappa from source galaxies ("s" as in "shear")

// ----------------------------------------------------------------------------
// Threading
// ----------------------------------------------------------------------------
// Thread loops in lambda_obs is not allowed. Most functions update static arrays 
// when varying lambda_obs in cluster_utils. With lambda_obs fixed, loops on 
// redshift bins can be threaded using the standard loop unrolling technique

// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// Correlation Functions (real Space) - Full Sky - bin average
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// nt = theta bin

// nl = lambda_obs bin, ni = cluster redshift bin, nj = source redshift bin
double w_gammat_cluster_tomo(const int nt, const int nl, const int ni, const int nj, 
const int limber);

// nl{1,2} = lambda_obs bins, n{i,j} = cluster redshift bins
double w_cc_tomo(const int nt, const int nl1, const int nl2, const int ni, const int nj, 
const int limber);

// nl = lambda_obs bin, ni = cluster redshift bin, nj = galaxy redshift bin
double w_cg_tomo(const int nt, const int nl, const int ni, const int nj, const int limber);

// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// Correlation Functions (real space) - flat sky
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------

// nl = lambda_obs bin, ni = cluster redshift bin, nj = source redshift bin
double w_gammat_cluster_tomo_flatsky(const double theta, const int nl, const int ni, 
const int nj, const int limber);

// nl{1,2} = lambda_obs bins, n{i,j} = cluster redshift bins
double w_cc_tomo_flatsky(const double theta, const int nl1, const int nl2, const int ni, 
const int nj, const int limber);

// nl = lambda_obs bin, ni = cluster redshift bin, nj = galaxy redshift bin
double w_cg_tomo_flatsky(double theta, int nl, int ni, int nj, int limber);

// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// Non Limber (Angular Power Spectrum)
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
void f_chi_for_Psi_cluster_cl(double *const chi, const int Nchi, double *const fchi, const int ni,
const int nl, const double zmin, const double zmax);

void f_chi_for_Psi_cluster_cl_RSD(double *const chi, int Nchi, double *const fchi, int ni, int nl,
const double zmin, const double zmax);

void f_chi_for_Psi_cluster_cl_Mag(double *const chi, const int Nchi, double *const fchi,
const int ni, const int nl, const double zmax);

// nl = lambda_obs bin, ni = cluster redshift bin, nj = galaxy redshift bin
void C_cg_tomo(int L, const int nl, const int ni, const int nj, double *const Cl, double dev,
const double tol);

// nl{1,2} = lambda_obs bins, n{i,j} = cluster redshift bins
void C_cc_tomo(int L, const int nl1, const int nl2, const int ni, const int nj, double *const Cl,
double dev, const double tol);

// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// Limber Approximation (Angular Power Spectrum)
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------

// nl = lambda_obs bin, ni = cluster redshift bin, nj = source redshift bin
double C_cs_tomo_limber_nointerp(const double l, const int nl, const int ni, const int nj, 
const int use_linear_ps, const int init_static_vars_only);

double C_cs_tomo_limber(const double l, const int nl, const int ni, const int nj);

// nl{1,2} = lambda_obs bins, n{i,j} = cluster redshift bins
double C_cc_tomo_limber_nointerp(const double l, const int nl1, const int nl2, const int ni, 
const int nj, const int use_linear_ps, const int init_static_vars_only);

double C_cc_tomo_limber(const double l, const int nl1, const int nl2, const int ni, 
const int nj);

// nl = lambda_obs bin, ni = cluster redshift bin, nj = galaxy redshift bin
double C_cg_tomo_limber_nointerp(const double l, const int nl, const int ni, const int nj, 
const int use_linear_ps, const int init_static_vars_only);

double C_cg_tomo_limber(const double l, const int nl, const int ni, const int nj);

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// cluster number counts
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// nl = lambda_obs bin

double binned_N_nointerp(const int nl, const int nz, const int interpolate_survey_area, 
const int init_static_vars_only);

double binned_N(const int nl, const int nz);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD