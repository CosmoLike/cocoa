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

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Correlation Functions (real space) - flat sky
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// nl = lambda_obs bin, ni = cluster redshift bin, nj = source redshift bin
double w_gammat_cluster_tomo_flatsky(double theta, int nl, int ni, int nj, int limber);

// nl{1,2} = lambda_obs bins, n{i,j} = cluster redshift bins
double w_cc_tomo_flatsky(double theta, int nl1, int nl2, int ni, int nj, int limber);

// nl = lambda_obs bin, ni = cluster redshift bin, nj = galaxy redshift bin
double w_cg_tomo_flatsky(double theta, int nl, int ni, int nj, int limber);

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Limber Approximation (Angular Power Spectrum)
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// nl = lambda_obs bin, ni = cluster redshift bin, nj = galaxy redshift bin
double C_cg_tomo_limber(double l, int nl, int ni, int nj);

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// cluster number counts
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// nl = lambda_obs bin, ni = cluster redshift bin
double projected_average_number_counts(int nl, int nz);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD