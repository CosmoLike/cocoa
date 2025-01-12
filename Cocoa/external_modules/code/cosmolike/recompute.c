#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "basics.h"
#include "recompute.h"
#include "redshift_spline.h"
#include "structs.h"

#include "log.c/src/log.h"



// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// RECOMPUTE
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------

int recompute_cosmo3D(cosmopara C)
{
  // COCOA: Random number between (0,1). Optimization as the 
  return fdiff(C.random, cosmology.random) ? 1 : 0;
  /*
  return (fdiff(C.Omega_m, cosmology.Omega_m) ||
          fdiff(C.Omega_v, cosmology.Omega_v) ||
          fdiff(C.Omega_nu, cosmology.Omega_nu) ||
          fdiff(C.h0, cosmology.h0) ||
          fdiff(C.random, cosmology.random) ||
          fdiff(C.Omega_b, cosmology.Omega_b) ||
          fdiff(C.sigma_8, cosmology.sigma_8)) ? 1 : 0;
  */
}

int recompute_zphot_shear(photozparams N)
{
  int res = 0;
  for (int i=0; i<MAX_SIZE_ARRAYS; i++)
  {
    if (fdiff(N.sigma_zphot_shear[i], photoz.sigma_zphot_shear[i]) ||
        fdiff(N.bias_zphot_shear[i], photoz.bias_zphot_shear[i]))
    {
      res = 1;
    }
  }
  return res;
}

int recompute_zphot_clustering(photozparams N) 
{
  int res = 0;
  for (int i=0; i<MAX_SIZE_ARRAYS; i++) 
  {
    if (fdiff(N.sigma_zphot_clustering[i], photoz.sigma_zphot_clustering[i]) ||
        fdiff(N.bias_zphot_clustering[i], photoz.bias_zphot_clustering[i])) 
    {
      res = 1;
    }
  }
  return res;
}

int recompute_IA(nuisancepara N) 
{
  if (fdiff(N.A_ia, nuisance.A_ia) || 
      fdiff(N.eta_ia, nuisance.eta_ia) ||
      fdiff(N.A2_ia, nuisance.A2_ia) || 
      fdiff(N.eta_ia_tt, nuisance.eta_ia_tt))
  {
    return 1;
  }
  for (int i=0; i<MAX_SIZE_ARRAYS; i++)
  {
    if (fdiff(N.A_z[i], nuisance.A_z[i]) ||
        fdiff(N.A2_z[i], nuisance.A2_z[i]) ||
        fdiff(N.b_ta_z[i], nuisance.b_ta_z[i]))
    {
      return 1;
    }
  }
  return 0;
}

int recompute_galaxies(galpara G, int i) 
{
  if (i < 0 || i > MAX_SIZE_ARRAYS -1)
  {
    log_fatal("invalid bin input ni = %d", i);
    exit(1);
  }
  return (fdiff(G.b[i], gbias.b[i]) || 
          fdiff(G.b2[i], gbias.b2[i]) || 
          fdiff(G.bs2[i], gbias.bs2[i])) ? 1 : 0;
}

int recompute_all_galaxies(galpara G) 
{
  for (int i=0; i<MAX_SIZE_ARRAYS; i++) 
  {
    if (recompute_galaxies(G, i))
    {
      return 1;
    }
  }  
  return 0;
}

int recompute_yhalo(ynuisancepara N)
{
  if (fdiff(N.gas_beta, ynuisance.gas_beta) || 
      fdiff(N.gas_lgM0, ynuisance.gas_lgM0) ||
      fdiff(N.gas_eps1, ynuisance.gas_eps1) || 
      fdiff(N.gas_eps2, ynuisance.gas_eps2) || 
      fdiff(N.gas_alpha, ynuisance.gas_alpha) || 
      fdiff(N.gas_A_star, ynuisance.gas_A_star) || 
      fdiff(N.gas_lgM_star, ynuisance.gas_lgM_star) || 
      fdiff(N.gas_sigma_star, ynuisance.gas_sigma_star) || 
      fdiff(N.gas_lgT_w, ynuisance.gas_lgT_w) || 
      fdiff(N.gas_f_H, ynuisance.gas_f_H) || 
      fdiff(N.gas_Gamma_KS, ynuisance.gas_Gamma_KS)) 
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

int recompute_clusters(cosmopara C, nuisancepara N)
{
  if (recompute_cosmo3D(C))
  {
    return 1;
  }
  for (int i=0; i<Cluster.N_MOR; i++)
  {
    if (fdiff(N.cluster_MOR[i], nuisance.cluster_MOR[i])) 
    {
      return 1;
    }
  }
  for (int i=0; i<Cluster.N_SF; i++)
  {
    if (fdiff(N.cluster_selection[i], nuisance.cluster_selection[i]))
    {
      return 1;
    }
  }
  return 0;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int recompute_shear(cosmopara C, nuisancepara N)
{
  return (recompute_cosmo3D(C) || 
          recompute_zphot_shear(N) || 
          recompute_IA(N)) ? 1 : 0;
}

int recompute_gs(cosmopara C, galpara G, nuisancepara N) 
{
  if (recompute_shear(C, N) || recompute_zphot_clustering(N))
  {
    return 1;
  }
  for (int i=0; i<MAX_SIZE_ARRAYS; i++) 
  {
    if (recompute_galaxies(G, i))
    {
      return 1;
    }
  }
  return 0;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int recompute_gg(cosmopara C, galpara G, nuisancepara N)
{
  if (recompute_cosmo3D(C) || recompute_zphot_clustering(N))
  {
    return 1;
  }
  for (int i=0; i<MAX_SIZE_ARRAYS; i++) 
  {
    if (recompute_galaxies(G, i))
    {
      return 1;
    }
  }
  return 0;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int recompute_ks(cosmopara C, nuisancepara N)
{
  return recompute_shear(C, N);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int recompute_gk(cosmopara C, galpara G, nuisancepara N)
{
  if (recompute_cosmo3D(C) || recompute_zphot_clustering(N))
  {
    return 1;
  }
  for (int i=0; i<MAX_SIZE_ARRAYS; i++) 
  {
    if (recompute_galaxies(G, i))
    {
      return 1;
    }
  }
  return 0;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int recompute_gy(cosmopara C, galpara G, nuisancepara N, ynuisancepara N2) 
{
  if (recompute_cosmo3D(C) || 
      recompute_zphot_clustering(N) || 
      recompute_yhalo(N2))
  {
    return 1;
  }
  for (int i=0; i<MAX_SIZE_ARRAYS; i++) 
  {
    if (recompute_galaxies(G, i))
    {
      return 1;
    }
  }
  return 0;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int recompute_ys(cosmopara C, nuisancepara N, ynuisancepara N2)
{
  return recompute_shear(C, N) || recompute_yhalo(N2);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int recompute_yy(cosmopara C, ynuisancepara N)
{
  return recompute_cosmo3D(C) || recompute_yhalo(N);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int recompute_ky(cosmopara C, ynuisancepara N)
{
  return recompute_cosmo3D(C) || recompute_yhalo(N);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int recompute_cc(cosmopara C, nuisancepara N)
{
  return recompute_clusters(C, N);
}

int recompute_cg(cosmopara C, galpara G, nuisancepara N)
{
  if (recompute_clusters(C, N) || recompute_zphot_clustering(N))
  {
    return 1;
  }
  for (int i=0; i<MAX_SIZE_ARRAYS; i++) 
  {
    if (recompute_galaxies(G, i))
    {
      return 1;
    }
  }
  return 0;
}

int recompute_cs(cosmopara C, galpara G, nuisancepara N)
{
  if (recompute_clusters(C, N) || recompute_zphot_shear(N) || recompute_IA(N))
  {
    return 1;
  }
  for (int i=0; i<MAX_SIZE_ARRAYS; i++) 
  {
    if (recompute_galaxies(G, i))
    {
      return 1;
    }
  }
  return 0;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int recompute_table(Ntab N)
{
 // COCOA: "recompute" can be too expensive if required to  
 // COCOA: check every element of the struct. Whenever 
 // COCOA: updating the struct, assign a new random number
  return fdiff(N.random, Ntable.random) ? 1 : 0;
  /*int res = 0;
  if ((N.N_a != Ntable.N_a) || 
      (N.N_ell != Ntable.N_ell) ||
      (N.Ntheta != Ntable.Ntheta) || 
      (N.N_ell_TATT != Ntable.N_ell_TATT) ||
      (N.high_def_integration != Ntable.high_def_integration)
     )
  {
    res = 1;
  }
  return res;
  */
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int recompute_redshift(redshiftparams N)
{
  // COCOA: Random number between (0,1). Optimization as the 
  // COCOA: function "recompute_redshift" is too expensive
  // COCOA: if required to check every element of the 2d arrays
  // COCOA: Whenever updating n(z), just create a new random num
  return fdiff(N.random, redshift.random) ? 1 : 0;

  /*
  if (redshift.clustering_zdist_table == NULL || 
      redshift.shear_zdist_table == NULL)
  {
    log_fatal("error: redshift struct not set");
    exit(1);
  }

  if ((N.shear_nbin != redshift.shear_nbin) ||
      (N.shear_nzbins != redshift.shear_nzbins) ||
      fdiff(N.shear_zdist_zmin_all, redshift.shear_zdist_zmin_all) ||
      fdiff(N.shear_zdist_zmax_all, redshift.shear_zdist_zmax_all) ||
      (N.clustering_nbin != redshift.clustering_nbin) ||
      (N.clustering_nzbins != redshift.clustering_nzbins) ||
      fdiff(N.clustering_zdist_zmin_all, redshift.clustering_zdist_zmin_all) ||
      fdiff(N.clustering_zdist_zmax_all, redshift.clustering_zdist_zmax_all))
  {
    return 1;
  }

  for(int i=0; i<MAX_SIZE_ARRAYS; i++)
  {
    if (fdiff(N.shear_zdist_zmin[i], redshift.shear_zdist_zmin[i]) ||
        fdiff(N.shear_zdist_zmax[i], redshift.shear_zdist_zmax[i]) ||
        fdiff(N.clustering_zdist_zmin[i], redshift.clustering_zdist_zmin[i]) ||
        fdiff(N.clustering_zdist_zmax[i], redshift.clustering_zdist_zmax[i]))
    {
      return 1;
    }

  }

  if (N.clustering_zdist_table != NULL)
  {
    for (int k=0; k<redshift.clustering_nzbins; k++) 
    { 
      for (int i=0; i<redshift.clustering_nbin+1; i++) 
      {
        if (fdiff(N.clustering_zdist_table[i][k], redshift.clustering_zdist_table[i][k]))
        {
          return 1;
        }
      }
    }
  }
  else
  {
    return 1;
  } 
  if (N.shear_zdist_table != NULL)
  {
    for (int k=0; k<redshift.shear_nzbins; k++) 
    { 
      for (int i=0; i<redshift.shear_nbin+1; i++) 
      {
        if (fdiff(N.shear_zdist_table[i][k], redshift.shear_zdist_table[i][k]))
        {
          return 1;
        }
      }
    }
  }
  else
  {
    return 1;
  }
  
  return 0;
  */
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------