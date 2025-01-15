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

/*
int recompute_clusters(cosmopara C, nuisanceparams N)
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
*/

/*
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int recompute_cc(cosmopara C, nuisancepara N)
{
  return recompute_clusters(C, N);
}

int recompute_cg(cosmopara C, galpara G, nuisancepara N)
{
  if (recompute_clusters(C, N))
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
  if (recompute_clusters(C, N) || recompute_IA(N))
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
*/

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------