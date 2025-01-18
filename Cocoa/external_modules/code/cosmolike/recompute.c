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

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------