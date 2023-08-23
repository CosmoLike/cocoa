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
  if (cosmology.is_cached == 1)
  {
    return 0;
  } 
  else
  {
    return (C.Omega_m != cosmology.Omega_m || C.Omega_v != cosmology.Omega_v ||
            C.Omega_nu != cosmology.Omega_nu || C.h0 != cosmology.h0 ||
            C.MGSigma != cosmology.MGSigma || C.MGmu != cosmology.MGmu ||
            C.random != cosmology.random || C.Omega_b != cosmology.Omega_b ||
            C.sigma_8 != cosmology.sigma_8) ? 1 : 0;
  }
}

int recompute_zphot_shear(nuisancepara N)
{
  static int photoz = -1;
  if (photoz != redshift.shear_photoz)
  {
    //photoz = redshift.shear_photoz;
    photoz = (redshift.shear_photoz == 7) ? 4 : redshift.shear_photoz;
    return 1;
  }
  if (redshift.shear_photoz != 3 && redshift.shear_photoz != 4)
  {
    return 0;
  }
  int res = 0;
  for (int i=0; i<tomo.shear_Nbin; i++)
  {
    if (N.sigma_zphot_shear[i] != nuisance.sigma_zphot_shear[i] ||
        N.bias_zphot_shear[i]  != nuisance.bias_zphot_shear[i])
    {
      res = 1;
    }
  }
  return res;
}

int recompute_zphot_clustering(nuisancepara N) 
{
  static int photoz = -1;
  if (photoz != redshift.clustering_photoz) 
  {
    //photoz = redshift.clustering_photoz;
    photoz = (redshift.clustering_photoz == 7) ? 4 : redshift.clustering_photoz;
    return 1;
  }
  if (redshift.clustering_photoz != 3 && redshift.clustering_photoz != 4) 
  {
    return 0;
  }
  int res = 0;
  for (int i=0; i<tomo.clustering_Nbin; i++) 
  {
    if (N.sigma_zphot_clustering[i] != nuisance.sigma_zphot_clustering[i] ||
        N.bias_zphot_clustering[i]  != nuisance.bias_zphot_clustering[i]) 
    {
      res = 1;
    }
  }
  return res;
}

int recompute_IA(nuisancepara N) 
{
  if ((N.A_ia      != nuisance.A_ia)   || 
      (N.eta_ia    != nuisance.eta_ia) ||
      (N.A2_ia     != nuisance.A2_ia)  || 
      (N.eta_ia_tt != nuisance.eta_ia_tt))
  {
    return 1;
  }
  for (int i=0; i<tomo.shear_Nbin; i++)
  {
    if (N.A_z[i] != nuisance.A_z[i])
    {
      return 1;
    }
    if (N.A2_z[i] != nuisance.A2_z[i])
    {
      return 1;
    }
    if (N.b_ta_z[i] != nuisance.b_ta_z[i])
    {
      return 1;
    }
  }
  return 0;
}

int recompute_galaxies(galpara G, int i) 
{
  if(i < 0 || i > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input ni = %d", i);
    exit(1);
  }
  return ((G.b[i] != gbias.b[i])   || 
          (G.b2[i] != gbias.b2[i]) || 
          (G.bs2[i] != gbias.bs2[i])) ? 1 : 0;
}

int recompute_all_galaxies(galpara G) 
{
  for(int i=0; i<tomo.clustering_Nbin; i++) 
  {
    if(recompute_galaxies(G, i))
    {
      return 1;
    }
  }  
  return 0;
}

int recompute_yhalo(ynuisancepara N)
{
  if (N.gas_beta != ynuisance.gas_beta || N.gas_lgM0 != ynuisance.gas_lgM0 ||
      N.gas_eps1 != ynuisance.gas_eps1 || N.gas_eps2 != ynuisance.gas_eps2 || 
      N.gas_alpha != ynuisance.gas_alpha || N.gas_A_star != ynuisance.gas_A_star || 
      N.gas_lgM_star != ynuisance.gas_lgM_star || N.gas_sigma_star != ynuisance.gas_sigma_star || 
      N.gas_lgT_w != ynuisance.gas_lgT_w || N.gas_f_H != ynuisance.gas_f_H || 
      N.gas_Gamma_KS != ynuisance.gas_Gamma_KS ) 
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
    if (N.cluster_MOR[i] != nuisance.cluster_MOR[i]) 
    {
      return 1;
    }
  }
  for (int i=0; i<Cluster.N_SF; i++)
  {
    if (N.cluster_selection[i] != nuisance.cluster_selection[i]) 
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
  return (recompute_cosmo3D(C) || recompute_zphot_shear(N) || recompute_IA(N)) ? 1 : 0;
}

int recompute_gs(cosmopara C, galpara G, nuisancepara N) 
{
  if (recompute_shear(C, N) || recompute_zphot_clustering(N))
  {
    return 1;
  }
  for(int i=0; i<tomo.clustering_Nbin; i++) 
  {
    if(recompute_galaxies(G, i))
    {
      return 1;
    }
  }
  return 0;
}

int recompute_gg(cosmopara C, galpara G, nuisancepara N)
{
  if (recompute_cosmo3D(C) || recompute_zphot_clustering(N))
  {
    return 1;
  }
  for(int i=0; i<tomo.clustering_Nbin; i++) 
  {
    if(recompute_galaxies(G, i))
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

int recompute_gk(cosmopara C, galpara G, nuisancepara N)
{
  if (recompute_cosmo3D(C) || recompute_zphot_clustering(N))
  {
    return 1;
  }
  for(int i=0; i<tomo.clustering_Nbin; i++) 
  {
    if(recompute_galaxies(G, i))
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
  if (recompute_cosmo3D(C) || recompute_zphot_clustering(N) || recompute_yhalo(N2))
  {
    return 1;
  }
  for(int i=0; i<tomo.clustering_Nbin; i++) 
  {
    if(recompute_galaxies(G, i))
    {
      return 1;
    }
  }
  return 0;
}

int recompute_ys(cosmopara C, nuisancepara N, ynuisancepara N2)
{
  return recompute_shear(C, N) || recompute_yhalo(N2);
}

int recompute_yy(cosmopara C, ynuisancepara N)
{
  return recompute_cosmo3D(C) || recompute_yhalo(N);
}

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
  if(recompute_clusters(C, N) || recompute_zphot_clustering(N))
  {
    return 1;
  }
  for(int i=0; i<tomo.clustering_Nbin; i++) 
  {
    if(recompute_galaxies(G, i))
    {
      return 1;
    }
  }
  return 0;
}

int recompute_cs(cosmopara C, galpara G, nuisancepara N)
{
  if(recompute_clusters(C, N) || recompute_zphot_shear(N) || recompute_IA(N))
  {
    return 1;
  }
  for(int i=0; i<tomo.shear_Nbin; i++) 
  {
    if(recompute_galaxies(G, i))
    {
      return 1;
    }
  }
  return 0;
}