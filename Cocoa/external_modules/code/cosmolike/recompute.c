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
// UPDATE
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------

void update_cosmopara(cosmopara *C) 
{
  C->Omega_m = cosmology.Omega_m;
  C->Omega_v = cosmology.Omega_v;
  C->Omega_nu = cosmology.Omega_nu;
  C->h0 = cosmology.h0;
  C->MGSigma = cosmology.MGSigma;
  C->MGmu = cosmology.MGmu;
  C->random = cosmology.random;
}

void update_galpara(galpara *G) 
{
  for (int i=0; i<tomo.clustering_Nbin; i++) 
  {
    if (gbias.b[i] > 0.2 && gbias.b[i] < 20) 
    {
      G->b[i] = gbias.b[i];
      G->b2[i] = gbias.b2[i];
      G->bs2[i] = gbias.bs2[i];
      G->cg[i] = gbias.cg[i];
    } else 
    {
      printf("lens bin %d: neither HOD nor linear bias set, exit\n", i);
      exit(1);
    }
  }
}

void update_nuisance(nuisancepara *N)
{
  N->A_ia = nuisance.A_ia;
  N->beta_ia = nuisance.beta_ia;
  N->eta_ia = nuisance.eta_ia;
  N->eta_ia_highz = nuisance.eta_ia_highz;
  N->A2_ia = nuisance.A2_ia;
  N->eta_ia_tt = nuisance.eta_ia_tt;
  
  for (int i=0; i<tomo.clustering_Nbin; i++) 
  {
    N->sigma_zphot_clustering[i] = nuisance.sigma_zphot_clustering[i];
    N->bias_zphot_clustering[i] = nuisance.bias_zphot_clustering[i];
  }
  
  for (int i=0; i<tomo.shear_Nbin; i++) 
  {
    N->sigma_zphot_shear[i] = nuisance.sigma_zphot_shear[i];
    N->bias_zphot_shear[i] = nuisance.bias_zphot_shear[i];
    N->A_z[i] = nuisance.A_z[i];
    N->A2_z[i] = nuisance.A2_z[i];
    N->b_ta_z[i] = nuisance.b_ta_z[i];
  }
  
  for (int i=0; i<nuisance.N_cluster_MOR; ++i) 
  {
    N->cluster_MOR[i] = nuisance.cluster_MOR[i];
  }
  for (int i=0; i<nuisance.N_cluster_selection; ++i) 
  {
    N->cluster_selection[i] = nuisance.cluster_selection[i];
  }
}

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
            C.random != cosmology.random) ? 1 : 0;
  }
}

int recompute_zphot_shear(nuisancepara N)
{
  static int photoz = -1;
  if (photoz != redshift.shear_photoz)
  {
    photoz = redshift.shear_photoz;
    return 1;
  }
  if (redshift.shear_photoz != 3 && redshift.shear_photoz != 4)
  {
    return 0;
  }
  int res = 0;
  for (int i = 0; i<tomo.shear_Nbin; i++)
  {
    if (N.sigma_zphot_shear[i] != nuisance.sigma_zphot_shear[i] ||
        N.bias_zphot_shear[i] != nuisance.bias_zphot_shear[i])
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
    photoz = redshift.clustering_photoz;
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
        N.bias_zphot_clustering[i] != nuisance.bias_zphot_clustering[i]) 
    {
      res = 1;
    }
  }
  return res;
}

int recompute_IA(nuisancepara N) 
{
  if (N.A_ia != nuisance.A_ia || N.eta_ia != nuisance.eta_ia ||
      N.A2_ia != nuisance.A2_ia || N.eta_ia_tt != nuisance.eta_ia_tt)
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
  if(i < -1 || i > tomo.clustering_Nbin -1)
  {
    log_fatal("invalid bin input ni = %d", i);
    exit(1);
  }
  if (i == -1)
  {
    return 0;
  }
  return (G.b[i] != gbias.b[i] || G.b2[i] != gbias.b2[i] || G.bs2[i] != gbias.bs2[i] ||
    G.cg[i] != gbias.cg[i]) ? 1 : 0;
}

int recompute_clusters(cosmopara C, nuisancepara N)
{
  if (recompute_cosmo3D(C))
  {
    return 1;
  }
  for (int i=0; i<nuisance.N_cluster_MOR; ++i)
  {
    if (N.cluster_MOR[i] != nuisance.cluster_MOR[i]) 
    {
      return 1;
    }
  }
  for (int i=0; i<nuisance.N_cluster_selection; ++i)
  {
    if (N.cluster_selection[i] != nuisance.cluster_selection[i]) 
    {
      return 1;
    }
  }
  return 0;
}

int recompute_shear(cosmopara C, nuisancepara N)
{
  return (recompute_cosmo3D(C) || recompute_zphot_shear(N) || recompute_IA(N)) ? 1 : 0;
}

int recompute_gs(cosmopara C, galpara G, nuisancepara N) 
{
  if (recompute_cosmo3D(C) || recompute_zphot_clustering(N) || recompute_zphot_shear(N) || recompute_IA(N))
  {
    return 1;
  }
  const int nggl_size = tomo.ggl_Npowerspectra; 
  for(int i=0; i<nggl_size; i++) 
  {
    const int ZLNZ = ZL(i);
    const int ZSNZ = ZS(i);
    if(recompute_galaxies(G, ZLNZ) || recompute_galaxies(G, ZSNZ))
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
  const int ncl_size = tomo.clustering_Nbin; 
  for(int i=0; i<ncl_size; i++) 
  {
    if(recompute_galaxies(G, i))
    {
      return 1;
    }
  }
  return 0;
}

int recompute_ks(cosmopara C, galpara G, nuisancepara N)
{
  if (recompute_cosmo3D(C) || recompute_zphot_shear(N) || recompute_IA(N))
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