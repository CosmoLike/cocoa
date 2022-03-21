#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "log.c/src/log.h"

#include "basics.h"
#include "bias.h"
#include "cosmo3D.h"
#include "redshift_spline.h"
#include "structs.h"

double W_kappa(double a, double fK, int nz) 
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  if(nz < -1 || nz > tomo.shear_Nbin - 1) 
  {
    log_fatal("invalid bin input ni = %d", nz);
    exit(1);
  }
  double wkappa = (1.5*cosmology.Omega_m*fK/a) * g_tomo(a, nz);
  if (cosmology.MGSigma != 0.) 
  {
    wkappa *= (1. + MG_Sigma(a));
  }
  return wkappa;
}

double W2_kappa(double a, double fK, int nz) 
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  if(nz < -1 || nz > tomo.shear_Nbin - 1) 
  {
    log_fatal("invalid bin input ni = %d", nz);
    exit(1);
  }
  const double tmp = (1.5*cosmology.Omega_m*fK/a);
  double wkappa = tmp*tmp*g2_tomo(a, nz);
  if(cosmology.MGSigma != 0)
  {
    const double MG = (1.0 + MG_Sigma(a));
    wkappa *= MG*MG;
  }
  return wkappa;
}

double W_mag(double a, double fK, int nz) 
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  if(nz < -1 || nz > tomo.clustering_Nbin - 1) 
  {
    log_fatal("invalid bin input ni = %d (max %d)", nz, tomo.clustering_Nbin);
    exit(1);
  }
  double wmag = (1.5 * cosmology.Omega_m * fK / a) * g_lens(a, nz);
  if (cosmology.MGSigma != 0) 
  {
    wmag *= (1. + MG_Sigma(a));
  }
  return wmag;
}

/*
double W_mag_cluster(double a, double fK, int nz, int nl)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  if(nz < -1 || nz > tomo.cluster_Nbin - 1) 
  {
    log_fatal("invalid bin input ni = %d (max %d)", nz, tomo.cluster_Nbin);
    exit(1);
  }
  if(nl < 0 || nl > Cluster.N200_Nbin - 1) 
  {
    log_fatal("invalid bin input ni = %d (max %d)", nl, Cluster.N200_Nbin);
    exit(1);
  } 
  
  double wmag = (1.5 * cosmology.Omega_m * fK / a)  * g_lens_cluster(a, nz, nl);
  if(cosmology.MGSigma != 0.)
  {
    wmag *= (1.0 + MG_Sigma(a));
  }
  return wmag;
}
*/

double W_gal(double a, int nz, double chi, double hoverh0) 
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  if(nz < -1 || nz > tomo.clustering_Nbin - 1) 
  {
    log_fatal("invalid bin input ni = %d", nz);
    exit(1);
  }
  const double z = 1. / a - 1;
  const double fK = f_K(chi);
  const double wgal = gbias.b1_function(z, nz) * pf_photoz(z, nz) * hoverh0;
  return wgal + gbias.b_mag[nz]*W_mag(a, fK, nz);
}

double W_source(double a, int nz, double hoverh0)
{
  if(nz < -1 || nz > tomo.shear_Nbin - 1) 
  {
    log_fatal("invalid bin input ni = %d", nz);
    exit(1);
  }
  return zdistr_photoz(1./a - 1., (int) nz) * hoverh0;
}

double f_rsd(double a) 
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  const double z = 1.0/a-1.0;
  return f_growth(z);
}

double W_RSD(double l, double a0, double a1, int nz) 
{
  if(!(a0>0) || !(a0<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  if(!(a1>0) || !(a1<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  if(nz < -1 || nz > tomo.clustering_Nbin - 1) 
  {
    log_fatal("invalid bin input ni = %d", nz);
    exit(1);
  }
  double w = (1 + 8. * l) / ((2. * l + 1.) * (2. * l + 1.)) *
    pf_photoz(1. / a0 - 1., nz) * hoverh0(a0) * f_rsd(a0);
  
  w -= 4. / (2 * l + 3.) * sqrt((2 * l + 1.) / (2 * l + 3.)) *
    pf_photoz(1. / a1 - 1., nz) * hoverh0(a1) * f_rsd(a1);
  
  return w;
}

double W_HOD(double a, int nz, double hoverh0) 
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  if(nz < -1 || nz > tomo.clustering_Nbin - 1) 
  {
    log_fatal("invalid bin input ni = %d", nz);
    exit(1);
  }
  const double z = 1. / a - 1.;
  return pf_photoz(z, nz) * hoverh0;
}

/*
double W_cluster(int nz, double a, double hoverh0)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  const double z = 1. / a - 1.;
  return zdistr_cluster(nz, z)*hoverh0;
}
*/

double W_k(double a, double fK) 
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  return (1.5*cosmology.Omega_m*fK/a)*g_cmb(a);
}