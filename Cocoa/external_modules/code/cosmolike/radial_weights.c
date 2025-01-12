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
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  if (nz < 0 || nz > redshift.shear_nbin - 1) 
  {
    log_fatal("invalid bin input ni = %d", nz);
    exit(1);
  }
  return (1.5*cosmology.Omega_m*fK/a) * g_tomo(a, nz);
}

double W2_kappa(double a, double fK, int nz) 
{
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  if (nz < 0 || nz > redshift.shear_nbin - 1) 
  {
    log_fatal("invalid bin input ni = %d", nz);
    exit(1);
  }
  const double tmp = (1.5*cosmology.Omega_m*fK/a);
  return tmp*tmp*g2_tomo(a, nz);
}

double W_mag(double a, double fK, int nz) 
{
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  if (nz < 0 || nz > redshift.clustering_nbin - 1) 
  {
    log_fatal("invalid bin input ni = %d (max %d)", nz, redshift.clustering_nbin);
    exit(1);
  }
  return (1.5 * cosmology.Omega_m * fK / a) * g_lens(a, nz);
}

double W_gal(double a, int ni, double hoverh0) 
{
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  if (ni < 0 || ni > redshift.clustering_nbin - 1) 
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }
  const double z = 1. / a - 1;
  return pf_photoz(z, ni) * hoverh0;
}

double W_source(double a, int ni, double hoverh0)
{
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  if (ni < 0 || ni > redshift.shear_nbin - 1) 
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }

  const double z = 1.0 / a - 1.0;
  return zdistr_photoz(z, ni) * hoverh0;
}

double f_rsd(double a) 
{
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  const double z = 1.0 / a - 1.0;
  return f_growth(z);
}

double W_RSD(double l, double a0, double a1, int ni) 
{
  if (!(a0>0) || !(a0<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  if (!(a1>0) || !(a1<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  if (ni < -1 || ni > redshift.clustering_nbin - 1) 
  {
    log_fatal("invalid bin input ni = %d", ni);
    exit(1);
  }
  double w = (1 + 8. * l) / ((2. * l + 1.) * (2. * l + 1.)) *
    pf_photoz(1. / a0 - 1., ni) * hoverh0(a0) * f_rsd(a0);
  
  w -= 4. / (2 * l + 3.) * sqrt((2 * l + 1.) / (2 * l + 3.)) *
    pf_photoz(1. / a1 - 1., ni) * hoverh0(a1) * f_rsd(a1);
  
  return w;
}

double W_k(double a, double fK) 
{
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  return (1.5*cosmology.Omega_m*fK/a)*g_cmb(a);
}

double W_y(double a) // efficiency weight function for Compton-y
{ // sigma_Th /(m_e*c^2) / a^2 , see Eq.D9 of 2005.00009.
  const double real_coverH0 = cosmology.coverH0 / cosmology.h0; // unit Mpc
  const double sigma_Th = 7.012e-74 / (real_coverH0*real_coverH0); // from Mpc^2 to (c/H0)^2
  const double E_e = 0.511*cosmology.h0*5.6131e-38;  // from MeV to [G(M_solar/h)^2/(c/H0)]
  return sigma_Th/(E_e*a*a); //  dim = [comoving L]^2 / [Energy], 
                             // units = [c/H0]^2 / [G(M_solar/h)^2/(c/H0)]
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

/*
double W_mag_cluster(double a, double fK, int nz, int nl)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  if(nz < -1 || nz > redshift.cluster_nbin - 1) 
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
  return wmag;
}
*/