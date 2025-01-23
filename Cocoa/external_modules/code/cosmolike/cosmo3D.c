#include <assert.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "log.c/src/log.h"

#include "basics.h"
#include "baryons.h"
#include "cosmo3D.h"
#include "structs.h"

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Background
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double chi(double io_a)
{
  struct chis r = chi_all(io_a);
  return r.chi;
}

double dchi_da(double io_a)
{
  struct chis r = chi_all(io_a);
  return r.dchida;
}

struct chis chi_all(double io_a)
{
  double out_dchi[4];
  const double redshift = 1.0/io_a - 1.0;

  int j = 0;
  {
    size_t ilo = 0;
    size_t ihi = cosmology.chi_nz - 1;
    while (ihi > ilo + 1)
    {
      size_t ll = (ihi + ilo)/2;
      if (cosmology.chi[0][ll] > redshift)
        ihi = ll;
      else
        ilo = ll;
    }
    j = ilo;
  }

  double dy = (redshift              - cosmology.chi[0][j])/
              (cosmology.chi[0][j+1] - cosmology.chi[0][j]);
  out_dchi[0] = cosmology.chi[1][j] + dy*(cosmology.chi[1][j+1] - cosmology.chi[1][j]);

  if (j>0)
  {
    const double up = (cosmology.chi[1][j+2] - cosmology.chi[1][j])/
                      (cosmology.chi[0][j+2] - cosmology.chi[0][j]);
    
    const double down = (cosmology.chi[1][j+1] - cosmology.chi[1][j-1])/
                        (cosmology.chi[0][j+1] - cosmology.chi[0][j-1]);
    out_dchi[1] = down + dy*(up-down);
  }
  else 
  {
    const double up = (cosmology.chi[1][j+2] - cosmology.chi[1][j])/
                      (cosmology.chi[0][j+2] - cosmology.chi[0][j]);
    
    const double down = (cosmology.chi[1][j+1] - cosmology.chi[1][j])/
                        (cosmology.chi[0][j+1] - cosmology.chi[0][j]);
    out_dchi[1] = down + dy*(up-down);
  }

  // convert from (Mpc/h) to (Mpc/h)/(c/H0=100)^3 (dimensioneless)
  out_dchi[1] = (out_dchi[1]/cosmology.coverH0);
  // convert from d\chi/dz to d\chi/da
  out_dchi[1] = out_dchi[1]/(io_a*io_a);

  struct chis result;
  result.chi = out_dchi[0]/cosmology.coverH0;
  result.dchida = out_dchi[1];
  return result;
}

double dchi_dz(double io_a)
{
  return (io_a*io_a)*dchi_da(io_a);
}

double hoverh0(double io_a)
{
  return 1.0/dchi_dz(io_a);
}

double hoverh0v2(double io_a, double dchida)
{
  return 1.0/((io_a*io_a)*dchida);
}

double a_chi(double io_chi)
{
  // convert from (Mpc/h)/(c/H0=100)^3 (dimensioneless) to (Mpc/h)
  io_chi *= cosmology.coverH0;

  int j = 0;
  {
    size_t ilo = 0;
    size_t ihi = cosmology.chi_nz-1;
    while (ihi > ilo + 1)
    {
      size_t ll = (ihi + ilo)/2;
      if (cosmology.chi[1][ll] > io_chi)
        ihi = ll;
      else
        ilo = ll;
    }
    j = ilo;
  }

  double dy = (io_chi                - cosmology.chi[1][j])/
              (cosmology.chi[1][j+1] - cosmology.chi[1][j]);
  double z  = cosmology.chi[0][j] + dy*(cosmology.chi[0][j+1] - cosmology.chi[0][j]);

  return 1.0/(1.0+z);
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Growth Factor
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double growfac(double io_a)
{
  return norm_growfac(io_a, true);
}

double norm_growfac(const double a, const bool normalize_z0)
{
  double growfact1;
  {
    const double z = 0.0;
    int j = 0;
    {
      size_t ilo = 0;
      size_t ihi = cosmology.G_nz-1;
      while (ihi>ilo+1) {
        size_t ll = (ihi+ilo)/2;
        if(cosmology.G[0][ll]>z)
          ihi = ll;
        else
          ilo = ll;
      }
      j = ilo;
    }
    const double dy = (z                   - cosmology.G[0][j])/
                      (cosmology.G[0][j+1] - cosmology.G[0][j]);
    
    growfact1 = cosmology.G[1][j] + dy*(cosmology.G[1][j+1] - cosmology.G[1][j]);
  }

  const double z = 1.0/a-1.0;

  int j = 0;
  {
    size_t ilo = 0;
    size_t ihi = cosmology.G_nz-1;
    while (ihi>ilo+1)
    {
      size_t ll = (ihi+ilo)/2;
      if(cosmology.G[0][ll]>z)
        ihi = ll;
      else
        ilo = ll;
    }
    j = ilo;
  }

  const double dy = (z                   - cosmology.G[0][j])/
                    (cosmology.G[0][j+1] - cosmology.G[0][j]);

  const double G = cosmology.G[1][j] + dy*(cosmology.G[1][j+1] - cosmology.G[1][j]);

  if(normalize_z0)
    return (G*a)/growfact1; // Growth D = G * a
  else
    return G*a; // Growth D = G * a
}

double f_growth(const double z)
{
  int j = 0;
  {
    size_t ilo = 0;
    size_t ihi = cosmology.G_nz-1;
    while (ihi>ilo+1)
    {
      size_t ll = (ihi+ilo)/2;
      if(cosmology.G[0][ll]>z)
        ihi = ll;
      else
        ilo = ll;
    }
    j = ilo;
  }

  const double dy = (z                   - cosmology.G[0][j])/
                    (cosmology.G[0][j+1] - cosmology.G[0][j]);
                    
  const double G = cosmology.G[1][j] + dy*(cosmology.G[1][j+1] - cosmology.G[1][j]);

  const double dlnGdlnz = ((cosmology.G[1][j+1] - cosmology.G[1][j])/
                           (cosmology.G[0][j+1] - cosmology.G[0][j]))*z/G;
  
  const double dlnGdlna = -dlnGdlnz*(1+z)/z;

  return 1 + dlnGdlna; // Growth D = G * a
}

struct growths norm_growfac_all(const double a, const bool normalize_z0)
{
  double growfact1;
  {
    const double z = 0.0;
    int j = 0;
    {
      size_t ilo = 0;
      size_t ihi = cosmology.G_nz-1;
      while (ihi>ilo+1) 
      {
        size_t ll = (ihi+ilo)/2;
        if(cosmology.G[0][ll]>z)
          ihi = ll;
        else
          ilo = ll;
      }
      j = ilo;
    }
    const double dy = (z                   - cosmology.G[0][j])/
                      (cosmology.G[0][j+1] - cosmology.G[0][j]);
    
    growfact1 = cosmology.G[1][j] + dy*(cosmology.G[1][j+1] - cosmology.G[1][j]);
  }

  const double z = 1.0/a-1.0;

  int j = 0;
  {
    size_t ilo = 0;
    size_t ihi = cosmology.G_nz - 1;
    while (ihi>ilo+1)
    {
      size_t ll = (ihi+ilo)/2;
      if(cosmology.G[0][ll]>z)
        ihi = ll;
      else
        ilo = ll;
    }
    j = ilo;
  }

  const double dy = (z                   - cosmology.G[0][j])/
                    (cosmology.G[0][j+1] - cosmology.G[0][j]);
                    
  const double G = cosmology.G[1][j] + dy*(cosmology.G[1][j+1] - cosmology.G[1][j]);

  const double dlnGdlnz = ((cosmology.G[1][j+1] - cosmology.G[1][j])/
                           (cosmology.G[0][j+1] - cosmology.G[0][j]))*z/G;
  
  const double dlnGdlna = -dlnGdlnz*(1+z)/z;

  struct growths Gf;
  Gf.f = 1 + dlnGdlna; // Growth D = G * a

  if(normalize_z0)
    Gf.D = (G*a)/growfact1; 
  else
    Gf.D = (G*a);

  return Gf;
}

struct growths growfac_all(double io_a)
{
  return norm_growfac_all(io_a, true);
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Power Spectrum
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double p_lin(const double k, const double a)
{
  // convert from (x/Mpc/h - dimensioneless) to h/Mpc with x = c/H0 (Mpc)
  const double log10k = log10(k/cosmology.coverH0);
  const double z = 1.0/a-1.0;

  // logk = cosmology.lnPL[0:nk,cosmology.lnPL_nz]
  // z    = cosmology.lnPL[cosmology.lnPL_nk,0:nz]
  
  int i = 0;
  {
    size_t ilo = 0;
    size_t ihi = cosmology.lnPL_nk-1;
    while (ihi>ilo+1) 
    {
      size_t ll = (ihi+ilo)/2;
      if(cosmology.lnPL[ll][cosmology.lnPL_nz] > log10k)
        ihi = ll;
      else
        ilo = ll;
    }
    i = ilo;
  }

  int j = 0;
  {
    size_t ilo = 0;
    size_t ihi = cosmology.lnPL_nz-1;
    while (ihi>ilo+1) 
    {
      size_t ll = (ihi+ilo)/2;
      if(cosmology.lnPL[cosmology.lnPL_nk][ll] > z)
        ihi = ll;
      else
        ilo = ll;
    }
    j = ilo;
  }

  double dx = (log10k                                 - cosmology.lnPL[i][cosmology.lnPL_nz])/
              (cosmology.lnPL[i+1][cosmology.lnPL_nz] - cosmology.lnPL[i][cosmology.lnPL_nz]);

  double dy = (z                                     - cosmology.lnPL[cosmology.lnPL_nk][j])/
              (cosmology.lnPL[cosmology.lnPL_nk][j+1]- cosmology.lnPL[cosmology.lnPL_nk][j]);

  const double out_lnP =    (1-dx)*(1-dy)*cosmology.lnPL[i][j]
                          + (1-dx)*dy*cosmology.lnPL[i][j+1]
                          + dx*(1-dy)*cosmology.lnPL[i+1][j]
                          + dx*dy*cosmology.lnPL[i+1][j+1];

  // convert from (Mpc/h)^3 to (Mpc/h)^3/(c/H0=100)^3 (dimensioneless)
  return exp(out_lnP)/(cosmology.coverH0*cosmology.coverH0*cosmology.coverH0);
}

double p_nonlin(const double k, const double a)
{
  // convert from (x/Mpc/h - dimensioneless) to h/Mpc with x = c/H0 (Mpc)
  const double log10k = log10(k/cosmology.coverH0);
  const double z = 1.0/a-1.0;

  // logk = cosmology.lnP[0:nk,cosmology.lnP_nz]
  // z    = cosmology.lnP[cosmology.lnP_nk,0:nz]
  
  int i = 0;
  {
    size_t ilo = 0;
    size_t ihi = cosmology.lnP_nk-1;
    while (ihi>ilo+1) 
    {
      size_t ll = (ihi+ilo)/2;
      if(cosmology.lnP[ll][cosmology.lnP_nz] > log10k)
        ihi = ll;
      else
        ilo = ll;
    }
    i = ilo;
  }

  int j = 0;
  {
    size_t ilo = 0;
    size_t ihi = cosmology.lnP_nz-1;
    while (ihi>ilo+1) 
    {
      size_t ll = (ihi+ilo)/2;
      if(cosmology.lnP[cosmology.lnP_nk][ll] > z)
        ihi = ll;
      else
        ilo = ll;
    }
    j = ilo;
  }

  double dx = (log10k                                 - cosmology.lnP[i][cosmology.lnP_nz])/
              (cosmology.lnP[i+1][cosmology.lnP_nz] - cosmology.lnP[i][cosmology.lnP_nz]);

  double dy = (z                                     - cosmology.lnP[cosmology.lnP_nk][j])/
              (cosmology.lnP[cosmology.lnP_nk][j+1]- cosmology.lnP[cosmology.lnP_nk][j]);

  const double out_lnP =    (1-dx)*(1-dy)*cosmology.lnP[i][j]
                          + (1-dx)*dy*cosmology.lnP[i][j+1]
                          + dx*(1-dy)*cosmology.lnP[i+1][j]
                          + dx*dy*cosmology.lnP[i+1][j+1];

  // convert from (Mpc/h)^3 to (Mpc/h)^3/(c/H0=100)^3 (dimensioneless)
  return exp(out_lnP)/(cosmology.coverH0*cosmology.coverH0*cosmology.coverH0);
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

double Pdelta(double io_kNL, double io_a) 
{
  double out_PK;
  static int P_type = -1;
  if (P_type == -1) 
  {
    if (strcmp(pdeltaparams.runmode,"linear") == 0) 
    {
      P_type = 3;
    }
  }
  switch (P_type) 
  {
    case 3:
      out_PK = p_lin(io_kNL, io_a);
      break;
    default:
      out_PK = p_nonlin(io_kNL, io_a);
      break;
  }
  return out_PK;
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

// calculating the angular diameter distance f_K
// BS01 2.4, 2.30: f_K is a radial function that, depending on the curvature of
// the Universe, is a trigonometric, linear, or hyperbolic function of chi
double f_K(double chi) 
{
  double K, K_h, f;
  K = (cosmology.Omega_m + cosmology.Omega_v - 1.);
  if (K > 1e-6) 
  { // open
    K_h = sqrt(K); // K in units H0/c see BS eq. 2.30
    f = 1. / K_h * sin(K_h * chi);
  } else if (K < -1e-6) 
  { // closed
    K_h = sqrt(-K);
    f = 1. / K_h * sinh(K_h * chi);
  } else 
  { // flat
    f = chi;
  }
  return f;
}

// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// Baryons
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------

// return P(k)_bary/P(k)_DMO from hydro sims ;
double PkRatio_baryons(double k_NL, double a)
{
  if (bary.is_Pk_bary == 0)
  {
    return 1.;
  } else
  {
    const double kintern = k_NL/cosmology.coverH0;
    double result;
    int status = gsl_interp2d_eval_extrap_e(bary.interp2d, bary.logk_bins,
      bary.a_bins, bary.log_PkR, log10(kintern), a, NULL, NULL, &result);
    if (status)
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
    return pow(10.0, result);
  }
}

// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// MODIFIED GRAVITY
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------

double MG_Sigma(double a __attribute__((unused))) {
  return 0.0;
}


double int_for_sigma2(double x, void* params) // inner integral
{
  double* ar = (double*) params;
  const double R = ar[0];
  const double a = ar[1];
  const double PK = p_lin(x/R, a);
  
  gsl_sf_result J1;
  int status = gsl_sf_bessel_j1_e(x, &J1);
  if (status) 
  {
    log_fatal(gsl_strerror(status));
    exit(1);
  }
  const double tmp = 3.0*J1.val/ar[0];
  return PK*tmp*tmp/(ar[0] * 2.0 * M_PI * M_PI);
}

double sigma2_nointerp(
    const double M, 
    const double a, 
    const int init_static_vars_only
  ) 
{
  static double cache[MAX_SIZE_ARRAYS];
  static gsl_integration_glfixed_table* w = NULL;

  if (w == NULL || fdiff(cache[0], Ntable.random))
  {
    const size_t szint = 500 + 500 * (Ntable.high_def_integration);
    if (w != NULL)  gsl_integration_glfixed_table_free(w);
    w = malloc_gslint_glfixed(szint);
    cache[0] = Ntable.random;
  }
  
  double ar[1] = {pow(0.75*M/(M_PI*cosmology.rho_crit*cosmology.Omega_m),1./3.)};
  const double xmin = 0;
  const double xmax = 14.1;

  double res;
  if (init_static_vars_only == 1)
    res = int_for_sigma2((xmin+xmax)/2.0, (void*) ar);
  else
  {
    gsl_function F;
    F.params = (void*) ar;
    F.function = int_for_sigma2;
    res = gsl_integration_glfixed(&F, xmin, xmax, w);
  }
  return res;
}

double sigma2(const double M) 
{
  static double cache[MAX_SIZE_ARRAYS];
  static double* table;
  static double lim[3];

  if (table == NULL || fdiff(cache[1], Ntable.random))
  {
    if (table != NULL) free(table);
    table = (double*) malloc(sizeof(double)*Ntable.N_M);
    lim[0] = log(limits.M_min);
    lim[1] = log(limits.M_max);
    lim[2] = (lim[1] - lim[0])/((double) Ntable.N_M - 1.0);
  } 
  if (fdiff(cache[0], cosmology.random) || fdiff(cache[1], Ntable.random))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    { // sigma2_nointerp(M, a, init_static_vars_only) 
      double init = sigma2_nointerp(exp(lim[0]), 1.0, 1);
    }
    #pragma GCC diagnostic pop
    
    #pragma omp parallel for
    for (int i=0; i<Ntable.N_M; i++) 
      table[i] = log(sigma2_nointerp(exp(lim[0] + i*lim[2]), 1.0, 0));
    
    cache[0] = cosmology.random;
    cache[1] = Ntable.random;
  }
  return exp(interpol1d(table, Ntable.N_M, lim[0], lim[1], lim[2], log(M)));
}


