#include <assert.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_integration.h>
#include "../cfftlog/cfftlog.h"
#include <fftw3.h>

#include "bias.h"
#include "basics.h"
#include "cfastpt/cfastpt.h"
#include "cosmo3D.h"
#include "cosmo2D.h"
#include "halo.h"
#include "IA.h"
#include "pt_cfastpt.h"
#include "radial_weights.h"
#include "redshift_spline.h"
#include "structs.h"

#include "log.c/src/log.h"
 
static int include_HOD_GX = 0; // 0 or 1

static int include_RSD_GS = 0; // 0 or 1 
static int include_RSD_GG = 1; // 0 or 1 
static int include_RSD_GK = 0; // 0 or 1
static int include_RSD_GY = 0; // 0 or 1

double beam_cmb(const double l)
{
  const double sigma = cmb.fwhm/sqrt(16.0*log(2.0));
  return exp(-l*(l+1.0)*sigma*sigma);
}

double w_pixel(const double ell)
{
  static int lbins;
  static double* cl_pixel =0;

  const int l = (int) ceil(ell);

  if (cl_pixel == 0)
  {
    FILE* ein = fopen(cmb.pathHealpixWinFunc, "r");
    if (ein == NULL)
    {
      log_info("File not open (%s) \n", cmb.pathHealpixWinFunc);
      lbins = 0;
      return 1.0;
    }

    lbins = line_count(cmb.pathHealpixWinFunc);
      
    cl_pixel = malloc(sizeof(double)*lbins);
            
    for (int i=0; i<lbins; i++)
    {
      int tmp;
      const int status = fscanf(ein, "%d %le\n", &tmp, &cl_pixel[i]);
      if (status)
      {
        log_debug("Read file failed (%s) \n", cmb.pathHealpixWinFunc);
        exit(1);
      }
    }

    fclose(ein);
  }

  return (lbins > 0) ? ((l < lbins) ? cl_pixel[l] : 0.0) : 1.0;
}

static int has_b2_galaxies()
{
  int res = 0;
  for (int i=0; i<redshift.clustering_nbin; i++) 
    if (nuisance.gb[1][i])
      res = 1;
  return res;
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Correlation Functions (real Space) - Full Sky - bin average
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double xi_pm_tomo(
    const int pm, 
    const int nt, 
    const int ni, 
    const int nj, 
    const int limber
  )
{  
  static double*** Glpm = NULL; //Glpm[0] = Gl+, Glpm[1] = Gl-
  static double** xipm = NULL;  //xipm[0] = xi+, xipm[1] = xi-
  static double cache[MAX_SIZE_ARRAYS];
  
  if (Ntable.Ntheta == 0)
  {
    log_fatal("Ntable.Ntheta not initialized");
    exit(1);
  }

  const int NSIZE = tomo.shear_Npowerspectra;

  if (Glpm == NULL || xipm == NULL || fdiff(cache[4], Ntable.random))
  {
    if (Glpm != NULL) free(Glpm);
    Glpm = (double***) malloc3d(2, Ntable.Ntheta, limits.LMAX);
    if (xipm != NULL) free(xipm);
    xipm = (double**) malloc2d(2, NSIZE*Ntable.Ntheta);
    
    double*** P = (double***) malloc3d(4, Ntable.Ntheta, limits.LMAX + 1);
    double** Pmin  = P[0]; double** Pmax  = P[1];
    double** dPmin = P[2]; double** dPmax = P[3];

    double xmin[Ntable.Ntheta];
    double xmax[Ntable.Ntheta];
    for (int i=0; i<Ntable.Ntheta; i++)
    { // Cocoa: dont thread (init of static variables inside set_bin_average)
      bin_avg r = set_bin_average(i, 0);
      xmin[i] = r.xmin;
      xmax[i] = r.xmax;
    }

    #pragma omp parallel for collapse(2)
    for (int i=0; i<Ntable.Ntheta; i++)
    {
      for (int l=0; l<(limits.LMAX+1); l++)
      {
        bin_avg r   = set_bin_average(i, l);
        Pmin[i][l]  = r.Pmin;
        Pmax[i][l]  = r.Pmax;
        dPmin[i][l] = r.dPmin;
        dPmax[i][l] = r.dPmax;
      }
    }

    const int lmin = 1;
    for (int i=0; i<Ntable.Ntheta; i++)
    {
      for (int l=0; l<lmin; l++)
      {
        Glpm[0][i][l] = 0.0;
        Glpm[1][i][l] = 0.0;
      }
    }

    #pragma omp parallel for collapse(2)
    for (int i=0; i<Ntable.Ntheta; i++)
    {
      for (int l=lmin; l<limits.LMAX; l++)
      {
        Glpm[0][i][l] = (2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(
          -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[i][l-1]-Pmax[i][l-1])
          -l*(l-1.)*(2.-l)/2 * (xmin[i]*Pmin[i][l]-xmax[i]*Pmax[i][l])
          +l*(l-1.)/(2.*l+1) * (Pmin[i][l+1]-Pmax[i][l+1])
          +(4-l)*(dPmin[i][l]-dPmax[i][l])
          +(l+2)*(xmin[i]*dPmin[i][l-1] - xmax[i]*dPmax[i][l-1] - Pmin[i][l-1] + Pmax[i][l-1])
          +2*(l-1)*(xmin[i]*dPmin[i][l] - xmax[i]*dPmax[i][l] - Pmin[i][l] + Pmax[i][l])
          -2*(l+2)*(dPmin[i][l-1]-dPmax[i][l-1])
        )/(xmin[i]-xmax[i]);

        Glpm[1][i][l] = (2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(
          -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[i][l-1]-Pmax[i][l-1])
          -l*(l-1.)*(2.-l)/2 * (xmin[i]*Pmin[i][l]-xmax[i]*Pmax[i][l])
          +l*(l-1.)/(2.*l+1)* (Pmin[i][l+1]-Pmax[i][l+1])
          +(4-l)*(dPmin[i][l]-dPmax[i][l])
          +(l+2)*(xmin[i]*dPmin[i][l-1] - xmax[i]*dPmax[i][l-1] - Pmin[i][l-1] + Pmax[i][l-1])
          -2*(l-1)*(xmin[i]*dPmin[i][l] - xmax[i]*dPmax[i][l] - Pmin[i][l] + Pmax[i][l])
          +2*(l+2)*(dPmin[i][l-1]-dPmax[i][l-1])
          )/(xmin[i]-xmax[i]);
      }
    }

    free(P);
  }

  if (fdiff(cache[0], cosmology.random) ||
      fdiff(cache[1], nuisance.random_photoz_shear) ||
      fdiff(cache[2], nuisance.random_ia) ||
      fdiff(cache[3], redshift.random_shear) ||
      fdiff(cache[4], Ntable.random))
  {
    if (limber == 1)
    {
      // Cl_EE = Cl[0], Cl_BB = Cl[1]
      double*** Cl = (double***) malloc3d(2, NSIZE, limits.LMAX);
      
      const int lmin = 1;
      for (int i=0; i<NSIZE; i++)
      {
        for (int l=0; l<lmin; l++)
        {
          Cl[0][i][l] = 0.0;
          Cl[1][i][l] = 0.0;
        }
      }

      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wunused-variable"
      { // init static variables inside the C_XY_limber function
        double ini = C_ss_tomo_limber(limits.LMIN_tab+1,Z1(0),Z2(0),1);
      }
      #pragma GCC diagnostic pop

      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<NSIZE; nz++) 
      {
        for (int l=lmin; l<limits.LMIN_tab; l++)
        {
          const int Z1NZ = Z1(nz);
          const int Z2NZ = Z2(nz);
          Cl[0][nz][l] = C_ss_tomo_limber_nointerp(l, Z1NZ, Z2NZ, 1, 0);
          Cl[1][nz][l] = C_ss_tomo_limber_nointerp(l, Z1NZ, Z2NZ, 0, 0);
        }
      }
      
      #pragma omp parallel for collapse(3)
      for (int i=0; i<2; i++) 
        for (int nz=0; nz<NSIZE; nz++) 
          for (int l=limits.LMIN_tab; l<limits.LMAX; l++)
            Cl[i][nz][l] = C_ss_tomo_limber(l, Z1(nz), Z2(nz), 1-i);
      
      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<NSIZE; nz++)
      {
        for (int i=0; i<Ntable.Ntheta; i++)
        {
          const int q = nz*Ntable.Ntheta + i;
          xipm[0][q] = 0;
          xipm[1][q] = 0;  
          for (int l=lmin; l<limits.LMAX; l++)
          {
            xipm[0][q] += Glpm[0][i][l] * (Cl[0][nz][l] + Cl[1][nz][l]);
            xipm[1][q] += Glpm[1][i][l] * (Cl[0][nz][l] - Cl[1][nz][l]);
          }
        }
      }
      free(Cl);
    }
    else
    {
      log_fatal("NonLimber not implemented");
      exit(1);
    }

    cache[0] = cosmology.random;
    cache[1] = nuisance.random_photoz_shear;
    cache[2] = nuisance.random_ia;
    cache[3] = redshift.random_shear;
    cache[4] = Ntable.random;
  }

  if (nt < 0 || nt > Ntable.Ntheta - 1)
  {
    log_fatal("error in selecting bin number nt = %d (max %d)", nt, Ntable.Ntheta);
    exit(1); 
  }
  
  if (ni < 0 || 
      ni > redshift.shear_nbin - 1 || 
      nj < 0 || 
      nj > redshift.shear_nbin - 1)
  {
    log_fatal("error in selecting bin number (ni, nj) = [%d,%d]", ni, nj);
    exit(1);
  }

  const int ntomo = N_shear(ni, nj);
  const int q = ntomo*Ntable.Ntheta + nt;
  if (q < 0 || q > NSIZE*Ntable.Ntheta - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  }
  return (pm > 0) ? xipm[0][q] : xipm[1][q];
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double w_gammat_tomo(const int nt, const int ni, const int nj, const int limber)
{
  static double** Pl = NULL;
  static double* w_vec = NULL;
  static double cache[MAX_SIZE_ARRAYS];

  if (Ntable.Ntheta == 0)
  {
    log_fatal("Ntable.Ntheta not initialized");
    exit(1);
  }

  const int NSIZE = tomo.ggl_Npowerspectra;

  if (Pl == NULL || w_vec == NULL || fdiff(cache[6], Ntable.random))
  {
    if (Pl != NULL) free(Pl);
    Pl = (double**) malloc2d(Ntable.Ntheta, limits.LMAX);;
    
    if (w_vec != NULL) free(w_vec);
    w_vec = (double*) calloc1d(NSIZE*Ntable.Ntheta);

    double*** P = (double***) malloc3d(2, Ntable.Ntheta, limits.LMAX + 1);
    double** Pmin  = P[0]; double** Pmax  = P[1];

    double xmin[Ntable.Ntheta];
    double xmax[Ntable.Ntheta];
    for (int i=0; i<Ntable.Ntheta; i++)
    { // Cocoa: dont thread (init of static variables inside set_bin_average)
      bin_avg r = set_bin_average(i,0);
      xmin[i] = r.xmin;
      xmax[i] = r.xmax;
    }

    #pragma omp parallel for collapse(2)
    for (int i=0; i<Ntable.Ntheta; i ++)
    {
      for (int l=0; l<(limits.LMAX+1); l++)
      {
        bin_avg r = set_bin_average(i, l);
        Pmin[i][l] = r.Pmin;
        Pmax[i][l] = r.Pmax;
      }
    }

    const int lmin = 1;
    for (int i=0; i<Ntable.Ntheta; i++)
      for (int l=0; l<lmin; l++)
        Pl[i][l] = 0.0;
    
    #pragma omp parallel for collapse(2)
    for (int i=0; i<Ntable.Ntheta; i++)
    {
      for (int l=lmin; l<limits.LMAX; l++)
      {
        Pl[i][l] = (2.*l+1)/(4.*M_PI*l*(l+1)*(xmin[i]-xmax[i]))
          *((l+2./(2*l+1.))*(Pmin[i][l-1]-Pmax[i][l-1])
          +(2-l)*(xmin[i]*Pmin[i][l]-xmax[i]*Pmax[i][l])
          -2./(2*l+1.)*(Pmin[i][l+1]-Pmax[i][l+1]));
      }
    }

    free(P);
  }

  if (fdiff(cache[0], cosmology.random) ||
      fdiff(cache[1], nuisance.random_photoz_shear) ||
      fdiff(cache[2], nuisance.random_photoz_clustering) ||
      fdiff(cache[3], nuisance.random_ia) ||
      fdiff(cache[4], redshift.random_shear) ||
      fdiff(cache[5], redshift.random_clustering) ||
      fdiff(cache[6], Ntable.random) ||
      fdiff(cache[7], nuisance.random_galaxy_bias))
  {    
    double** Cl = (double**) malloc2d(NSIZE, limits.LMAX);

    const int lmin = 1;
    for (int i=0; i<NSIZE; i++)
      for (int l=0; l<lmin; l++)
        Cl[i][l] = 0.0;
      
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    { // init static variables inside the C_XY_limber function
      double ini = C_gs_tomo_limber(limits.LMIN_tab + 1, ZL(0), ZS(0));
    }
    #pragma GCC diagnostic pop

    if (limber == 1)
    {
      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<NSIZE; nz++)
        for (int l=lmin; l<limits.LMIN_tab; l++)
          Cl[nz][l] = C_gs_tomo_limber_nointerp(l, ZL(nz), ZS(nz), 0);
      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<NSIZE; nz++)
        for (int l=limits.LMIN_tab; l<limits.LMAX; l++)      
          Cl[nz][l] = C_gs_tomo_limber(l, ZL(nz), ZS(nz));
    }
    else
    {
      log_fatal("NonLimber not implemented");
      exit(1);
    }
      
    #pragma omp parallel for collapse(2)
    for (int nz=0; nz<NSIZE; nz++)
    {
      for (int i=0; i<Ntable.Ntheta; i++)
      {
        const int q = nz*Ntable.Ntheta+i;
        w_vec[q] = 0;
        for (int l=lmin; l<limits.LMAX; l++)
          w_vec[q] += Pl[i][l]*Cl[nz][l];
      }
    }

    free(Cl);

    cache[0] = cosmology.random;
    cache[1] = nuisance.random_photoz_shear;
    cache[2] = nuisance.random_photoz_clustering;
    cache[3] = nuisance.random_ia;
    cache[4] = redshift.random_shear;
    cache[5] = redshift.random_clustering;
    cache[6] = Ntable.random;
    cache[7] = nuisance.random_galaxy_bias;
  }

  if (nt < 0 || nt > Ntable.Ntheta - 1)
  {
    log_fatal("error in selecting bin number nt = %d (max %d)", nt, Ntable.Ntheta);
    exit(1); 
  }
  if (ni < -1 || 
      ni > redshift.clustering_nbin - 1 || 
      nj < -1 || 
      nj > redshift.shear_nbin - 1)
  {
    log_fatal("error in selecting bin number (ni, nj) = [%d,%d]", ni, nj);
    exit(1);
  }
  
  const int q = N_ggl(ni, nj)*Ntable.Ntheta + nt;
  if (q < 0 || q > NSIZE*Ntable.Ntheta - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  }

  return w_vec[q];
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double w_gg_tomo(const int nt, const int ni, const int nj, const int limber)
{
  static double** Pl = NULL;
  static double* w_vec = NULL;
  static double cache[MAX_SIZE_ARRAYS];

  if (Ntable.Ntheta == 0)
  {
    log_fatal("Ntable.Ntheta not initialized");
    exit(1);
  }

  const int NSIZE = tomo.clustering_Npowerspectra;

  if (Pl == NULL || 
      w_vec == NULL || 
      fdiff(cache[3], Ntable.random))
  {
    if (Pl != NULL) free(Pl);
    Pl = (double**) malloc2d(Ntable.Ntheta, limits.LMAX);;

    if (w_vec != NULL) free(w_vec);
    w_vec = (double*) calloc1d(NSIZE*Ntable.Ntheta);

    double*** P = (double***) malloc3d(2, Ntable.Ntheta, limits.LMAX + 1);
    double** Pmin  = P[0]; double** Pmax  = P[1];

    double xmin[Ntable.Ntheta];
    double xmax[Ntable.Ntheta];
    for (int i=0; i<Ntable.Ntheta; i ++)
    { // Cocoa: dont thread (init of static variables inside set_bin_average)
      bin_avg r = set_bin_average(i,0);
      xmin[i] = r.xmin;
      xmax[i] = r.xmax;
    }

    #pragma omp parallel for collapse(2)
    for (int i=0; i<Ntable.Ntheta; i++)
    {
      for (int l=0; l<(limits.LMAX+1); l++)
      {
        bin_avg r = set_bin_average(i,l);
        Pmin[i][l] = r.Pmin;
        Pmax[i][l] = r.Pmax;
      }
    }

    const int lmin = 1;
    for (int i=0; i<Ntable.Ntheta; i++)
      for (int l=0; l<lmin; l++)
        Pl[i][l] = 0.0;

    #pragma omp parallel for collapse(2)
    for (int i=0; i<Ntable.Ntheta; i++)
    {
      for (int l=lmin; l<limits.LMAX; l++)
      { 
        const double tmp = (1.0/(xmin[i] - xmax[i]))*(1. / (4.0 * M_PI));
        Pl[i][l] = tmp*(Pmin[i][l + 1] - Pmax[i][l + 1] 
                        - Pmin[i][l - 1] + Pmax[i][l - 1]);
      }
    }

    free(P);
  }

  if (fdiff(cache[0], cosmology.random) || 
      fdiff(cache[1], nuisance.random_photoz_clustering) ||
      fdiff(cache[2], redshift.random_clustering) ||
      fdiff(cache[3], Ntable.random) ||
      fdiff(cache[4], nuisance.random_galaxy_bias))
  {        
    double** Cl = (double**) malloc2d(NSIZE, limits.LMAX);
    
    const int lmin = 1;
    for (int i=0; i<NSIZE; i++)
      for (int l=0; l<lmin; l++)
        Cl[i][l] = 0.0;
    
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    { // init static variables inside the C_XY_limber function
      double ini = C_gg_tomo_limber(limits.LMIN_tab + 1, 0, 0);
    }
    #pragma GCC diagnostic pop
    
    if (limber == 1)
    {
      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<NSIZE; nz++)
        for (int l=lmin; l<limits.LMIN_tab; l++)
          Cl[nz][l] = C_gg_tomo_limber_nointerp(l, nz, nz, 0);
      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<NSIZE; nz++)
        for (int l=limits.LMIN_tab; l<limits.LMAX; l++)
          Cl[nz][l] = C_gg_tomo_limber( l, nz, nz);
    }
    else
    {
      for (int nz=0; nz<NSIZE; nz++) // NONLIMBER PART
      { 
        const int L = 1;
        const double tolerance = 0.01;      // required fractional accuracy in C(l)
        const double dev = 10. * tolerance; // will be diff  exact vs Limber init to
                                            // large value in order to start while loop
        const int Z1 = nz; // cross redshift bin not supported so not using ZCL1(k)
        const int Z2 = nz; // cross redshift bin not supported so not using ZCL2(k)
        
        C_cl_tomo(L, Z1, Z2, Cl[nz], dev, tolerance);
      }
      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<NSIZE; nz++) // LIMBER PART
        for (int l=limits.LMAX_NOLIMBER+1; l<limits.LMAX; l++)
          Cl[nz][l] = C_gg_tomo_limber( l, nz, nz);
    }

    #pragma omp parallel for collapse(2)
    for (int nz=0; nz<NSIZE; nz++)
    {
      for (int i=0; i<Ntable.Ntheta; i++)
      {
        const int q = nz*Ntable.Ntheta + i;
        w_vec[q] = 0;
        for (int l=lmin; l<limits.LMAX; l++)
          w_vec[q] += Pl[i][l]*Cl[nz][l];
      }
    }

    free(Cl);

    cache[0] = cosmology.random;
    cache[1] = nuisance.random_photoz_clustering;
    cache[2] = redshift.random_clustering;
    cache[3] = Ntable.random;
    cache[4] = nuisance.random_galaxy_bias;
  }

  if (nt < 0 || nt > Ntable.Ntheta - 1)
  {
    log_fatal("error in selecting bin number nt = %d (max %d)", nt, Ntable.Ntheta);
    exit(1); 
  }
  if (ni < -1 || 
      ni > redshift.clustering_nbin - 1 || 
      nj < -1 || 
      nj > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number (ni, nj) = [%d,%d]", ni, nj);
    exit(1);
  }
  if (ni != nj)
  {
    log_fatal("ni != nj tomography not supported");
    exit(1);
  }
  
  const int q = ni * Ntable.Ntheta + nt;
  if (q  < 0 || q > NSIZE*Ntable.Ntheta - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  }  
  
  return w_vec[q];
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double w_gk_tomo(const int nt, const int ni, const int limber)
{
  static double** Pl = NULL;
  static double* w_vec = NULL;
  static double cache[MAX_SIZE_ARRAYS];

  if (Ntable.Ntheta == 0)
  {
    log_fatal("Ntable.Ntheta not initialized");
    exit(1);
  }

  const int NSIZE = redshift.clustering_nbin;
  
  if (Pl == NULL || w_vec == NULL || fdiff(cache[3], Ntable.random))
  {
    if (Pl != NULL) free(Pl);
    Pl = (double**) malloc2d(Ntable.Ntheta, limits.LMAX);;

    if (w_vec != NULL) free(w_vec);
    w_vec = calloc1d(NSIZE*Ntable.Ntheta);

    double*** P = (double***) malloc3d(2, Ntable.Ntheta, limits.LMAX + 1);
    double** Pmin  = P[0]; double** Pmax  = P[1];

    double xmin[Ntable.Ntheta];
    double xmax[Ntable.Ntheta];
    for (int i=0; i<Ntable.Ntheta; i++)
    { // Cocoa: dont thread (init of static variables inside set_bin_average)
      bin_avg r = set_bin_average(i,0);
      xmin[i] = r.xmin;
      xmax[i] = r.xmax;
    }

    #pragma omp parallel for collapse(2)
    for (int i=0; i<Ntable.Ntheta; i++)
    {
      for (int l=0; l<limits.LMAX; l++)
      {
        bin_avg r = set_bin_average(i,l);
        Pmin[i][l] = r.Pmin;
        Pmax[i][l] = r.Pmax;
      }
    }

    const int lmin = 1;
    for (int i=0; i<Ntable.Ntheta; i++)
      for (int l=0; l<lmin; l++)
        Pl[i][l] = 0.0;
    
    #pragma omp parallel for collapse(2)
    for (int i=0; i<Ntable.Ntheta; i++)
    {
      for (int l=lmin; l<limits.LMAX; l++)
      {
        const double tmp = (1.0/(xmin[i] - xmax[i]))*(1.0 / (4.0 * M_PI));
        Pl[i][l] = tmp*(Pmin[i][l + 1] - Pmax[i][l + 1] 
                        - Pmin[i][l - 1] + Pmax[i][l - 1]);
      }
    }

    free(P);
  }

  if (fdiff(cache[0], cosmology.random) ||
      fdiff(cache[1], nuisance.random_photoz_clustering) ||
      fdiff(cache[2], redshift.random_clustering) ||
      fdiff(cache[3], Ntable.random) ||
      fdiff(cache[4], nuisance.random_galaxy_bias))
  { 
    double** Cl = (double**) malloc2d(NSIZE, limits.LMAX);

    const int lmin = 1;
    for (int i=0; i<NSIZE; i++)
      for (int l=0; l<lmin; l++)
        Cl[i][l] = 0.0;

    if (limber == 1)
    {
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wunused-variable"
      { // init static variables inside the C_XY_limber function
        const double l = limits.LMIN_tab+1;      
        double ini = C_gk_tomo_limber(l, 0)*beam_cmb(l)*w_pixel(l);
      }
      #pragma GCC diagnostic pop
      
      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<NSIZE; nz++)
        for (int l=lmin; l<limits.LMIN_tab; l++)
          Cl[nz][l] = C_gk_tomo_limber_nointerp(l, nz, 0)*beam_cmb(l)*w_pixel(l);

      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<NSIZE; nz++)
        for (int l=limits.LMIN_tab; l<limits.LMAX; l++)
          Cl[nz][l] = C_gk_tomo_limber(l, nz)*beam_cmb(l)*w_pixel(l);

      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<NSIZE; nz++)
      {
        for (int i=0; i<Ntable.Ntheta; i++)
        {
          w_vec[nz*Ntable.Ntheta+i] = 0;
          for (int l=lmin; l<limits.LMAX; l++)
            w_vec[nz*Ntable.Ntheta+i] += Pl[i][l]*Cl[nz][l];
        }
      }
    }
    else
    {
      log_fatal("NonLimber not implemented");
      exit(1);
    }

    free(Cl);

    cache[0] = cosmology.random;
    cache[1] = nuisance.random_photoz_clustering;
    cache[2] = redshift.random_clustering;
    cache[3] = Ntable.random;
    cache[4] = nuisance.random_galaxy_bias;
  }

  if (ni < -1 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d (max %d)", ni, redshift.clustering_nbin);
    exit(1); 
  } 
  if (nt < 0 || nt > Ntable.Ntheta - 1)
  {
    log_fatal("error in selecting bin number nt = %d (max %d)", nt, Ntable.Ntheta);
    exit(1); 
  } 
  
  const int q = ni * Ntable.Ntheta + nt;
  if (q < 0 || q > NSIZE*Ntable.Ntheta - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  }  
  return w_vec[q];
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double w_ks_tomo(const int nt, const int ni, const int limber)
{
  static double** Pl = NULL;
  static double* w_vec = NULL;
  static double cache[MAX_SIZE_ARRAYS];
  
  if (Ntable.Ntheta == 0)
  {
    log_fatal("Ntable.Ntheta not initialized");
    exit(1);
  }

  const int NSIZE = redshift.shear_nbin;

  if (Pl == NULL || w_vec == NULL || fdiff(cache[4], Ntable.random))
  {
    if (Pl != NULL) free(Pl);
    Pl = (double**) malloc2d(Ntable.Ntheta, limits.LMAX);;

    if (w_vec != NULL) free(w_vec);
    w_vec = calloc1d(NSIZE*Ntable.Ntheta);
    
    double*** P = (double***) malloc3d(2, Ntable.Ntheta, limits.LMAX + 1);
    double** Pmin  = P[0]; double** Pmax  = P[1];

    double xmin[Ntable.Ntheta];
    double xmax[Ntable.Ntheta];
    for (int i=0; i<Ntable.Ntheta; i++)
    { // Cocoa: dont thread (init of static variables inside set_bin_average)
      bin_avg r = set_bin_average(i,0);
      xmin[i] = r.xmin;
      xmax[i] = r.xmax;
    }

    #pragma omp parallel for collapse(2)
    for (int i=0; i<Ntable.Ntheta; i++)
    {
      for (int l=0; l<limits.LMAX; l++)
      {
        bin_avg r = set_bin_average(i,l);
        Pmin[i][l] = r.Pmin;
        Pmax[i][l] = r.Pmax;
      }
    }

    const int lmin = 1;
    for (int i=0; i<Ntable.Ntheta; i++)
      for (int l=0; l<lmin; l++)
        Pl[i][0] = 0.0;

    #pragma omp parallel for collapse(2)
    for (int i=0; i<Ntable.Ntheta; i++)
    {
      for (int l=lmin; l<limits.LMAX; l++)
      {
        Pl[i][l] = (2.*l+1)/(4.*M_PI*l*(l+1)*(xmin[i]-xmax[i]))
          *((l+2./(2*l+1.))*(Pmin[i][l-1]-Pmax[i][l-1])
          +(2-l)*(xmin[i]*Pmin[i][l]-xmax[i]*Pmax[i][l])
          -2./(2*l+1.)*(Pmin[i][l+1]-Pmax[i][l+1]));
      }
    }
    
    free(P);
  }

  if (fdiff(cache[0], cosmology.random) ||
      fdiff(cache[1], nuisance.random_photoz_shear) ||
      fdiff(cache[2], nuisance.random_ia) ||
      fdiff(cache[3], redshift.random_shear) || 
      fdiff(cache[4], Ntable.random))
  {
    double** Cl = (double**) malloc2d(NSIZE, limits.LMAX);

    const int lmin = 1;
    for (int i=0; i<NSIZE; i++)
      for (int l=0; l<lmin; l++)
        Cl[i][l] = 0.0;
      
    if (limber == 1)
    {      
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wunused-variable"
      { // init static variables inside the C_XY_limber function        
        const double l = limits.LMIN_tab+1;
        double ini = C_ks_tomo_limber(l, 0)*beam_cmb(l)*w_pixel(l);
      }
      #pragma GCC diagnostic pop
      
      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<redshift.shear_nbin; nz++)
        for (int l=lmin; l<limits.LMIN_tab; l++)
          Cl[nz][l] = C_ks_tomo_limber_nointerp(l, nz, 0)*beam_cmb(l)*w_pixel(l);
         
      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<NSIZE; nz++)
        for (int l=limits.LMIN_tab; l<limits.LMAX; l++)   
          Cl[nz][l] = C_ks_tomo_limber(l, nz)*beam_cmb(l)*w_pixel(l);
   
      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<NSIZE; nz++)
      {
        for (int i=0; i<Ntable.Ntheta; i++)
        {
          w_vec[nz*Ntable.Ntheta+i] = 0;
          for (int l=lmin; l<limits.LMAX; l++)
            w_vec[nz*Ntable.Ntheta+i] += Pl[i][l]*Cl[nz][l];
        }
      }
    } 
    else
    {
      log_fatal("NonLimber not implemented");
      exit(1);
    }

    free(Cl);
    cache[0] = cosmology.random;
    cache[1] = nuisance.random_photoz_shear;
    cache[2] = nuisance.random_ia;
    cache[3] = redshift.random_shear; 
    cache[4] = Ntable.random;
  }

  if (nt < 0 || nt > Ntable.Ntheta - 1)
  {
    log_fatal("error in selecting bin number nt = %d (max %d)", nt, Ntable.Ntheta);
    exit(1); 
  }
  if (ni < -1 || ni > redshift.shear_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d (max %d)", ni, redshift.shear_nbin);
    exit(1);
  }
  
  const int q = ni * Ntable.Ntheta + nt;
  if (q  < 0 || q > NSIZE*Ntable.Ntheta - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  }  
  
  return w_vec[q];
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Limber Approximation (Angular Power Spectrum)
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double int_for_C_ss_tomo_limber(double a, void* params)
{
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;

  const int n1 = (int) ar[0]; // first source bin 
  const int n2 = (int) ar[1]; // second source bin 
  if (n1 < 0 || 
      n1 > redshift.shear_nbin - 1 || 
      n2 < 0 || 
      n2 > redshift.shear_nbin - 1)
  {
    log_fatal("error in selecting bin number (ni, nj) = [%d,%d]", n1, n2);
    exit(1);
  }
  const double l = ar[2];
  const int EE = (int) ar[3];

  const double ell = l + 0.5;
  struct chis chidchi = chi_all(a);
  const double growfac_a = growfac(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  
  const double WK1 = W_kappa(a, fK, n1);
  const double WK2 = W_kappa(a, fK, n2);
  const double WS1 = W_source(a, n1, hoverh0);
  const double WS2 = W_source(a, n2, hoverh0);
  const double PK  = Pdelta(k,a);

  const double ell4 = ell*ell*ell*ell; // correction (1812.05995 eqs 74-79)
  const double ell_prefactor = l*(l - 1.)*(l + 1.)*(l + 2.)/ell4; 

  double ans;

  switch(like.IA_MODEL) 
  {
    case IA_MODEL_TATT:
    { 
      get_FPT_IA();

      const double lnk = log(k);
      const double g4 = growfac_a*growfac_a*growfac_a*growfac_a;
      
      double IA_AX[2];
      IA_A1_Z1Z2(a, growfac_a, n1, n2, IA_AX);
      const double C11 = IA_AX[0];
      const double C12 = IA_AX[1];
      IA_A2_Z1Z2(a, growfac_a, n1, n2, IA_AX);
      const double C21 = IA_AX[0];
      const double C22 = IA_AX[1];
      IA_BTA_Z1Z2(a, growfac_a, n1, n2, IA_AX);
      const double bta1 = IA_AX[0];
      const double bta2 = IA_AX[1];
      
      double lim[3];
      lim[0] = log(FPTIA.k_min);
      lim[1] = log(FPTIA.k_max);
      lim[2] = (lim[1] - lim[0])/FPTIA.N;

      if (EE == 1)
      {     
        const double tt = (lnk<lim[0] || lnk>lim[1]) ? 0.0 : 
          g4*interpol1d(FPTIA.tab[0], FPTIA.N, lim[0], lim[1], lim[2], lnk);
        
        const double ta_dE1 = (lnk<lim[0] || lnk>lim[1]) ? 0.0 : 
          g4*interpol1d(FPTIA.tab[2], FPTIA.N, lim[0], lim[1], lim[2], lnk);
        
        const double ta_dE2 = (lnk<lim[0] || lnk>lim[1]) ? 0.0 : 
          g4*interpol1d(FPTIA.tab[3], FPTIA.N, lim[0], lim[1], lim[2], lnk);
        
        const double ta = (lnk<lim[0] || lnk>lim[1]) ? 0.0 : 
          g4*interpol1d(FPTIA.tab[4], FPTIA.N, lim[0], lim[1], lim[2], lnk);
        
        const double mixA = (lnk<lim[0] || lnk>lim[1]) ? 0.0 : 
          g4*interpol1d (FPTIA.tab[6], FPTIA.N, lim[0], lim[1], lim[2], lnk);
        
        const double mixB = (lnk<lim[0] || lnk>lim[1]) ? 0.0 : 
          g4*interpol1d(FPTIA.tab[7], FPTIA.N, lim[0], lim[1], lim[2], lnk);
        
        const double mixEE = (lnk<lim[0] || lnk>lim[1]) ? 0.0 : 
          g4*interpol1d(FPTIA.tab[8], FPTIA.N, lim[0], lim[1], lim[2], lnk);
        
        ans = WK1*WK2*PK 
              - WS1*WK2*(C11*PK + C11*bta1*(ta_dE1+ta_dE2) - 5*C21*(mixA+mixB))
              - WS2*WK1*(C12*PK + C12*bta2*(ta_dE1+ta_dE2) - 5*C22*(mixA+mixB))
              + WS1*WS2*(C11*C12*PK 
                         + C11*C12*(bta1*bta2*ta + (bta1 + bta2)*(ta_dE1+ta_dE2))
                         - 5*(C11*C22 + C12*C21)*(mixA+mixB)
                         - 5*(C11*bta1*C22 + C12*bta2*C21)*mixEE
                         + 25*C21*C22*tt);
      }
      else  
      {        
        const double tt = (lnk<lim[0] || lnk>lim[1]) ? 0.0 : 
          g4*interpol1d(FPTIA.tab[1],FPTIA.N, lim[0], lim[1], lim[2], lnk);
        
        const double ta = (lnk<lim[0] || lnk>lim[1]) ? 0.0 : 
          g4*interpol1d(FPTIA.tab[5],FPTIA.N, lim[0], lim[1], lim[2], lnk);
        
        const double mix = (lnk<lim[0] || lnk>lim[1]) ? 0.0 : 
          g4*interpol1d(FPTIA.tab[9],FPTIA.N, lim[0], lim[1], lim[2], lnk);
        
        ans = WS1*WS2*(C11*C12*bta1*bta2*ta 
                       - 5*(C11*bta1*C22+C12*bta2*C21)*mix 
                       + 25*C21*C22*tt);
      }
      break;
    }
    case IA_MODEL_NLA:
    {
      if (EE == 1)
      { 
        double IA_A1[2];
        IA_A1_Z1Z2(a, growfac_a, n1, n2, IA_A1);
        const double C11 = IA_A1[0];
        const double C12 = IA_A1[1];
        ans =   WK1*WK2*PK 
              - WS1*WK2*C11*PK 
              - WS2*WK1*C12*PK
              + WS1*WS2*C11*C12*PK;
      }
      else
      {
        ans = 0.0;
      }
      break;
    }
    default:
    {
      log_fatal("like.IA_MODEL = %d not supported", like.IA_MODEL);
      exit(1);
    }
  }

  return ans*(chidchi.dchida/(fK*fK))*ell_prefactor;
}

double C_ss_tomo_limber_nointerp(
    const double l, 
    const int ni, 
    const int nj, 
    const int EE, 
    const int init
  )
{
  static double cache[MAX_SIZE_ARRAYS];
  static gsl_integration_glfixed_table* w = NULL;
  
  if (ni < -1 || 
      ni > redshift.shear_nbin -1 || 
      nj < -1 || 
      nj > redshift.shear_nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  if (w == NULL || fdiff(cache[0], Ntable.random))
  {
    const size_t szint = 60 + 50 * (Ntable.high_def_integration);
    if (w != NULL)  gsl_integration_glfixed_table_free(w);
    w = malloc_gslint_glfixed(szint);
    cache[0] = Ntable.random;
  }

  double ar[4] = {ni, nj, l, EE};
  const double amin = fmax(amin_source(ni), amin_source(nj));
  const double amax = fmin(amax_source(ni), amax_source(nj));;
  if (!(amin>0) || !(amin<1) || !(amax>0) || !(amax<1)) 
  {
    log_fatal("0 < amin/amax < 1 not true");
    exit(1);
  }

  double res = 0.0;
  if (init == 1)
    res = int_for_C_ss_tomo_limber(amin, (void*) ar);
  else
  {
    gsl_function F;
    F.params = (void*) ar;
    F.function = int_for_C_ss_tomo_limber;
    res = gsl_integration_glfixed(&F, amin, amax, w);
  }
  return res;     
}

double C_ss_tomo_limber(const double l, const int ni, const int nj, const int EE)
{
  static double cache[MAX_SIZE_ARRAYS];
  static double*** table;
  static double lim[3];
  static int nell;
  
  if (table == NULL || fdiff(cache[4], Ntable.random))
  {
    nell   = Ntable.N_ell;
    lim[0] = log(fmax(limits.LMIN_tab - 1., 1.0));
    lim[1] = log(limits.LMAX + 1);
    lim[2] = (lim[1] - lim[0]) / ((double) nell - 1.);

    if (table != NULL) free(table);
    table = (double***) malloc3d(2, tomo.shear_Npowerspectra, nell);
  }
  
  if (fdiff(cache[0], cosmology.random) ||
      fdiff(cache[1], nuisance.random_photoz_shear) ||
      fdiff(cache[2], nuisance.random_ia) ||
      fdiff(cache[3], redshift.random_shear) ||
      fdiff(cache[4], Ntable.random))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable" 
    {
      double ini = C_ss_tomo_limber_nointerp(exp(lim[0]),Z1(0),Z2(0),1,1);
      ini = C_ss_tomo_limber_nointerp(exp(lim[0]),Z1(0),Z2(0),0,1);
    }
    #pragma GCC diagnostic pop
    
    #pragma omp parallel for collapse(2)
    for (int k=0; k<tomo.shear_Npowerspectra; k++)
    {  
      for (int i=0; i<nell; i++)
      { 
        const double l = exp(lim[0]+i*lim[2]);
        const double Z1NZ = Z1(k);
        const double Z2NZ = Z2(k);
        table[0][k][i] = C_ss_tomo_limber_nointerp(l,Z1NZ,Z2NZ,1,0);
        table[1][k][i] = C_ss_tomo_limber_nointerp(l,Z1NZ,Z2NZ,0,0);
      }
    }
    cache[0] = cosmology.random;
    cache[1] = nuisance.random_photoz_shear;
    cache[2] = nuisance.random_ia;
    cache[3] = redshift.random_shear;
    cache[4] = Ntable.random;
  }
  
  if (ni < 0 || 
      ni > redshift.shear_nbin - 1 || 
      nj < 0 ||
      nj > redshift.shear_nbin - 1)
  {
    log_fatal("error in selecting bin number (ni, nj) = [%d,%d]", ni, nj);
    exit(1);
  }
  
  const double lnl = log(l);
  if (lnl < lim[0])
    log_warn("l = %e < l_min = %e. Extrapolation adopted", l, exp(lim[0]));
  if (lnl > lim[1])
    log_warn("l = %e > l_max = %e. Extrapolation adopted", l, exp(lim[1]));
  
  const int q = N_shear(ni, nj);
  if (q < 0 || q > tomo.shear_Npowerspectra - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  }
  return  interpol1d((EE == 1) ? table[0][q] : table[1][q], 
                     nell, lim[0], lim[1], lim[2], lnl);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_gs_tomo_limber(double a, void* params)
{
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;
  const int nl = (int) ar[0];
  const int ns = (int) ar[1];
  if (nl < 0 || 
      nl > redshift.clustering_nbin - 1 || 
      ns < 0 || 
      ns > redshift.shear_nbin - 1)
  {
    log_fatal("error in selecting bin number (nl, ns) = [%d,%d]", nl, ns);
    exit(1);
  }
  const double l = ar[2];
  const int nonlinear_bias = ar[3];
  
  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double g4 = growfac_a*growfac_a*growfac_a*growfac_a;
  const double ell = l + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double z = 1.0/a - 1.0;
  const double PK = Pdelta(k,a);

  const double b1 = gb1(z, nl);
  const double bmag = gbmag(z, nl);

  const double WK = W_kappa(a, fK, ns);
  const double WGAL = W_gal(a, nl, hoverh0);
  const double WMAG = W_mag(a, fK, nl);
  const double WS   = W_source(a, ns, hoverh0);

  const double ell_prefactor = l*(l + 1.)/(ell*ell); // correction (1812.05995 eqs 74-79)
  const double tmp = (l - 1.)*l*(l + 1.)*(l + 2.);   // correction (1812.05995 eqs 74-79)
  const double ell_prefactor2 = (tmp > 0) ? sqrt(tmp)/(ell*ell) : 0.0;

  double ans;

  switch(like.IA_MODEL)
  {
    case IA_MODEL_TATT:
    {
      if (include_HOD_GX == 1)
      {
        log_fatal("HOD NOT IMPLEMENTED");
        exit(1);
      }

      get_FPT_IA();
      const double lnk = log(k);
      double lim[3];
      lim[0] = log(FPTIA.k_min);
      lim[1] = log(FPTIA.k_max);
      lim[2] = (lim[1] - lim[0])/FPTIA.N;

      const double mixA = (lnk<lim[0] || lnk>lim[1]) ? 0.0 : 
        g4*interpol1d(FPTIA.tab[6], FPTIA.N, lim[0], lim[1], lim[2], lnk);
      
      const double mixB = (lnk<lim[0] || lnk>lim[1]) ? 0.0 :
        g4*interpol1d(FPTIA.tab[7], FPTIA.N, lim[0], lim[1], lim[2], lnk);
      
      const double ta_dE1 = (lnk<lim[0] || lnk>lim[1]) ? 0.0 :
        g4*interpol1d(FPTIA.tab[2], FPTIA.N, lim[0], lim[1], lim[2], lnk);
      
      const double ta_dE2 = (lnk<lim[0] || lnk>lim[1]) ? 0.0 :
        g4*interpol1d(FPTIA.tab[3], FPTIA.N, lim[0], lim[1], lim[2], lnk);

      double WRSD = 0.0;
      if (include_RSD_GS == 1)
      {
        const double chi_0 = f_K(ell/k);
        const double chi_1 = f_K((ell+1.)/k);
        const double a_0 = a_chi(chi_0);
        const double a_1 = a_chi(chi_1);
        WRSD = W_RSD(ell, a_0, a_1, nl);
      }

      double oneloop = 0.0;
      if (1 == nonlinear_bias)
      { 
        get_FPT_bias();
        lim[0] = log(FPTbias.k_min);
        lim[1] = log(FPTbias.k_max);
        lim[2] = (lim[1] - lim[0])/FPTbias.N;

        const double d1d2 = (lnk<lim[0] || lnk>lim[1]) ? 0.0 :
          interpol1d(FPTbias.tab[0], FPTbias.N, lim[0], lim[1], lim[2], lnk);
        
        const double d1s2 = (lnk<lim[0] || lnk>lim[1]) ? 0.0 :
          interpol1d(FPTbias.tab[2], FPTbias.N, lim[0], lim[1], lim[2], lnk);
        
        const double d1d3 = (lnk<lim[0] || lnk>lim[1]) ? 0.0 :
          interpol1d(tab_d1d3, FPTbias.N, lim[0], lim[1], lim[2], lnk);

        const double b2 = gb2(z, nl);
        const double bs2 = gbs2(z, nl);
        const double b3 = gb3(z, nl);

        oneloop = 0.5*g4*(b2 * d1d2 + bs2 * d1s2 + b3 * d1d3);
      }

      const double C1ZS  = IA_A1_Z1(a, growfac_a, ns);
      const double btazs = IA_BTA_Z1(a, growfac_a, ns);
      const double C2ZS  = IA_A2_Z1(a, growfac_a, ns);

      ans = (WK-WS*C1ZS)*((WGAL*b1+WMAG*ell_prefactor*bmag+WRSD)*PK+WGAL*oneloop) 
            - WS*(WGAL*b1 + WMAG*ell_prefactor*bmag)*(C1ZS*btazs*(ta_dE1+ta_dE2) 
                                                      - 5*C2ZS*(mixA+mixB));
      break;
    }
    case IA_MODEL_NLA:
    {
      if (include_HOD_GX == 1)
      {
        log_fatal("HOD NOT IMPLEMENTED");
        exit(1);
      }

      double WRSD = 0.0;
      if (include_RSD_GS == 1)
      {
        const double chi_0 = f_K(ell/k);
        const double chi_1 = f_K((ell+1.)/k);
        const double a_0 = a_chi(chi_0);
        const double a_1 = a_chi(chi_1);
        WRSD = W_RSD(ell, a_0, a_1, nl);
      }

      double oneloop = 0.0;
      if (1 == nonlinear_bias)
      {
        get_FPT_bias();
        const double lnk = log(k);
        double lim[3];
        lim[0] = log(FPTbias.k_min);
        lim[1] = log(FPTbias.k_max);
        lim[2] = (lim[1] - lim[0])/FPTbias.N;

        const double d1d2 = (lnk<lim[0] || lnk>lim[1]) ? 0.0 :
          interpol1d(FPTbias.tab[0], FPTbias.N, lim[0], lim[1], lim[2], lnk);
        
        const double d1s2 = (lnk<lim[0] || lnk>lim[1]) ? 0.0 :
          interpol1d(FPTbias.tab[2], FPTbias.N, lim[0], lim[1], lim[2], lnk);
        
        const double d1d3 = (lnk<lim[0] || lnk>lim[1]) ? 0.0 :
          interpol1d(tab_d1d3, FPTbias.N, lim[0], lim[1], lim[2], lnk);

        const double b2 = gb2(z, nl);
        const double bs2 = gbs2(z, nl);
        const double b3 = gb3(z, nl);

        oneloop = 0.5*g4*(b2*d1d2 + bs2*d1s2 + b3*d1d3);
      }
      
      const double C1ZS = IA_A1_Z1(a, growfac_a, ns);

      ans = (WK-WS*C1ZS)*((WGAL*b1+WMAG*ell_prefactor*bmag+WRSD)*PK+WGAL*oneloop);
      break;
    }
    default:
    {
      log_fatal("like.IA_MODEL = %d not supported", like.IA_MODEL);
      exit(1);
    }
  }
  return ans*(chidchi.dchida/(fK*fK))*ell_prefactor2;
}

double C_gs_tomo_limber_nointerp(
    const double l, 
    const int nl, 
    const int ns,
    const int init
  )
{
  static double cache[MAX_SIZE_ARRAYS];
  static gsl_integration_glfixed_table* w = NULL;
  
  if (nl < -1 || 
      nl > redshift.clustering_nbin -1 || 
      ns < -1 || 
      ns > redshift.shear_nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", nl, ns);
    exit(1);
  }

  if (w == NULL || fdiff(cache[0], Ntable.random))
  {
    const size_t szint = 125 + 50 * (Ntable.high_def_integration);
    if (w != NULL)  gsl_integration_glfixed_table_free(w);
    w = malloc_gslint_glfixed(szint);
    cache[0] = Ntable.random;
  }

  double ar[4] = {(double) nl, (double) ns, l, has_b2_galaxies()};
  
  const double amin = amin_lens(nl);
  const double amax = amax_lens(nl);
  
  if (!(amin>0) || !(amin<1) || !(amax>0) || !(amax<1)) 
  {
    log_fatal("0 < amin/amax < 1 not true");
    exit(1);
  }
  if (!(amin < amax))
  {
    log_fatal("amin < amax not true");
    exit(1);
  }

  double res;
  if (init == 1)
    res = int_for_C_gs_tomo_limber(amin, (void*) ar);
  else
  {    
    gsl_function F;
    F.params = (void*) ar;
    F.function = int_for_C_gs_tomo_limber;
    res = gsl_integration_glfixed(&F, amin, amax, w);
  }
  return res;
}

double C_gs_tomo_limber(const double l, const int ni, const int nj)
{
  static double cache[MAX_SIZE_ARRAYS];
  static double** table = NULL;
  static int nell;
  static double lim[3];

  if (table == NULL || fdiff(cache[6], Ntable.random)) 
  {
    nell   = Ntable.N_ell;
    lim[0] = log(fmax(limits.LMIN_tab, 1.0));
    lim[1] = log(limits.LMAX + 1);
    lim[2] = (lim[1] - lim[0]) / ((double) nell - 1.0);

    if (table != NULL) free(table);
    table = (double**) malloc2d(tomo.ggl_Npowerspectra, nell);
  }

  if (fdiff(cache[0], cosmology.random) ||
      fdiff(cache[1], nuisance.random_photoz_shear) ||
      fdiff(cache[2], nuisance.random_photoz_clustering) ||
      fdiff(cache[3], nuisance.random_ia) ||
      fdiff(cache[4], redshift.random_shear) ||
      fdiff(cache[5], redshift.random_clustering) ||
      fdiff(cache[6], Ntable.random) ||
      fdiff(cache[7], nuisance.random_galaxy_bias))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    for (int k=0; k<tomo.ggl_Npowerspectra; k++)
    { // init static variables inside the C_XY_limber_nointerp function     
      double ini = C_gs_tomo_limber_nointerp(exp(lim[0]), ZL(k), ZS(k), 1);
    }
    #pragma GCC diagnostic pop

    #pragma omp parallel for collapse(2)
    for (int k=0; k<tomo.ggl_Npowerspectra; k++)
      for (int i=0; i<nell; i++)      
        if (test_zoverlap(ZL(k), ZS(k)))
          table[k][i] = C_gs_tomo_limber_nointerp(exp(lim[0]+i*lim[2]),ZL(k),ZS(k),0);
        else
          table[k][i] = 0.0;

    cache[0] = cosmology.random;
    cache[1] = nuisance.random_photoz_shear;
    cache[2] = nuisance.random_photoz_clustering;
    cache[3] = nuisance.random_ia;
    cache[4] = redshift.random_shear;
    cache[5] = redshift.random_clustering;
    cache[6] = Ntable.random;
    cache[7] = nuisance.random_galaxy_bias;
  }
  
  const int q = N_ggl(ni, nj);
  if (q < 0 || q > tomo.ggl_Npowerspectra - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  }

  double res = 0.0;
  if (test_zoverlap(ni, nj))
  {
    const double lnl = log(l);
    if (lnl < lim[0])
      log_warn("l = %e < l_min = %e. Extrapolation adopted", l, exp(lim[0]));
    if (lnl > lim[1])
      log_warn("l = %e > l_max = %e. Extrapolation adopted", l, exp(lim[1]));
    res = interpol1d(table[q], nell, lim[0], lim[1], lim[2], lnl);
  }
  return res;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_gg_tomo_limber(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;

  const int ni = (int) ar[0];
  const int nj = (int) ar[1];
  const double l  = ar[2];
  const int use_linear_ps = (int) ar[3];
  const int nonlinear_bias = ar[4];


  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double ell = l + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell / fK;
  const double z = 1.0/a - 1.0;

  const double b1i   = gb1(z, ni);
  const double bmagi = gbmag(z, ni);
  const double WGALi = W_gal(a, ni, hoverh0);
  const double WMAGi = W_mag(a, fK, ni);
  
  const double b1j   = b1i;    // We assume ni = nj;
  const double bmagj = bmagi;
  const double WGALj = WGALi;
  const double WMAGj = WMAGi;

  const double ell_prefactor = l*(l+1.)/(ell*ell); // prefactor correction (1812.05995 eqs 74-79)
  
  double res = 1.0;
  
  if (include_HOD_GX == 1)
  {
    if (include_RSD_GG == 1)
    {
      log_fatal("RSD not implemented with (HOD = TRUE)");
      exit(1);
    }
    else
      res *= WGALi*WGALj;
    res *= p_gg(k, a, ni, nj);
  }
  else
  {
    if(include_RSD_GG == 1)
    {
      const double chi_a_min = chi(limits.a_min);
      const double chi_0 = f_K(ell/k);
      const double chi_1 = f_K((ell + 1.0)/k);
      if (chi_1 > chi_a_min)  return 0;
      
      const double a_0 = a_chi(chi_0);
      const double a_1 = a_chi(chi_1);
      const double WRSDi =  W_RSD(ell, a_0, a_1, ni);

      res *= (WGALi*b1i + WMAGi*ell_prefactor*bmagi + WRSDi);
      res *= (WGALi*b1i + WMAGi*ell_prefactor*bmagi + WRSDi);
    }
    else
    {
      res *= (WGALi*b1i + WMAGi*ell_prefactor*bmagi);
      res *= (WGALi*b1i + WMAGi*ell_prefactor*bmagi);
    }
    res *= (use_linear_ps ? p_lin(k,a) : Pdelta(k,a));
  }

  double oneloop = 0.0;
  if (1 == nonlinear_bias && 0 == use_linear_ps)
  {
    get_FPT_bias();
    const double lnk = log(k);
    double lim[3];
    lim[0] = log(FPTbias.k_min);
    lim[1] = log(FPTbias.k_max);
    lim[2] = (lim[1] - lim[0])/FPTbias.N;

    const double d1d2 = (lnk<lim[0] || lnk>lim[1]) ? 0.0 :
      interpol1d(FPTbias.tab[0], FPTbias.N, lim[0], lim[1], lim[2], lnk);
    
    const double d2d2 = (lnk<lim[0] || lnk>lim[1]) ? 0.0 :
      interpol1d(FPTbias.tab[1], FPTbias.N, lim[0], lim[1], lim[2], lnk);
    
    const double d1s2 = (lnk<lim[0] || lnk>lim[1]) ? 0.0 :
      interpol1d(FPTbias.tab[2], FPTbias.N, lim[0], lim[1], lim[2], lnk);
    
    const double d2s2 = (lnk<lim[0] || lnk>lim[1]) ? 0.0 :
      interpol1d(FPTbias.tab[3], FPTbias.N, lim[0], lim[1], lim[2], lnk);
    
    const double s2s2 = (lnk<lim[0] || lnk>lim[1]) ? 0.0 :
      interpol1d(FPTbias.tab[4], FPTbias.N, lim[0], lim[1], lim[2], lnk);
    
    const double d1d3 = (lnk<lim[0] || lnk>lim[1]) ? 0.0 :
      interpol1d(tab_d1d3, FPTbias.N, lim[0], lim[1], lim[2], lnk);
    
    const double s4 = 0.; // PT_sigma4(k);

    const double growfac_a = growfac(a);
    const double g4 = growfac_a*growfac_a*growfac_a*growfac_a;
    const double b2 = gb2(z, ni);
    const double bs2 = gbs2(z, ni);
    const double b3 = gb3(z, ni);
    
    oneloop *= WGALi*WGALi;
    oneloop *= g4*(b1i*b2*d1d2 + 0.25*b2*b2 * (d2d2 - 2.*s4) +
      b1i*bs2*d1s2 + 0.5*b2*bs2 * (d2s2 - 4. / 3.*s4) +
      0.25*bs2*bs2* (s2s2 - 8. / 9. * s4) + b1i*b3*d1d3);
  }
  return (res +  oneloop)*chidchi.dchida/(fK*fK);
}

double C_gg_tomo_limber_linpsopt_nointerp(
    const double l, 
    const int ni, 
    const int nj,
    const int use_linear_ps,
    const int init
  ) // We need that for limber calculation
{
  static double cache[MAX_SIZE_ARRAYS];
  static gsl_integration_glfixed_table* w = NULL;
  
  if (ni < 0 || 
      ni > redshift.clustering_nbin - 1 || 
      nj < 0 || 
      nj > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number (ni, nj) = [%d,%d]", ni, nj);
    exit(1);
  }
  if (ni != nj)
  {
    log_fatal("cross-tomography (ni,nj) = (%d,%d) bins not supported", ni, nj);
    exit(1);
  }

  if (w == NULL || fdiff(cache[0], Ntable.random))
  {
    const size_t szint = 200 + 50 * (Ntable.high_def_integration);
    if (w != NULL)  gsl_integration_glfixed_table_free(w);
    w = malloc_gslint_glfixed(szint);
    cache[0] = Ntable.random;
  }

  double ar[5] = {ni, nj, l, use_linear_ps, has_b2_galaxies()};
  
  const double amin = amin_lens(ni);
  const double amax = amax_lens(ni);
  
  if (!(amin>0) || !(amin<1) || !(amax>0) || !(amax<1)) 
  {
    log_fatal("0 < amin/amax < 1 not true");
    exit(1);
  }
  if (!(amin < amax))
  {
    log_fatal("amin < amax not true");
    exit(1);
  }

  double res = 0.0;
  if (init == 1)
    int_for_C_gg_tomo_limber(amin, (void*) ar);
  else
  {
    gsl_function F;
    F.params = (void*) ar;
    F.function = int_for_C_gg_tomo_limber;
    res = gsl_integration_glfixed(&F, amin, amax, w);
  }
  return res;
}

double C_gg_tomo_limber_nointerp(
    const double l, 
    const int ni, 
    const int nj,
    const int init
  )
{
  return  C_gg_tomo_limber_linpsopt_nointerp(l, ni, nj, 0, init);
}

double C_gg_tomo_limber(const double l, const int ni, const int nj)
{ // cross redshift bin not supported
  static double cache[MAX_SIZE_ARRAYS];
  static double** table = NULL;
  static int nell;
  static int NSIZE;
  static double lnlmin;
  static double lnlmax;
  static double dlnl;
  
  if (table == NULL || fdiff(cache[3], Ntable.random))
  {
    nell   = Ntable.N_ell;
    NSIZE  = redshift.clustering_nbin;
    lnlmin = log(fmax(limits.LMIN_tab, 1.0));
    lnlmax = log(limits.LMAX + 1);
    dlnl   = (lnlmax - lnlmin) / ((double) nell - 1.0);

    if (table != NULL) free(table);
    table = (double**) malloc2d(NSIZE, nell);
  }

  if (fdiff(cache[0], cosmology.random) ||
      fdiff(cache[1], nuisance.random_photoz_clustering) ||
      fdiff(cache[2], redshift.random_clustering) ||
      fdiff(cache[3], Ntable.random) ||
      fdiff(cache[4], nuisance.random_galaxy_bias))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    for (int k=0; k<NSIZE; k++) 
    { // init static variables inside the C_XY_limber_nointerp function
      const int init = 1; // TRUE
      double ini = C_gg_tomo_limber_nointerp(exp(lnlmin), k, k, 1);
    }
    #pragma GCC diagnostic pop
    
    #pragma omp parallel for collapse(2)
    for (int k=0; k<NSIZE; k++)  
      for (int p=0; p<nell; p++)
        table[k][p] = C_gg_tomo_limber_nointerp(exp(lnlmin+p*dlnl), k, k, 0);

    cache[0] = cosmology.random;
    cache[1] = nuisance.random_photoz_clustering;
    cache[2] = redshift.random_clustering;
    cache[3] = Ntable.random;
    cache[4] = nuisance.random_galaxy_bias;
  }
  
  if (ni < 0 || 
      ni > redshift.clustering_nbin - 1 || 
      nj < 0 || 
      nj > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number (ni, nj) = [%d,%d]", ni, nj);
    exit(1);
  }

  if (ni != nj) 
  {
    log_fatal("cross-tomography not supported");
    exit(1);
  }
  
  const double lnl = log(l);
  if (lnl < lnlmin)
    log_warn("l = %e < l_min = %e. Extrapolation adopted", l, exp(lnlmin));
  if (lnl > lnlmax)
    log_warn("l = %e > l_max = %e. Extrapolation adopted", l, exp(lnlmax));
  
  const int q = ni; // cross redshift bin not supported; not using N_CL(ni, nj)
  if (q < 0 || q > NSIZE - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  }  
  return interpol1d(table[q], nell, lnlmin, lnlmax, dlnl, lnl);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_gk_tomo_limber(double a, void* params)
{
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;

  const int nl = (int) ar[0];
  if (nl < 0 || nl > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", nl);
    exit(1);
  }
  const double l = ar[1];
  const int nonlinear = ar[2];

  const double ell = l + 0.5;
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double z = 1./a - 1.;

  const double b1   = gb1(z, nl);
  const double bmag = gbmag(z, nl);

  const double WK = W_k(a, fK);
  const double WGAL = W_gal(a, nl, hoverh0);
  const double WMAG = W_mag(a, fK, nl);

  const double ell_prefactor = l*(l + 1.)/(ell*ell); // prefactor correction (1812.05995 eqs 74-79)

  double res = WK; 

  if (include_HOD_GX == 1)
  {
    if (include_RSD_GK == 1)
    {
      log_fatal("RSD not implemented with (HOD = TRUE)");
      exit(1);
    }
    else
    { 
      res *= WGAL;
    }
    res *= p_gm(k, a, nl);
  }
  else
  {
    if (include_RSD_GK == 1)
    {
      const double chi_0 = f_K(ell/k);
      const double chi_1 = f_K((ell+1.)/k);
      const double a_0 = a_chi(chi_0);
      const double a_1 = a_chi(chi_1);
      const double WRSD = W_RSD(ell, a_0, a_1, nl);

      res *= WGAL*b1 + WMAG*ell_prefactor*bmag + WRSD;
    }
    else
    {
      res *= WGAL*b1 + WMAG*ell_prefactor*bmag;
    }
    const double PK = Pdelta(k,a);
    res *= PK;
  }

  double oneloop = WK;
  if (1 == nonlinear)
  {
    get_FPT_bias();
    const double growfac_a = growfac(a);
    const double g4 = growfac_a*growfac_a*growfac_a*growfac_a;

    const double lnk = log(k);
    double lim[3];
    lim[0] = log(FPTbias.k_min);
    lim[1] = log(FPTbias.k_max);
    lim[2] = (lim[1] - lim[0])/FPTbias.N;

    const double d1d2 = (lnk<lim[0] || lnk>lim[1]) ? 0.0 :
      interpol1d(FPTbias.tab[0], FPTbias.N, lim[0], lim[1], lim[2], lnk);
    
    const double d1s2 = (lnk<lim[0] || lnk>lim[1]) ? 0.0 :
      interpol1d(FPTbias.tab[0], FPTbias.N, lim[0], lim[1], lim[2], lnk);

    const double b2 = gb2(z, nl);
    const double bs2 = gbs2(z, nl);
    
    oneloop *= WGAL;
    oneloop *= g4*(0.5*b2*d1d2 + 0.5*bs2*d1s2);
  }
  return ((res + oneloop)*chidchi.dchida/(fK*fK))*ell_prefactor;
}

double C_gk_tomo_limber_nointerp(const double l, const int ni, const int init)
{
  static double cache[MAX_SIZE_ARRAYS];
  static gsl_integration_glfixed_table* w = NULL;
  
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }

  if (w == NULL || fdiff(cache[0], Ntable.random))
  {
    const size_t szint = 200 + 50 * (Ntable.high_def_integration);
    if (w != NULL)  gsl_integration_glfixed_table_free(w);
    w = malloc_gslint_glfixed(szint);
    cache[0] = Ntable.random;
  }

  double ar[3] = {(double) ni, l, has_b2_galaxies()};
  
  const double amin = amin_lens(ni);
  const double amax = amax_lens(ni);
  if (!(amin>0) || !(amin<1) || !(amax>0) || !(amax<1)) 
  {
    log_fatal("0 < amin/amax < 1 not true");
    exit(1);
  }
  if (!(amin < amax))
  {
    log_fatal("amin < amax not true");
    exit(1);
  }

  double res = 0.0;
  if (init == 1)
    res = int_for_C_gk_tomo_limber(amin, (void*) ar);
  else
  {
    gsl_function F;
    F.params = (void*) ar;
    F.function = int_for_C_gk_tomo_limber;
    res = gsl_integration_glfixed(&F, amin, amax, w);
  }
  return res;
}

double C_gk_tomo_limber(const double l, const int ni)
{
  static double cache[MAX_SIZE_ARRAYS];
  static double** table = NULL;
  static double lim[3];

  if (table == NULL || fdiff(cache[3], Ntable.random))
  {
    lim[0] = log(fmax(limits.LMIN_tab, 1.0));
    lim[1] = log(limits.LMAX + 1);
    lim[2] = (lim[1] - lim[0])/((double) Ntable.N_ell - 1.0);

    if (table != NULL) free(table);
    table = (double**) malloc2d(redshift.clustering_nbin, Ntable.N_ell);
  }

  if (fdiff(cache[0], cosmology.random) ||
      fdiff(cache[1], nuisance.random_photoz_clustering) ||
      fdiff(cache[2], redshift.random_clustering) ||
      fdiff(cache[3], Ntable.random) ||
      fdiff(cache[4], nuisance.random_galaxy_bias))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    { // init static variables inside the C_XY_limber_nointerp function      
      double ini = C_gk_tomo_limber_nointerp(exp(lim[0]), 0, 1);
    }
    #pragma GCC diagnostic pop
    
    #pragma omp parallel for collapse(2)
    for (int k=0; k<redshift.clustering_nbin; k++)
      for (int i=0; i<Ntable.N_ell; i++)  
        table[k][i] = C_gk_tomo_limber_nointerp(exp(lim[0] + i*lim[2]), k, 0);
    
    cache[0] = cosmology.random;
    cache[1] = nuisance.random_photoz_clustering;
    cache[2] = redshift.random_clustering;
    cache[3] = Ntable.random;
    cache[4] = nuisance.random_galaxy_bias;
  }
  
  if (ni < -1 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }

  const double lnl = log(l);
  if (lnl < lim[0])
    log_warn("l = %e < l_min = %e. Extrapolation adopted", l, exp(lim[0]));
  if (lnl > lim[1])
    log_warn("l = %e > l_max = %e. Extrapolation adopted", l, exp(lim[1]));
  
  const int q =  ni; 
  if (q < 0 || q > redshift.clustering_nbin - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  }
  return interpol1d(table[q], Ntable.N_ell, lim[0], lim[1], lim[2], lnl);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_ks_tomo_limber(double a, void* params)
{
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }

  double* ar = (double*) params;
  const int ni = (int) ar[0];
  if (ni < -1 || ni > redshift.shear_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }
  const double l = ar[1];  
  const double ell = l + 0.5;  
  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double PK = Pdelta(k,a);

  const double WK1 = W_kappa(a, fK, ni);
  const double WK2 = W_k(a, fK);

  const double ell_prefactor1 = l*(l + 1.)/(ell*ell); // prefactor correction (1812.05995 eqs 74-79)
  const double tmp = (l - 1.)*l*(l + 1.)*(l + 2.);    // prefactor correction (1812.05995 eqs 74-79)
  const double ell_prefactor2 = (tmp > 0) ? sqrt(tmp)/(ell*ell) : 0.0; 

  const double A_Z1 = IA_A1_Z1(a, growfac_a, ni);
  const double WS1  = W_source(a, ni, hoverh0) * A_Z1;

  const double res = (WK1 - WS1)*WK2;

  return (res*PK*chidchi.dchida/(fK*fK))*ell_prefactor1*ell_prefactor2;
}

double C_ks_tomo_limber_nointerp(const double l, const int ni, const int init)
{
  static double cache[MAX_SIZE_ARRAYS];
  static gsl_integration_glfixed_table* w = NULL;  

  if (ni < 0 || ni > redshift.shear_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }

  if (w == NULL || fdiff(cache[0], Ntable.random))
  {
    const size_t szint = 200 + 50 * (Ntable.high_def_integration);
    if (w != NULL)  gsl_integration_glfixed_table_free(w);
    w = malloc_gslint_glfixed(szint);
    cache[0] = Ntable.random;
  }
  
  double ar[3] = {(double) ni, l};
  const double amin = amin_source(ni);
  const double amax = amax_source(ni);
  if (!(amin>0) || !(amin<1) || !(amax>0) || !(amax<1)) 
  {
    log_fatal("0 < amin/amax < 1 not true");
    exit(1);
  }
 
  double res = 0.0;
  if (init == 1)
    res = int_for_C_ks_tomo_limber(amin, (void*) ar);
  else
  {
    gsl_function F;
    F.params = (void*) ar;
    F.function = int_for_C_ks_tomo_limber;
    res = gsl_integration_glfixed(&F, amin, amax, w);
  }
  return res;  
}

double C_ks_tomo_limber(double l, int ni)
{
  static double cache[MAX_SIZE_ARRAYS];
  static double** table = NULL;
  static double lim[3];

  if (table == NULL || fdiff(cache[4], Ntable.random))
  {
    if (table != NULL) free(table);
    table = (double**) malloc2d(redshift.shear_nbin, Ntable.N_ell);

    lim[0] = log(fmax(limits.LMIN_tab, 1.0));
    lim[1] = log(limits.LMAX + 1);
    lim[2] = (lim[1] - lim[0])/((double) Ntable.N_ell - 1.0);
  }

  if (fdiff(cache[0], cosmology.random) ||
      fdiff(cache[1], nuisance.random_photoz_shear) ||
      fdiff(cache[2], nuisance.random_ia) ||
      fdiff(cache[3], redshift.random_shear) ||
      fdiff(cache[4], Ntable.random))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    { // init static variables inside the C_XY_limber_nointerp function
      double ini = C_ks_tomo_limber_nointerp(exp(lim[0]), 0, 1);
    }
    #pragma GCC diagnostic pop
    
    #pragma omp parallel for collapse(2)
    for (int k=0; k<redshift.shear_nbin; k++)
      for (int i=0; i<Ntable.N_ell; i++)
        table[k][i] = C_ks_tomo_limber_nointerp(exp(lim[0] + i*lim[2]), k, 0);
 
    cache[0] = cosmology.random;
    cache[1] = nuisance.random_photoz_shear;
    cache[2] = nuisance.random_ia;
    cache[3] = redshift.random_shear;
    cache[4] = Ntable.random;
  } 
  
  if (ni < 0 || ni > redshift.shear_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d (max %d)", ni, redshift.shear_nbin);
    exit(1);
  }

  const double lnl = log(l);
  if (lnl < lim[0])
    log_warn("l = %e < l_min = %e. Extrapolation adopted", l, exp(lim[0]));
  if (lnl > lim[1])
    log_warn("l = %e > l_max = %e. Extrapolation adopted", l, exp(lim[1]));
  
  const int q =  ni; 
  if (q < 0 || q > redshift.shear_nbin - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  }
  return interpol1d(table[q], Ntable.N_ell, lim[0], lim[1], lim[2], lnl);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_kk_limber(double a, void* params)
{
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  
  double* ar = (double*) params;
  const double l = ar[0];
  
  struct chis chidchi = chi_all(a);
  const double ell = l + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double WK = W_k(a, fK);
  const double PK = Pdelta(k,a);
  
  const double ell_prefactor = l*(l + 1.0)/(ell*ell);  // prefac correction (1812.05995 eqs 74-79)

  return WK*WK*PK*(chidchi.dchida/(fK*fK))*ell_prefactor*ell_prefactor;
}

double C_kk_limber_nointerp(const double l, const int init)
{
  static double cache[MAX_SIZE_ARRAYS];
  static gsl_integration_glfixed_table* w = NULL;
  
  if (w == NULL || fdiff(cache[0], Ntable.random))
  {
    const size_t szint = 200 + 50 * (Ntable.high_def_integration);
    if (w != NULL)  gsl_integration_glfixed_table_free(w);
    w = malloc_gslint_glfixed(szint);
    cache[0] = Ntable.random;
  }

  double ar[1] = {l};
  const double amin = limits.a_min*(1. + 1.e-5);
  const double amax = 0.99999;
  
  double res = 0.0;
  if (init == 1)
    res = int_for_C_kk_limber(amin, (void*) ar);
  else
  {
    gsl_function F;
    F.params = (void*) ar;
    F.function = int_for_C_kk_limber;
    res = gsl_integration_glfixed(&F, amin, amax, w);
  }
  return res;
}

double C_kk_limber(const double l)
{
  static double cache[MAX_SIZE_ARRAYS];
  static double* table = NULL;
  static double lim[3];

  if (table == NULL || fdiff(cache[1], Ntable.random))
  {
    lim[0] = log(fmax(limits.LMIN_tab, 1.0));
    lim[1] = log(limits.LMAX + 1);
    lim[2] = (lim[1] - lim[0])/((double) Ntable.N_ell - 1.0);

    if (table != NULL) free(table);
    table = (double*) malloc1d(Ntable.N_ell);
  }

  if (fdiff(cache[0], cosmology.random) || fdiff(cache[1], Ntable.random))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    { // init static variables inside the C_XY_limber_nointerp function
      double ini = C_kk_limber_nointerp(exp(lim[0]), 1);
    }
    #pragma GCC diagnostic pop
    
    #pragma omp parallel for
    for (int i=0; i<Ntable.N_ell; i++)
      table[i] = C_kk_limber_nointerp(exp(lim[0] + i*lim[2]), 0);
  
    cache[0] = cosmology.random;
    cache[1] = Ntable.random;
  }

  const double lnl = log(l);
  if (lnl < lim[0])
    log_warn("l = %e < l_min = %e. Extrapolation adopted", l, exp(lim[0]));
  if (lnl > lim[1])
    log_warn("l = %e > l_max = %e. Extrapolation adopted", l, exp(lim[1]));
  
  return interpol1d(table, Ntable.N_ell, lim[0], lim[1], lim[2], lnl);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_gy_tomo_limber(double a, void* params)
{
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;
  const int nl = (int) ar[0];
  if (nl < 0 || nl > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", nl);
    exit(1);
  }
  const double l = ar[1];

  const double ell = l + 0.5;
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double z = 1./a - 1.;

  const double b1   = gb1(z, nl);
  const double bmag = gbmag(z, nl);

  const double WY = W_y(a);
  const double WGAL = W_gal(a, nl, hoverh0);
  const double WMAG = W_mag(a, fK, nl);

  const double ell_prefactor = l*(l + 1.)/(ell*ell); // prefactor correction (1812.05995 eqs 74-79)

  double res = WY;

  if (include_HOD_GX == 1)
  {
    if (include_RSD_GY == 1)
    {
      log_fatal("RSD not implemented with (HOD = TRUE)");
      exit(1);
    }
    else
    { 
      log_fatal("(HOD = TRUE) not implemented");
      exit(1);
    }
  }
  else
  {
    if (include_RSD_GY == 1)
    {
      log_fatal("RSD not implemented");
      exit(1);
    }
    else
      res *= WGAL*b1 + WMAG*ell_prefactor*bmag;

    const double PK = p_my(k, a);
    res *= PK;
  }
  return res*chidchi.dchida/(fK*fK);
}

double C_gy_tomo_limber_nointerp(const double l, const int ni, const int init)
{
  static double cache[MAX_SIZE_ARRAYS];
  static gsl_integration_glfixed_table* w = NULL;

  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }

  if (w == NULL || fdiff(cache[0], Ntable.random))
  {
    const size_t szint = 200 + 50 * (Ntable.high_def_integration);
    if (w != NULL)  gsl_integration_glfixed_table_free(w);
    w = malloc_gslint_glfixed(szint);
    cache[0] = Ntable.random;
  }

  double ar[2] = {(double) ni, l};
  const double amin = amin_lens(ni);
  const double amax = 0.99999;

  double res = 0.0;
  if (init == 1)
    res = int_for_C_gy_tomo_limber(amin, (void*) ar);
  else
  {
    gsl_function F;
    F.params = (void*) ar;
    F.function = int_for_C_gy_tomo_limber;
    res =  gsl_integration_glfixed(&F, amin, amax, w);
  }
  return res;
}

double C_gy_tomo_limber(double l, int ni)
{
  static double cache[MAX_SIZE_ARRAYS];
  static double** table = NULL;
  static double lim[3];

  if (table == NULL || fdiff(cache[4], Ntable.random))
  {
    lim[0] = log(fmax(limits.LMIN_tab, 1.0));
    lim[1] = log(limits.LMAX + 1);
    lim[2]   = (lim[1] - lim[0])/((double) Ntable.N_ell - 1.0);

    if (table != NULL) free(table);
    table = (double**) malloc2d(redshift.clustering_nbin, Ntable.N_ell);
  }

  if (fdiff(cache[1], cosmology.random) || 
      fdiff(cache[2], nuisance.random_photoz_clustering) ||
      fdiff(cache[3], redshift.random_clustering) ||
      fdiff(cache[4], Ntable.random) ||
      fdiff(cache[5], nuisance.random_galaxy_bias))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    { // init static variables inside the C_XY_limber_nointerp function
      double ini = C_gy_tomo_limber_nointerp(exp(lim[0]), 0, 1);
    }
    #pragma GCC diagnostic pop
    
    #pragma omp parallel for collapse(2)
    for (int k=0; k<redshift.clustering_nbin; k++)
      for (int i=0; i<Ntable.N_ell; i++)
        table[k][i]= C_gy_tomo_limber_nointerp(exp(lim[0] + i*lim[2]), k, 0);

    cache[1] = cosmology.random;
    cache[2] = nuisance.random_photoz_clustering;
    cache[3] = redshift.random_clustering;
    cache[4] = Ntable.random;
    cache[5] = nuisance.random_galaxy_bias;
  }

  const int q =  ni; 
  if (q < 0 || q > redshift.clustering_nbin - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  } 
  
  const double lnl = log(l);
  if (lnl < lim[0])
    log_warn("l = %e < l_min = %e. Extrapolation adopted", l, exp(lim[0]));
  if (lnl > lim[1])
    log_warn("l = %e > l_max = %e. Extrapolation adopted", l, exp(lim[1]));

  return interpol1d(table[q], Ntable.N_ell, lim[0], lim[1], lim[2], lnl);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_ys_tomo_limber(double a, void* params)
{
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }

  double* ar = (double*) params;
  const int ni = (int) ar[0];
  if (ni < 0 || ni > redshift.shear_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }
  const double l = ar[1];
  
  const double ell = l + 0.5;
  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double fK = f_K(chidchi.chi);
  const double k  = ell/fK;
  
  const double PK = p_my(k, a);

  const double WK1 = W_kappa(a, fK, ni);
  const double WY  = W_y(a);

  const double tmp = (l - 1.0)*l*(l + 1.0)*(l + 2.0); // prefactor correction (1812.05995 eqs 74-79)
  const double ell_prefactor2 = (tmp > 0) ? sqrt(tmp)/(ell*ell) : 0.0;

  const double A_Z1 = IA_A1_Z1(a, growfac_a, ni);
  const double WS1  = W_source(a, ni, hoverh0) * A_Z1;

  const double res = (WK1 - WS1)*WY;

  return res*PK*(chidchi.dchida/(fK*fK))*ell_prefactor2;
}

double C_ys_tomo_limber_nointerp(const double l, const int ni, const int init)
{
  static double cache[MAX_SIZE_ARRAYS];
  static gsl_integration_glfixed_table* w = NULL;

  if (ni < 0 || ni > redshift.shear_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  } 

  if (w == NULL || fdiff(cache[0], Ntable.random))
  {
    const size_t szint = 200 + 50 * (Ntable.high_def_integration);
    if (w != NULL)  gsl_integration_glfixed_table_free(w);
    w = malloc_gslint_glfixed(szint);
    cache[0] = Ntable.random;
  }

  double ar[2] = {(double) ni, l};
  const double amin = amin_source(ni);
  const double amax = 0.99999;

  double res = 0.0;
  if (init == 1)
    res = int_for_C_ys_tomo_limber(amin, (void*) ar);
  else
  {
    gsl_function F;
    F.params = (void*) ar;
    F.function = int_for_C_ys_tomo_limber;
    res =  gsl_integration_glfixed(&F, amin, amax, w);
  }
  return res;
}

double C_ys_tomo_limber(double l, int ni)
{
  static double cache[MAX_SIZE_ARRAYS];
  static double** table = NULL;
  static double lim[3];

  if (table == NULL || fdiff(cache[4], Ntable.random))
  {
    if (table != NULL) free(table);
    table = (double**) malloc2d(redshift.shear_nbin, Ntable.N_ell);

    lim[0] = log(fmax(limits.LMIN_tab, 1.0));
    lim[1] = log(limits.LMAX + 1);
    lim[2] = (lim[1] - lim[0])/((double) Ntable.N_ell - 1.0);
  }

  if (fdiff(cache[0], cosmology.random) ||
      fdiff(cache[1], nuisance.random_photoz_shear) ||
      fdiff(cache[2], nuisance.random_ia) ||
      fdiff(cache[3], redshift.random_shear) ||
      fdiff(cache[4], Ntable.random))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    { // init static variables inside the C_XY_limber_nointerp function
      double ini = C_ys_tomo_limber_nointerp(exp(lim[0]), 0, 1);
    }
    #pragma GCC diagnostic pop
    
    #pragma omp parallel for collapse(2)
    for (int k=0; k<redshift.shear_nbin; k++)
      for (int i=0; i<Ntable.N_ell; i++)
        table[k][i] = C_ys_tomo_limber_nointerp(exp(lim[0] + i*lim[2]), k, 0);
          
    cache[0] = cosmology.random;
    cache[1] = nuisance.random_photoz_shear;
    cache[2] = nuisance.random_ia;
    cache[3] = redshift.random_shear;
    cache[4] = Ntable.random;
  }
  
  if (ni < 0 || ni > redshift.shear_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d (max %d)", ni, 
      redshift.shear_nbin);
    exit(1);
  }

  const int q =  ni; 
  if (q < 0 || q > redshift.shear_nbin - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  } 
  
  const double lnl = log(l);
  if (lnl < lim[0])
    log_warn("l = %e < l_min = %e. Extrapolation adopted", l, exp(lim[0]));
  if (lnl > lim[1])
    log_warn("l = %e > l_max = %e. Extrapolation adopted", l, exp(lim[1]));

  return interpol1d(table[q], Ntable.N_ell, lim[0], lim[1], lim[2], lnl);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_ky_limber(double a, void* params)
{
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }

  double *ar = (double*) params;
  const double l = ar[0];

  const double ell = l + 0.5;
  struct chis chidchi = chi_all(a);
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  
  const double PK = p_my(k, a);
  const double WK = W_k(a, fK);
  const double WY = W_y(a);

  const double ell_prefactor = l*(l + 1.0)/(ell*ell); // prefactor correction (1812.05995 eqs 74-79)

  return (WK*WY*PK*chidchi.dchida/(fK*fK))*ell_prefactor;
}

double C_ky_limber_nointerp(const double l, const int init)
{
  static double cache[MAX_SIZE_ARRAYS];
  static gsl_integration_glfixed_table* w = NULL;
  
  if (w == NULL || fdiff(cache[0], Ntable.random))
  {
    const size_t szint = 200 + 50 * (Ntable.high_def_integration);
    if (w != NULL)  gsl_integration_glfixed_table_free(w);
    w = malloc_gslint_glfixed(szint);
    cache[0] = Ntable.random;
  }

  double ar[1] = {l};
  const double amin = limits.a_min_hm;
  const double amax = 1.0 - 1.e-5;

  double res = 0.0;
  if (init == 1)
    res = int_for_C_ky_limber(amin, (void*) ar);
  else
  {
    gsl_function F;
    F.params = (void*) ar;
    F.function = int_for_C_ky_limber;
    res =  gsl_integration_glfixed(&F, amin, amax, w);
  }
  return res;
}

double C_ky_limber(double l)
{
  static double cache[MAX_SIZE_ARRAYS];
  static double* table = NULL;
  static double lim[3];

  if (table == NULL || fdiff(cache[1], Ntable.random))
  {
    if (table != NULL) free(table);
    table = (double*) malloc1d(Ntable.N_ell);

    lim[0] = log(fmax(limits.LMIN_tab, 1.0));
    lim[1] = log(limits.LMAX + 1);
    lim[2] = (lim[1] - lim[0])/((double) Ntable.N_ell - 1.0);
  }

  if (fdiff(cache[0], cosmology.random) || fdiff(cache[1], Ntable.random))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    { // init static variables inside the C_XY_limber_nointerp function
      double ini = C_ky_limber_nointerp(exp(lim[0]), 1);
    }
    #pragma GCC diagnostic pop
    
    #pragma omp parallel for
    for (int i=0; i<Ntable.N_ell; i++)
      table[i] = log(C_ky_limber_nointerp(exp(lim[0] + i*lim[2]), 0));

    cache[0] = cosmology.random;
    cache[1] = Ntable.random;
  }
  
  const double lnl = log(l);
  if (lnl < lim[0])
    log_warn("l = %e < l_min = %e. Extrapolation adopted", l, exp(lim[0]));
  if (lnl > lim[1])
    log_warn("l = %e > l_max = %e. Extrapolation adopted", l, exp(lim[1]));
  
  return interpol1d(table, Ntable.N_ell, lim[0], lim[1], lim[2], lnl);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_yy_limber(double a, void *params)
{
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;
  const double l = ar[0];

  const double ell = l + 0.5;
  struct chis chidchi = chi_all(a);
  const double fK = f_K(chidchi.chi);
  const double k  = ell/fK;

  const double PK = p_yy(k, a);
  
  const double WY = W_y(a);

  return WY*WY*PK*chidchi.dchida/(fK*fK);
}

double C_yy_limber_nointerp(const double l, const int init)
{
  static double cache[MAX_SIZE_ARRAYS];
  static gsl_integration_glfixed_table* w = NULL;
  
  if (w == NULL || fdiff(cache[0], Ntable.random))
  {
    const size_t szint = 200 + 50 * (Ntable.high_def_integration);
    if (w != NULL)  gsl_integration_glfixed_table_free(w);
    w = malloc_gslint_glfixed(szint);
    cache[0] = Ntable.random;
  }

  double ar[2] = {l};
  const double amin = limits.a_min;
  const double amax = 1.0 - 1.e-5;

  double res = 0.0;
  if (init == 1)
    res = int_for_C_yy_limber(amin, (void*) ar);
  else
  {
    gsl_function F;
    F.params = (void*) ar;
    F.function = int_for_C_yy_limber;
    res = gsl_integration_glfixed(&F, amin, amax, w);
  }
  return res;
}

double C_yy_limber(double l)
{
  static double cache[MAX_SIZE_ARRAYS];
  static double* table = NULL;
  static double lim[3];

  if (table == NULL || fdiff(cache[1], Ntable.random))
  {
    lim[0] = log(fmax(limits.LMIN_tab, 1.0));
    lim[1] = log(limits.LMAX + 1.0);
    lim[2] = (lim[1] - lim[0])/((double) Ntable.N_ell - 1.0);

    if (table != NULL) free(table);
    table = (double*) malloc1d(Ntable.N_ell);
  }

  if (fdiff(cache[0], cosmology.random) || fdiff(cache[1], Ntable.random))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    { // init static variables inside the C_XY_limber_nointerp function
      double ini = C_yy_limber_nointerp(exp(lim[0]), 1);
    }
    #pragma GCC diagnostic pop
    
    #pragma omp parallel for
    for (int i=0; i<Ntable.N_ell; i++)
      table[i] = C_yy_limber_nointerp(exp(lim[0] + i*lim[2]), 0);

    cache[0] = cosmology.random;
    cache[1] = Ntable.random;
  }

  const double lnl = log(l);
  if (lnl < lim[0])
    log_warn("l = %e < l_min = %e. Extrapolation adopted", l, exp(lim[0]));
  if (lnl > lim[1])
    log_warn("l = %e > l_max = %e. Extrapolation adopted", l, exp(lim[1]));
    
  return interpol1d(table, Ntable.N_ell, lim[0], lim[0], lim[2], lnl);
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Non Limber (Angular Power Spectrum)
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

void C_cl_tomo(
    int L, 
    const int ni, 
    const int nj, 
    double *const Cl, 
    double dev, 
    double tol
  )
{
  if (ni < -1 || ni > redshift.clustering_nbin - 1 || 
      nj < -1 || nj > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number (ni, nj) = [%d,%d]", ni, nj);
    exit(1);
  }
  if (ni != nj)
  {
    log_fatal("Cocoa disabled cross-spectrum w_gg");
    exit(1);
  }
    
  const double real_coverH0 = cosmology.coverH0/cosmology.h0;
  const double chi_min = chi(1./(1.0 + 0.002))*real_coverH0; // DIMENSIONELESS
  const double chi_max = chi(1./(1.0 + 4.0))*real_coverH0;   // DIMENSIONELESS
  const double dlnchi = log(chi_max/chi_min) / ((double) Ntable.NL_Nchi - 1.0);
  const double dlnk = dlnchi;

  int* ell_ar = (int*) malloc(sizeof(int)*Ntable.NL_Nell_block);
  double*** Fk1 = (double***) malloc3d(3, Ntable.NL_Nell_block, Ntable.NL_Nchi); // (Fk1, Fk1_Mag, k1)    
  double** f1_chi = (double**) malloc2d(4, Ntable.NL_Nchi); // (f1, f1_RSD, f1_MAG, chi)

  #pragma omp parallel for
  for (int i=0; i<Ntable.NL_Nchi; i++)
  {
    f1_chi[3][i] = chi_min * exp(dlnchi * i); 
    const double a = a_chi(f1_chi[3][i]/real_coverH0);
    const double z = 1. / a - 1.;
    
    if (z < redshift.clustering_zdist_zmin[ni] || 
        z > redshift.clustering_zdist_zmax[ni])
    { 
      f1_chi[0][i] = 0.;
      f1_chi[1][i] = 0.;
      f1_chi[2][i] = 0.;
    }
    else
    {
      const double pf = pf_photoz(z,ni);
      const double hoverh0_a = hoverh0(a);
      const double fK = f_K(f1_chi[3][i]/real_coverH0); 
      const double WM = W_mag(a, fK, ni);
      struct growths growfac_a = growfac_all(a);
      const double D = growfac_a.D;
      const double f = growfac_a.f;

      f1_chi[0][i] = gb1(z, ni)*f1_chi[3][i]*pf*D*hoverh0_a/real_coverH0;
      
      f1_chi[1][i] = -f1_chi[3][i]*pf*D*f*hoverh0_a/real_coverH0;

      f1_chi[2][i] = (WM/fK/(real_coverH0*real_coverH0)) * D; // [Mpc^-2] 
    }
  }

  config cfg;
  cfg.nu = 1.;
  cfg.c_window_width = 0.25;
  cfg.derivative = 0;
  cfg.N_pad = 200;
  cfg.N_extrap_low = 0;
  cfg.N_extrap_high = 0;

  config cfg_RSD;
  cfg_RSD.nu = 1.01;
  cfg_RSD.c_window_width = 0.25;
  cfg_RSD.derivative = 2;
  cfg_RSD.N_pad = 500;
  cfg_RSD.N_extrap_low = 0;
  cfg_RSD.N_extrap_high = 0;

  config cfg_Mag;
  cfg_Mag.nu = 1.;
  cfg_Mag.c_window_width = 0.25;
  cfg_Mag.derivative = 0;
  cfg_Mag.N_pad = 500;
  cfg_Mag.N_extrap_low = 0;
  cfg_Mag.N_extrap_high = 0;

  int i_block = 0;

  while ((fabs(dev) > tol) && (L < limits.LMAX_NOLIMBER))
  {
    for (int i=0; i<Ntable.NL_Nell_block; i++)
      ell_ar[i] = i + i_block * Ntable.NL_Nell_block; 

    i_block++;
  
    if (L >= limits.LMAX_NOLIMBER - Ntable.NL_Nell_block) break; //Xiao: break before memory leak in next iteration
  
    L = i_block * Ntable.NL_Nell_block - 1;

    cfftlog_ells(f1_chi[3], f1_chi[0], Ntable.NL_Nchi, &cfg, ell_ar, 
      Ntable.NL_Nell_block, Fk1[2], Fk1[0]);
    
    cfftlog_ells_increment(f1_chi[3], f1_chi[1], Ntable.NL_Nchi, &cfg_RSD, 
      ell_ar, Ntable.NL_Nell_block, Fk1[2], Fk1[0]);   
    
    cfftlog_ells(f1_chi[3], f1_chi[2], Ntable.NL_Nchi, &cfg_Mag, ell_ar, 
      Ntable.NL_Nell_block, Fk1[2], Fk1[1]);    
    
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
    { // init static variables inside the C_XY_limber_nointerp function
      double ini = p_lin(Fk1[2][0][0]*real_coverH0, 1.0);
      ini = C_gg_tomo_limber_nointerp((double) 100, 0, 0, 1);
    }
    #pragma GCC diagnostic pop
    #pragma GCC diagnostic pop

    #pragma omp parallel for
    for (int i=0; i<Ntable.NL_Nell_block; i++)
    {
      const double ell_prefactor = ell_ar[i] * (ell_ar[i] + 1.);

      double cl_temp = 0.;
      for (int j=0; j<Ntable.NL_Nchi; j++)
      {
        Fk1[0][i][j] += 
          gbmag(0.0, ni)*ell_prefactor*Fk1[1][i][j]/(Fk1[2][i][j]*Fk1[2][i][j]);
        
        const double k1cH0 = Fk1[2][i][j] * real_coverH0;
        
        cl_temp += Fk1[0][i][j]*Fk1[0][i][j]*(k1cH0*k1cH0*k1cH0)*p_lin(k1cH0,1);
      }
      
      Cl[ell_ar[i]] = cl_temp * dlnk * 2. / M_PI + 
         C_gg_tomo_limber_linpsopt_nointerp((double) ell_ar[i], ni, nj, 0, 0)
        -C_gg_tomo_limber_linpsopt_nointerp((double) ell_ar[i], ni, nj, 1, 0);
    }
    dev = Cl[L]/C_gg_tomo_limber_nointerp(L, ni, nj, 0) - 1;
  }
  L++;

  Cl[limits.LMAX_NOLIMBER] = C_gg_tomo_limber(limits.LMAX_NOLIMBER, ni, nj);
  
  #pragma omp parallel for
  for (int l=L; l<limits.LMAX_NOLIMBER; l++)
    Cl[l] = (l > limits.LMIN_tab) ? C_gg_tomo_limber(l, ni, nj) :
                                    C_gg_tomo_limber_nointerp(l, ni, nj, 0);
  free(Fk1);
  free(f1_chi);
  free(ell_ar);
}
