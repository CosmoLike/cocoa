#include <assert.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../cfftlog/cfftlog.h"
#include <fftw3.h>

#include "bias.h"
#include "basics.h"
#include "cosmo3D.h"
#include "cosmo2D.h"
#include "halo.h"
#include "IA.h"
#include "pt_cfastpt.h"
#include "radial_weights.h"
#include "recompute.h"
#include "redshift_spline.h"
#include "structs.h"

#include "log.c/src/log.h"

static int GSL_WORKSPACE_SIZE = 250;
static int use_linear_ps_limber = 0; // 0 or 1 

static int include_HOD_GX = 0; // 0 or 1

static int include_RSD_GS = 0; // 0 or 1 
static int include_RSD_GG = 1; // 0 or 1 
static int include_RSD_GK = 0; // 0 or 1
static int include_RSD_GY = 0; // 0 or 1

double beam_cmb(const double l)
{
  const double sigma = cmb.fwhm/sqrt(16.0*log(2.0));
  double norm = 1.0;
  if(l < like.lmin_kappacmb || l > like.lmax_kappacmb)
  {
    norm = 0.0;
  }
  return exp(-l*(l+1.0)*sigma*sigma)*norm;
}

double w_pixel(const double ell)
{
  static int lbins;
  static double *cl_pixel =0;
  FILE* ein;
  int l = (int) ceil(ell);

  if (cl_pixel == 0)
  {
    ein = fopen(cmb.pathHealpixWinFunc, "r");
    if (ein == NULL)
    {
      log_info("Can not open file %s, ignore healpix window function\n",
        cmb.pathHealpixWinFunc);
      lbins = 0;
      cl_pixel = malloc(sizeof(double));
      cl_pixel[0] = 1.0;
    }
    else
    {
      fclose(ein);
      lbins = line_count(cmb.pathHealpixWinFunc);
      cl_pixel = malloc(sizeof(double)*lbins);
      ein = fopen(cmb.pathHealpixWinFunc, "r");
      for (int i = 0; i < lbins; i++)
      {
        int tmp;
        double tmp2;
        fscanf(ein, "%d %le\n", &tmp, &tmp2);
        cl_pixel[i] = tmp2;
      }
      fclose(ein);
    }
  }
  if (lbins>0){return (l < lbins) ? cl_pixel[l] : 0.0;}
  else{return 1.0;}
}

double C_gk_tomo_limber_nointerp_wrapper(double l, int ni, int use_linear_ps,
const int init_static_vars_only)
{
  return C_gk_tomo_limber_nointerp(l, ni, use_linear_ps, init_static_vars_only)*beam_cmb(l)*w_pixel(l);
}

double C_gk_tomo_limber_wrapper(double l, int ni)
{
  return C_gk_tomo_limber(l, ni)*beam_cmb(l)*w_pixel(l);
}

double C_ks_tomo_limber_nointerp_wrapper(double l, int ni, int use_linear_ps, 
const int init_static_vars_only)
{
  return C_ks_tomo_limber_nointerp(l, ni, use_linear_ps, init_static_vars_only)*beam_cmb(l)*w_pixel(l);
}

double C_ks_tomo_limber_wrapper(double l, int ni)
{
  return C_ks_tomo_limber(l, ni)*beam_cmb(l)*w_pixel(l);
}

static int has_b2_galaxies()
{
  int res = 0;
  for (int i=0; i<tomo.clustering_Nbin; i++) 
  {
    if (gbias.b2[i])
    {
      res = 1;
    }
  }
  return res;
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Correlation Functions (real Space) - Full Sky - bin average
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double xi_pm_tomo(const int pm, const int nt, const int ni, const int nj, const int limber)
{
  if (like.Ntheta == 0)
  {
    log_fatal("like.Ntheta not initialized");
    exit(1);
  }

  static double** Glplus = 0;
  static double** Glminus = 0;
  static double* xi_vec_plus = 0;
  static double* xi_vec_minus = 0;
  static cosmopara C;
  static nuisancepara N;
  
  const int nell = limits.LMAX;
  const int ntheta = like.Ntheta;
  const int NSIZE = tomo.shear_Npowerspectra;

  if (Glplus == 0)
  {
    Glplus = (double**) malloc(sizeof(double*)*ntheta);
    for (int i=0; i<ntheta; i++) 
    {
      Glplus[i] = (double*) calloc(nell, sizeof(double));
    }
    Glminus = (double**) malloc(sizeof(double*)*ntheta);
    for (int i=0; i<ntheta; i++) 
    {
      Glminus[i] = (double*) calloc(nell, sizeof(double));
    }
    xi_vec_plus = (double*) calloc(NSIZE*ntheta, sizeof(double));
    xi_vec_minus = (double*) calloc(NSIZE*ntheta, sizeof(double));

    double xmin[ntheta];
    double xmax[ntheta];
    for (int i=0; i<ntheta; i++)
    { // Cocoa: dont thread (init of static variables inside set_bin_average)
      bin_avg r = set_bin_average(i,0);
      xmin[i] = r.xmin;
      xmax[i] = r.xmax;
    }

    double** Pmin = (double**) malloc(sizeof(double)*ntheta);
    double** Pmax = (double**) malloc(sizeof(double)*ntheta);
    double** dPmin = (double**) malloc(sizeof(double)*ntheta);
    double** dPmax = (double**) malloc(sizeof(double)*ntheta);
    for (int i=0; i<ntheta; i ++)
    {
      Pmin[i] = (double*) calloc((nell + 1), sizeof(double));
      Pmax[i] = (double*) calloc((nell + 1), sizeof(double));
      dPmin[i] = (double*) calloc((nell + 1), sizeof(double));
      dPmax[i] = (double*) calloc((nell + 1), sizeof(double));
    }
    #pragma omp parallel for collapse(2)
    for (int i=0; i<ntheta; i ++)
    {
      for (int l=0; l<nell; l++)
      {
        bin_avg r = set_bin_average(i, l);
        Pmin[i][l] = r.Pmin;
        Pmax[i][l] = r.Pmax;
        dPmin[i][l] = r.dPmin;
        dPmax[i][l] = r.dPmax;
      }
    }
    #pragma omp parallel for collapse(2)
    for (int i=0; i<ntheta; i ++)
    {
      for (int l=1; l<nell; l++)
      {
        Glplus[i][l] = (2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(
          -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[i][l-1]-Pmax[i][l-1])
          -l*(l-1.)*(2.-l)/2 * (xmin[i]*Pmin[i][l]-xmax[i]*Pmax[i][l])
          +l*(l-1.)/(2.*l+1) * (Pmin[i][l+1]-Pmax[i][l+1])
          +(4-l)*(dPmin[i][l]-dPmax[i][l])
          +(l+2)*(xmin[i]*dPmin[i][l-1] - xmax[i]*dPmax[i][l-1] - Pmin[i][l-1] + Pmax[i][l-1])
          +2*(l-1)*(xmin[i]*dPmin[i][l] - xmax[i]*dPmax[i][l] - Pmin[i][l] + Pmax[i][l])
          -2*(l+2)*(dPmin[i][l-1]-dPmax[i][l-1])
        )/(xmin[i]-xmax[i]);

        Glminus[i][l] = (2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(
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
    for (int i=0; i<ntheta; i++)
    {
      free(Pmin[i]);
      free(Pmax[i]);
      free(dPmin[i]);
      free(dPmax[i]);
    }
    free(Pmin);
    free(Pmax);
    free(dPmin);
    free(dPmax);
  }
  if (recompute_shear(C, N))
  {
    if (limber == 1)
    {
      if (like.IA == 5 || like.IA == 6)
      { // NEW TATT MODELING
        double** Cl_EE = (double**) malloc(NSIZE*sizeof(double*));
        double** Cl_BB = (double**) malloc(NSIZE*sizeof(double*));
        for (int i = 0; i<NSIZE; i++)
        {
          Cl_EE[i] = calloc(nell, sizeof(double));
          Cl_BB[i] = calloc(nell, sizeof(double));
        }
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-variable"
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
        { // init the functions C_ss_tomo_TATT_EE/BB_limber
          // only compute BB if the TATT parameters allow for B-mode terms
          const int nz = 0;
          const int Z1NZ = Z1(nz);
          const int Z2NZ = Z2(nz);

          const int BM = (nuisance.b_ta_z[0] || nuisance.b_ta_z[Z1NZ] ||
                          nuisance.b_ta_z[Z2NZ] || nuisance.A2_ia ||
                          nuisance.A2_z[Z1NZ] || nuisance.A2_z[Z2NZ]) ? 1 : 0;
          double x = C_ss_tomo_TATT_EE_limber(limits.LMIN_tab + 1, Z1(0), Z2(0)); 
          x = (BM == 1) ? C_ss_tomo_TATT_BB_limber(limits.LMIN_tab + 1, Z1(0), Z2(0)) : 0.0;
        }
        #pragma GCC diagnostic pop
        #pragma GCC diagnostic pop
        #pragma omp parallel for collapse(2)
        for (int nz=0; nz<NSIZE; nz++) 
        {
          for (int l=2; l<nell; l++)
          {
            const int Z1NZ = Z1(nz);
            const int Z2NZ = Z2(nz);
          
            // only compute BB if the TATT parameters allow for B-mode terms
            const int BM = (nuisance.b_ta_z[0] || nuisance.b_ta_z[Z1NZ] ||
                          nuisance.b_ta_z[Z2NZ] || nuisance.A2_ia ||
                          nuisance.A2_z[Z1NZ] || nuisance.A2_z[Z2NZ]) ? 1 : 0;

            Cl_EE[nz][l] = (l > limits.LMIN_tab) ?
              C_ss_tomo_TATT_EE_limber((double) l, Z1NZ, Z2NZ) :
              C_ss_tomo_TATT_EE_limber_nointerp((double) l, Z1NZ, Z2NZ, 0);
        
            Cl_BB[nz][l] = (BM == 1) ? (l > limits.LMIN_tab) ? 
              C_ss_tomo_TATT_BB_limber((double) l, Z1NZ, Z2NZ) :
              C_ss_tomo_TATT_BB_limber_nointerp((double) l, Z1NZ, Z2NZ, 0) : 0.0;
          }
        }
        #pragma omp parallel for collapse(2)
        for (int nz=0; nz<NSIZE; nz++)
        {
          for (int i=0; i<ntheta; i++)
          {
            const int q = nz*ntheta + i;
            xi_vec_plus[q] = 0;
            xi_vec_minus[q] = 0;
            for (int l=2; l<nell; l++)
            {
              xi_vec_plus[q] += Glplus[i][l] * (Cl_EE[nz][l] + Cl_BB[nz][l]);
              xi_vec_minus[q] += Glminus[i][l] * (Cl_EE[nz][l] - Cl_BB[nz][l]);
            }
          }
        }
        for (int nz=0; nz<NSIZE; nz++)
        {
          free(Cl_EE[nz]);
          free(Cl_BB[nz]);
        }
        free(Cl_EE);
        free(Cl_BB);
      }
      else
      {
        double** Cl = malloc(NSIZE*sizeof(double*));
        for (int i=0; i<NSIZE; i++)
        {
          Cl[i] = calloc(nell, sizeof(double));
        }
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-variable"
        {
          double init_static_vars_only = C_ss_tomo_limber(limits.LMIN_tab, Z1(0), Z2(0)); 
        }
        #pragma GCC diagnostic pop
        #pragma omp parallel for collapse(2)
        for (int nz=0; nz<NSIZE; nz++)
        {
          for (int l=2; l<nell; l++)
          {
            const int Z1NZ = Z1(nz);
            const int Z2NZ = Z2(nz);
            Cl[nz][l] = (l > limits.LMIN_tab) ?
              C_ss_tomo_limber((double) l, Z1NZ, Z2NZ) :
              C_ss_tomo_limber_nointerp((double) l, Z1NZ, Z2NZ, use_linear_ps_limber, 0);
          }
        }
        #pragma omp parallel for collapse(2)
        for (int nz=0; nz<NSIZE; nz++)
        {
          for (int i=0; i<ntheta; i++)
          {
            const int q = nz*ntheta + i;
            xi_vec_plus[q] = 0;
            xi_vec_minus[q] = 0;
            for (int l=2; l<nell; l++)
            {
              xi_vec_plus[q] += Glplus[i][l]*Cl[nz][l];
              xi_vec_minus[q] += Glminus[i][l]*Cl[nz][l];
            }
          }
        }
        for (int nz=0; nz<NSIZE; nz++)
        {
          free(Cl[nz]);
        }
        free(Cl);
      }
    }
    else
    {
      log_fatal("NonLimber not implemented");
      exit(1);
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  if (nt < 0 || nt > like.Ntheta - 1)
  {
    log_fatal("error in selecting bin number nt = %d (max %d)", nt, like.Ntheta);
    exit(1); 
  }
  if (ni < 0 || ni > tomo.shear_Nbin - 1 || nj < 0 || nj > tomo.shear_Nbin - 1)
  {
    log_fatal("error in selecting bin number (ni, nj) = [%d,%d]", ni, nj);
    exit(1);
  }

  const int ntomo = N_shear(ni, nj);
  const int q = ntomo*ntheta + nt;
  if (q < 0 || q > NSIZE*ntheta - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  }
  return (pm > 0) ? xi_vec_plus[q] : xi_vec_minus[q];
}

// ---------------------------------------------------------------------------

double w_gammat_tomo(const int nt, const int ni, const int nj, const int limber)
{
  if (like.Ntheta == 0)
  {
    log_fatal("like.Ntheta not initialized");
    exit(1);
  }

  static double** Pl = 0;
  static double* w_vec = 0;
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
    
  const int nell = limits.LMAX;
  const int ntheta = like.Ntheta;
  const int NSIZE = tomo.ggl_Npowerspectra;

  if (Pl == 0)
  {    
    Pl = (double**) malloc(sizeof(double*)*ntheta);
    for (int i=0; i<ntheta; i++) 
    {
      Pl[i] = (double*) calloc(nell, sizeof(double));
    }
    w_vec = (double*) calloc(NSIZE*ntheta, sizeof(double)); 
    double xmin[ntheta];
    double xmax[ntheta];
    // Cocoa: dont thread (init of static variables inside set_bin_average)
    for (int i=0; i<ntheta; i++)
    {
      bin_avg r = set_bin_average(i,0);
      xmin[i] = r.xmin;
      xmax[i] = r.xmax;
    }

    double** Pmin = (double**) malloc(sizeof(double)*ntheta);
    double** Pmax = (double**) malloc(sizeof(double)*ntheta);
    for (int i=0; i<ntheta; i++)
    {
      Pmin[i] = (double*) calloc((nell + 1), sizeof(double));
      Pmax[i] = (double*) calloc((nell + 1), sizeof(double));
    }
    #pragma omp parallel for collapse(2)
    for (int i=0; i<ntheta; i ++)
    {
      for (int l=0; l<nell; l++)
      {
        bin_avg r = set_bin_average(i, l);
        Pmin[i][l] = r.Pmin;
        Pmax[i][l] = r.Pmax;
      }
    }
    #pragma omp parallel for collapse(2)
    for (int i=0; i<ntheta; i++)
    {
      for (int l=1; l<nell; l++)
      {
        Pl[i][l] = (2.*l+1)/(4.*M_PI*l*(l+1)*(xmin[i]-xmax[i]))
          *((l+2./(2*l+1.))*(Pmin[i][l-1]-Pmax[i][l-1])
          +(2-l)*(xmin[i]*Pmin[i][l]-xmax[i]*Pmax[i][l])
          -2./(2*l+1.)*(Pmin[i][l+1]-Pmax[i][l+1]));
      }
    }
    for (int i=0; i<ntheta; i++)
    {
      free(Pmin[i]);
      free(Pmax[i]);
    }
    free(Pmin);
    free(Pmax);
  }
  if (recompute_gs(C, G, N))
  {
    double** Cl = malloc(NSIZE*sizeof(double*));
    for (int nz = 0; nz<NSIZE; nz++)
    {
      Cl[nz] = calloc(nell, sizeof(double));
    }
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init_static_vars_only = C_gs_tomo_limber(limits.LMIN_tab, ZL(0), ZS(0));
    }
    #pragma GCC diagnostic pop
    if (limber == 1)
    {
      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<NSIZE; nz++)
      {
        for (int l=2; l<nell; l++)
        {
          const int ZLNZ = ZL(nz);
          const int ZSNZ = ZS(nz);
          Cl[nz][l] = (l > limits.LMIN_tab) ?
            C_gs_tomo_limber((double) l, ZLNZ, ZSNZ) :
            C_gs_tomo_limber_nointerp((double) l, ZLNZ, ZSNZ, use_linear_ps_limber, 0);
        }
      }
    }
    else
    {
      for (int nz=0; nz<NSIZE; nz++) // NONLIMBER PART
      { 
        const int L = 1;
        const double tolerance = 0.0075;    // required fractional accuracy in C(l)
        const double dev = 10. * tolerance; // will be diff exact vs Limber init to large
                                            // value in order to start while loop
        const int ZLNZ = ZL(nz);
        const int ZSNZ = ZS(nz);
        C_gl_tomo(L, ZLNZ, ZSNZ, Cl[nz], dev, tolerance);
      }
      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<NSIZE; nz++) // LIMBER PART
      {
        for (int l=limits.LMAX_NOLIMBER+1; l<nell; l++)
        {
          const int ZLNZ = ZL(nz);
          const int ZSNZ = ZS(nz);
          Cl[nz][l] = C_gs_tomo_limber((double) l, ZLNZ, ZSNZ);
        }
      }
    }
    #pragma omp parallel for collapse(2)
    for (int nz=0; nz<NSIZE; nz++)
    {
      for (int i=0; i<ntheta; i++)
      {
        const int q = nz*ntheta+i;
        w_vec[q] = 0;
        for (int l=1; l<nell; l++)
        {
          w_vec[q] += Pl[i][l]*Cl[nz][l];
        }
      }
    }
    for (int nz=0; nz<NSIZE; nz++)
    {
      free(Cl[nz]);
    }
    free(Cl);
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  if (nt < 0 || nt > like.Ntheta - 1)
  {
    log_fatal("error in selecting bin number nt = %d (max %d)", nt, like.Ntheta);
    exit(1); 
  }
  if (ni < -1 || ni > tomo.clustering_Nbin - 1 || nj < -1 || nj > tomo.shear_Nbin - 1)
  {
    log_fatal("error in selecting bin number (ni, nj) = [%d,%d]", ni, nj);
    exit(1);
  }
  const int q = N_ggl(ni, nj)*ntheta + nt;
  if (q < 0 || q > NSIZE*ntheta - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  }
  return w_vec[q];
}

// ---------------------------------------------------------------------------

double w_gg_tomo(const int nt, const int ni, const int nj, const int limber)
{
  if (like.Ntheta == 0)
  {
    log_fatal("like.Ntheta not initialized");
    exit(1);
  }

  static double** Pl = 0;
  static double* w_vec = 0;
  static cosmopara C;
  static nuisancepara N;
  static galpara G;

  const int nell = limits.LMAX;
  const int ntheta = like.Ntheta;
  const int NSIZE = tomo.clustering_Npowerspectra;

  if (Pl == 0)
  {
    Pl = (double**) malloc(sizeof(double*)*ntheta);
    for (int i=0; i<ntheta; i++) 
    {
      Pl[i] = (double*) calloc(nell, sizeof(double));
    }
    w_vec = (double*) calloc(NSIZE*ntheta, sizeof(double));
    double xmin[ntheta];
    double xmax[ntheta];
    // Cocoa: dont thread (init of static variables inside set_bin_average)
    for (int i=0; i<ntheta; i ++)
    {
      bin_avg r = set_bin_average(i,0);
      xmin[i] = r.xmin;
      xmax[i] = r.xmax;
    }
    double** Pmin = (double**) malloc(sizeof(double)*ntheta);
    double** Pmax = (double**) malloc(sizeof(double)*ntheta);
    for (int i=0; i<ntheta; i ++)
    {
      Pmin[i] = (double*) calloc((nell + 1), sizeof(double));
      Pmax[i] = (double*) calloc((nell + 1), sizeof(double));
    }
    #pragma omp parallel for
    for (int i=0; i<ntheta; i ++)
    {
      for (int l=0; l<nell; l++)
      {
        bin_avg r = set_bin_average(i,l);
        Pmin[i][l] = r.Pmin;
        Pmax[i][l] = r.Pmax;
      }
      Pl[i][0] = 1.0;
      const double tmp = (1.0/(xmin[i] - xmax[i]))*(1. / (4.0 * M_PI));
      for (int l=1; l<nell; l++)
      {
        Pl[i][l] = tmp*(Pmin[i][l + 1] - Pmax[i][l + 1] - Pmin[i][l - 1] + Pmax[i][l - 1]);
      }
    }
    for (int i=0; i<ntheta; i ++)
    {
      free(Pmin[i]);
      free(Pmax[i]);
    }
    free(Pmin);
    free(Pmax);
  }

  if (recompute_gg(C, G, N))
  {    
    double** Cl = malloc(NSIZE*sizeof(double*));
    for (int nz=0; nz<NSIZE; nz++)
    {
      Cl[nz] = calloc(nell, sizeof(double));
    }
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init_static_vars_only = C_gg_tomo_limber(limits.LMIN_tab, 0, 0);
    }
    #pragma GCC diagnostic pop
    if (limber == 1)
    {
      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<NSIZE; nz++)
      {
        for (int l=1; l<nell; l++)
        {
          const int q = nz;
          const int Z1 = nz; // cross redshift bin not supported so not using ZCL1(k)
          const int Z2 = nz; // cross redshift bin not supported so not using ZCL2(k)
          Cl[q][l] = (l > limits.LMIN_tab) ?
            C_gg_tomo_limber((double) l, Z1, Z2) :
            C_gg_tomo_limber_nointerp((double) l, Z1, Z2, use_linear_ps_limber, 0);
        }
      }
    }
    else
    {
      for (int nz=0; nz<NSIZE; nz++) // NONLIMBER PART
      { 
        const int L = 1;
        const double tolerance = 0.01;    // required fractional accuracy in C(l)
        const double dev = 10. * tolerance; // will be diff  exact vs Limber init to
                                            // large value in order to start while loop
        const int Z1 = nz; // cross redshift bin not supported so not using ZCL1(k)
        const int Z2 = nz; // cross redshift bin not supported so not using ZCL2(k)
        C_cl_tomo(L, Z1, Z2, Cl[nz], dev, tolerance);
      }
      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<NSIZE; nz++) // LIMBER PART
      {
        for (int l=limits.LMAX_NOLIMBER+1; l<nell; l++)
        {
          const int q = nz;
          const int Z1 = nz; // cross redshift bin not supported so not using ZCL1(k)
          const int Z2 = nz; // cross redshift bin not supported so not using ZCL2(k)
          Cl[q][l] = C_gg_tomo_limber((double) l, Z1, Z2);
        }
      }
    }
    #pragma omp parallel for collapse(2)
    for (int nz=0; nz<NSIZE; nz++)
    {
      for (int i=0; i<ntheta; i++)
      {
        const int q = nz*ntheta + i;
        w_vec[q] = 0;
        for (int l=1; l<nell; l++)
        {
          w_vec[q] += Pl[i][l]*Cl[nz][l];
        }
      }
    }
    for (int nz=0; nz<NSIZE; nz++)
    {
      free(Cl[nz]);
    }
    free(Cl);
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  if (nt < 0 || nt > like.Ntheta - 1)
  {
    log_fatal("error in selecting bin number nt = %d (max %d)", nt, like.Ntheta);
    exit(1); 
  }
  if (ni < -1 || ni > tomo.clustering_Nbin - 1 || nj < -1 || nj > tomo.clustering_Nbin - 1)
  {
    log_fatal("error in selecting bin number (ni, nj) = [%d,%d]", ni, nj);
    exit(1);
  }
  if (ni != nj)
  {
    log_fatal("ni != nj tomography not supported");
    exit(1);
  }
  const int q = ni * ntheta + nt;
  if (q  < 0 || q > NSIZE*ntheta - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  }  
  return w_vec[q];
}

// ---------------------------------------------------------------------------

double w_gk_tomo(const int nt, const int ni, const int limber)
{
  if (like.Ntheta == 0)
  {
    log_fatal("like.Ntheta not initialized");
    exit(1);
  }

  static double** Pl =0;
  static double* w_vec = 0;
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  
  const int nell = limits.LMAX;
  const int ntheta = like.Ntheta;
  const int NSIZE = tomo.clustering_Nbin;
  
  if (Pl == 0)
  {
    Pl = (double**) malloc(sizeof(double*)*ntheta);
    for (int i=0; i<ntheta; i++) 
    {
      Pl[i] = (double*) calloc(nell, sizeof(double));
    }
    w_vec = (double*) calloc(NSIZE*ntheta, sizeof(double));
    double xmin[ntheta];
    double xmax[ntheta];
    // Cocoa: dont thread (init of static variables inside set_bin_average)
    for (int i=0; i<ntheta; i++)
    {
      bin_avg r = set_bin_average(i,0);
      xmin[i] = r.xmin;
      xmax[i] = r.xmax;
    }
    double** Pmin = (double**) malloc(sizeof(double)*ntheta);
    double** Pmax = (double**) malloc(sizeof(double)*ntheta);
    for (int i=0; i<ntheta; i ++)
    {
      Pmin[i] = (double*) calloc((nell + 1), sizeof(double));
      Pmax[i] = (double*) calloc((nell + 1), sizeof(double));
    }
    #pragma omp parallel for
    for (int i=0; i<ntheta; i++)
    {
      for (int l=0; l<nell; l++)
      {
        bin_avg r = set_bin_average(i,l);
        Pmin[i][l] = r.Pmin;
        Pmax[i][l] = r.Pmax;
      }
      const double tmp = (1.0/(xmin[i] - xmax[i]))*(1.0 / (4.0 * M_PI));
      for (int l=1; l<nell; l++)
      {
        Pl[i][l] = tmp*(Pmin[i][l + 1] - Pmax[i][l + 1] - Pmin[i][l - 1] + Pmax[i][l - 1]);
      }
    }
    for (int i=0; i<ntheta; i++)
    {
      free(Pmin[i]);
      free(Pmax[i]);
    }
    free(Pmin);
    free(Pmax);
  }
  if (recompute_gk(C, G, N))
  {
    if (limber == 1)
    {
      double** Cl = malloc(NSIZE*sizeof(double*));
      for (int nz = 0; nz<NSIZE; nz++)
      {
        Cl[nz] = calloc(nell, sizeof(double));
      }
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wunused-variable"
      {
        double init_static_vars_only = C_gk_tomo_limber_wrapper(limits.LMIN_tab + 1, 0);
      }
      #pragma GCC diagnostic pop
      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<NSIZE; nz++)
      {
        for (int l=1; l<nell; l++)
        {
          Cl[nz][l] = (l > limits.LMIN_tab) ?
            C_gk_tomo_limber_wrapper((double) l, nz) :
            C_gk_tomo_limber_nointerp_wrapper((double) l, nz, use_linear_ps_limber, 0);
        }
      }
      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<NSIZE; nz++)
      {
        for (int i=0; i<ntheta; i++)
        {
          w_vec[nz*ntheta+i] = 0;
          for (int l=1; l<nell; l++)
          {
            w_vec[nz*ntheta+i] += Pl[i][l]*Cl[nz][l];
          }
        }
      }
      for (int nz=0; nz<NSIZE; nz++)
      {
        free(Cl[nz]);
      }
      free(Cl);
    }
    else
    {
      log_fatal("NonLimber not implemented");
      exit(1);
    }
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  if (ni < -1 || ni > tomo.clustering_Nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d (max %d)", ni, tomo.clustering_Nbin);
    exit(1); 
  } 
  if (nt < 0 || nt > like.Ntheta - 1)
  {
    log_fatal("error in selecting bin number nt = %d (max %d)", nt, like.Ntheta);
    exit(1); 
  } 
  const int q = ni * ntheta + nt;
  if (q < 0 || q > NSIZE*ntheta - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  }  
  return w_vec[q];
}

// ---------------------------------------------------------------------------

double w_ks_tomo(const int nt, const int ni, const int limber)
{
  if (like.Ntheta ==0)
  {
    log_fatal("like.Ntheta not initialized");
    exit(1);
  }

  static double** Pl =0;
  static double* w_vec =0;
  static cosmopara C;
  static nuisancepara N;
  
  const int nell = limits.LMAX;
  const int ntheta = like.Ntheta;
  const int NSIZE = tomo.shear_Nbin;

  if (Pl == 0)
  {
    Pl = (double**) calloc(ntheta, sizeof(double*));
    for (int i=0; i<ntheta; i++) 
    {
      Pl[i] = (double*) calloc(nell, sizeof(double));
    }
    w_vec = (double*) calloc(NSIZE*ntheta, sizeof(double));
    double xmin[ntheta];
    double xmax[ntheta];
    // Cocoa: dont thread (init of static variables inside set_bin_average)
    for (int i=0; i<ntheta; i++)
    {
      bin_avg r = set_bin_average(i,0);
      xmin[i] = r.xmin;
      xmax[i] = r.xmax;
    }
    #pragma omp parallel for
    for (int i=0; i<ntheta; i++)
    {
      double* Pmin = (double*) calloc((nell + 1), sizeof(double));
      double* Pmax = (double*) calloc((nell + 1), sizeof(double));
      for (int l=0; l<nell; l++)
      {
        bin_avg r = set_bin_average(i,l);
        Pmin[l] = r.Pmin;
        Pmax[l] = r.Pmax;
      }
      for (int l=1; l<nell; l++)
      {
        Pl[i][l] = (2.*l+1)/(4.*M_PI*l*(l+1)*(xmin[i]-xmax[i]))
          *((l+2./(2*l+1.))*(Pmin[l-1]-Pmax[l-1])
          +(2-l)*(xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
          -2./(2*l+1.)*(Pmin[l+1]-Pmax[l+1]));
      }
      free(Pmin);
      free(Pmax);
    }
  }
  if (recompute_ks(C, N))
  {
    if (limber == 1)
    {
      double** Cl = malloc(tomo.shear_Nbin*sizeof(double*));
      for (int nz=0; nz<tomo.shear_Nbin; nz++)
      {
        Cl[nz] = calloc(nell, sizeof(double));
      }
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wunused-variable"
      {
        double init_static_vars_only = C_ks_tomo_limber_wrapper(limits.LMIN_tab + 1, 0.0);
      }
      #pragma GCC diagnostic pop
      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<tomo.shear_Nbin; nz++)
      {
        for (int l=1; l<nell; l++)
        {
          Cl[nz][l] = (l > limits.LMIN_tab) ?
            C_ks_tomo_limber_wrapper(l, nz) :
            C_ks_tomo_limber_nointerp_wrapper(l, nz, use_linear_ps_limber, 0);
        }
      }
      #pragma omp parallel for collapse(2)
      for (int nz=0; nz<tomo.shear_Nbin; nz++)
      {
        for (int i=0; i<ntheta; i++)
        {
          w_vec[nz*ntheta+i] = 0;
          for (int l=2; l<nell; l++)
          {
            w_vec[nz*ntheta+i] += Pl[i][l]*Cl[nz][l];
          }
        }
      }
      for (int nz=0; nz<tomo.shear_Nbin; nz++)
      {
        free(Cl[nz]);
      }
      free(Cl);
    } 
    else
    {
      log_fatal("NonLimber not implemented");
      exit(1);
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  if (nt < 0 || nt > like.Ntheta - 1)
  {
    log_fatal("error in selecting bin number nt = %d (max %d)", nt, like.Ntheta);
    exit(1); 
  }
  if (ni < -1 || ni > tomo.shear_Nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d (max %d)", ni, tomo.shear_Nbin);
    exit(1);
  }
  const int q = ni * ntheta + nt;
  if (q  < 0 || q > NSIZE*ntheta - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  }  
  return w_vec[q];
}

// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
// Correlation Functions (real space) - flat sky
// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------

double xi_pm_tomo_flatsky(int pm, double theta, int ni, int nj, int limber)
{
  static cosmopara C;
  static nuisancepara N;
  static double** table;

  const int NSIZE = tomo.shear_Npowerspectra;
  const int ntheta = Ntable.N_thetaH;
  const double l_min = limits.LMIN_hankel;
  const double l_max = limits.LMAX_hankel;
  const double lnlmax = log(l_max);
  const double lnlmin = log(l_min);
  const double dlnl = (lnlmax - lnlmin)/(1.0*ntheta);
  const double lnrc = 0.5*(lnlmax + lnlmin);
  const double nc = ntheta/2 + 1;
  
  const double lnthetamin = (nc-ntheta+1)*dlnl-lnrc;
  const double lnthetamax = nc*dlnl - lnrc;
  const double dlntheta = (lnthetamax - lnthetamin)/((double) ntheta);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*2*NSIZE);  // 2 NSIZE = {xi+, xi-}
    for (int i=0; i<2*NSIZE; i++) 
    {
      table[i] = (double*) calloc(ntheta, sizeof(double));
    }
  }
  if (recompute_shear(C, N))
  {
    typedef fftw_complex fftwZ;

    if (limber != 1)
    {
      log_fatal("NonLimber not implemented");
      exit(1);
    }
    else
    {
      if (like.IA == 5 || like.IA == 6)
      { // NEW TATT MODELING
        log_fatal("Limber && TATT not implemented");
        exit(1);
      }
      else
      {
        fftwZ** flP = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
        for (int j=0; j<NSIZE; j++)
        {
          flP[j] = (fftwZ*) fftw_malloc((Ntable.N_thetaH/2+1)*sizeof(fftwZ));
        }
        {
          double** lP = (double**) malloc(sizeof(double*)*NSIZE);
          fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
          for (int j=0; j<NSIZE; j++)
          {
            lP[j] = (double*) malloc(ntheta*sizeof(double));
            plan[j] = fftw_plan_dft_r2c_1d(ntheta,lP[j],flP[j],FFTW_ESTIMATE);
          }
          
          // --------------------------------------------------------------------------------
          // Power spectrum on logarithmic bins (begins)
          // --------------------------------------------------------------------------------
          if (like.IA == 5 || like.IA == 6)
          { // NEW TATT MODELING
            log_fatal("Limber && TATT not implemented");
            exit(1);
          }
          else
          {
            #pragma GCC diagnostic push
            #pragma GCC diagnostic ignored "-Wunused-variable"
            {
              double init_static_vars_only = C_ss_tomo_limber(limits.LMIN_tab + 1, Z1(0), Z2(0));
            }
            #pragma GCC diagnostic pop
            #pragma omp parallel for collapse(2)
            for (int k=0; k<NSIZE; k++)
            {
              for (int p=0; p<ntheta; p++)
              {
                const int Z1NZ = Z1(k);
                const int Z2NZ = Z2(k);
                const double l = exp(lnrc+(p - nc)*dlnl);
                lP[k][p] = (l > limits.LMIN_tab) ?
                  l*C_ss_tomo_limber(l, Z1NZ, Z2NZ) :
                  l*C_ss_tomo_limber_nointerp(l, Z1NZ, Z2NZ, use_linear_ps_limber, 0);
              }
            }
          }
          // --------------------------------------------------------------------------------
          // Power spectrum on logarithmic bins (ends)
          // --------------------------------------------------------------------------------
          #pragma omp parallel for
          for (int j=0; j<NSIZE; j++)
          {
            fftw_execute(plan[j]); // Execute FFTW in parallel (thread-safe)
          }
          for (int j = 0; j < NSIZE; j++)
          {
            fftw_free(lP[j]);
            fftw_destroy_plan(plan[j]);
          }
          free(lP);
          free(plan);
        }
        
        double*** lP = (double***) malloc(sizeof(double**)*NSIZE);
        fftwZ*** kernel = (fftwZ***) malloc(sizeof(fftwZ**)*NSIZE);
        fftwZ*** conv = (fftwZ***) malloc(sizeof(fftwZ**)*NSIZE);
        fftw_plan** plan = (fftw_plan**) malloc(sizeof(fftw_plan*)*NSIZE);
        double*** tab = (double***) malloc(sizeof(double**)*NSIZE);
        for (int j=0; j<NSIZE; j++)
        {
          lP[j] = (double**) malloc(sizeof(double*)*2);
          kernel[j] = (fftwZ**) malloc(sizeof(fftwZ*)*2);
          conv[j] = (fftwZ**) malloc(sizeof(fftwZ*)*2);
          plan[j] = (fftw_plan*) malloc(sizeof(fftw_plan)*2);
          tab[j] = (double**) malloc(sizeof(double*)*2);
          tab[j][0] = (double*) malloc(sizeof(double)*ntheta);
          tab[j][1] = (double*) malloc(sizeof(double)*ntheta);
          for (int m=0; m<2; m++)
          {
            const int COVSZ = (ntheta/2+1);
            lP[j][m] = (double*) malloc(ntheta*sizeof(double));
            kernel[j][m] = (fftwZ*) fftw_malloc(sizeof(fftwZ));
            conv[j][m] = (fftwZ*) fftw_malloc(COVSZ*sizeof(fftwZ));
            plan[j][m] = fftw_plan_dft_c2r_1d(ntheta,conv[j][m],lP[j][m],FFTW_ESTIMATE);
          }
        }
        #pragma omp parallel for
        for (int j=0; j<NSIZE; j++)
        {
          for (int m=0; m<2; m++)
          {
            double arg[2];
            arg[0] = 0; // bias
            arg[1] = (m == 0 ? 0 : 4); // order of Bessel function

            // perform the convolution, negative sign for kernel (complex conj.)
            for(int i=0; i<(ntheta/2+1); i++)
            {
              const double k = 2*M_PI*i/(dlnl*ntheta);
              hankel_kernel_FT(k, kernel[j][m], arg, 2);
              conv[j][m][i][0] = flP[j][i][0]*(kernel[j][m][0][0]) - flP[j][i][1]*(kernel[j][m][0][1]);
              conv[j][m][i][1] = flP[j][i][1]*(kernel[j][m][0][0]) + flP[j][i][0]*(kernel[j][m][0][1]);
            }

            // force Nyquist- and 0-frequency-components to be double
            conv[j][m][0][1] = 0;
            conv[j][m][ntheta/2][1] = 0;

            fftw_execute(plan[j][m]);

            for(int k=0; k<ntheta; k++)
            {
              const double t = exp((nc-k)*dlnl-lnrc); // theta=1/l
              tab[j][m][ntheta-k-1] = lP[j][m][k]/(t*2*M_PI*ntheta);
            }
          }
          for (int k=0; k<ntheta; k++)
          {
            table[2*j][k] = tab[j][0][k];
            table[2*j+1][k] = tab[j][1][k];
          }
        }
        for (int j=0; j<NSIZE; j++)
        {
          for (int m=0; m<2; m++)
          {
            fftw_free(lP[j][m]);
            fftw_free(kernel[j][m]);
            fftw_free(conv[j][m]);
            fftw_destroy_plan(plan[j][m]);
          }
          fftw_free(flP[j]);
          free(lP[j]);
          free(kernel[j]);
          free(conv[j]);
          free(plan[j]);
          free(tab[j][0]);
          free(tab[j][1]);
          free(tab[j]);
        }
        free(flP);
        free(lP);
        free(kernel);
        free(conv);
        free(plan);
        free(tab);
        // --------------------------------------------------------------------
        // Cocoa: code extracted (& adapted) from xipm_via_hankel (ends)
        // --------------------------------------------------------------------
      }
    } 
    update_cosmopara(&C);
    update_nuisance(&N);
  }

  if(ni < -1 || ni > tomo.shear_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  const int q = 2*N_shear(ni, nj) + (1 - pm)/2;
  if(q > 2*NSIZE - 1)
  {
    log_fatal("error in selecting bin number, q = %d and max bin = %d", q, 2*NSIZE);
    exit(1);
  }

  const double lntheta = log(theta);
  if (lntheta < lnthetamin || lntheta > lnthetamax)
  {
    const double theta = exp(lntheta);
    const double theta_min = exp(lnthetamin);
    const double theta_max = exp(lnthetamax);
    log_fatal("theta = %e outside look-up table range [%e, %e]", theta, theta_min, theta_max);
    exit(1);
  }

  return interpol(table[q], ntheta, lnthetamin, lnthetamax, dlntheta, lntheta, 0, 0);
}

// ---------------------------------------------------------------------------

double w_gammat_tomo_flatsky(const double theta, const int ni, const int nj, const int limber)
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double** table;

  const int NSIZE = tomo.ggl_Npowerspectra;
  const int ntheta = Ntable.N_thetaH;
  const double l_min = limits.LMIN_hankel;
  const double l_max = limits.LMAX_hankel;
  const double lnlmax = log(l_max);
  const double lnlmin = log(l_min);
  const double dlnl = (lnlmax - lnlmin)/((double) ntheta - 1.0);
  const double lnrc = 0.5*(lnlmax + lnlmin);
  const double nc = ntheta/2 + 1;
  const double lnthetamin = (nc - ntheta + 1)*dlnl - lnrc;
  const double lnthetamax = nc*dlnl - lnrc;
  const double dlntheta = (lnthetamax - lnthetamin)/((double) ntheta);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for (int i=0; i<NSIZE; i++) 
    {
      table[i] = (double*) malloc(sizeof(double)*ntheta);
    }
  }
  if (recompute_gs(C, G, N)) 
  {
    typedef fftw_complex fftwZ;

    fftwZ** flP = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE); // go to log-Fourier-space
    for (int j=0; j<NSIZE; j++)
    {
      flP[j] = fftw_malloc((ntheta/2 + 1)*sizeof(fftwZ));
    }
    {
      double** lP = (double**) malloc(sizeof(double*)*NSIZE);
      fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
      for (int j=0; j<NSIZE; j++)
      {
        lP[j] = (double*) malloc(ntheta*sizeof(double));
        plan[j] = fftw_plan_dft_r2c_1d(ntheta, lP[j], flP[j], FFTW_ESTIMATE);
      }
      // --------------------------------------------------------------------------------
      // Power spectrum on logarithmic bins (begins) 
      // --------------------------------------------------------------------------------
      gsl_spline** fCL_NL = NULL; 
      if (limber != 1)
      { 
        fCL_NL = (gsl_spline**) malloc(sizeof(gsl_spline*)*NSIZE);
        double** Cl_NL = (double**) malloc(sizeof(double*)*NSIZE);
        for (int i=0; i<NSIZE; i++)
        {
          Cl_NL[i] = calloc(limits.LMAX_NOLIMBER, sizeof(double));
        }
        double* ll = calloc(limits.LMAX_NOLIMBER, sizeof(double));
        for (int i=0; i<limits.LMAX_NOLIMBER; i++)
        {
          ll[i] = i;
        }

        const int L = 1;
        const double tolerance = 0.0075;    // required fractional accuracy in C(l)
        const double dev = 10. * tolerance; // will be diff exact vs Limber init to large
                                            // value in order to start while loop
        for (int k=0; k<NSIZE; k++)
        { // Cocoa: no threading allowed here - (fftw allocation @C_gl_tomo)
          const int Z1 = ZL(k);
          const int Z2 = ZS(k);
          C_gl_tomo(L, Z1, Z2, Cl_NL[k], dev, tolerance);

          const gsl_interp_type* T = gsl_interp_linear;
          fCL_NL[k] = gsl_spline_alloc(T, limits.LMAX_NOLIMBER);
          if (fCL_NL[k] == NULL)
          {
            log_fatal("fail allocation");
            exit(1);
          }
        }
        #pragma omp parallel for
        for (int k=0; k<NSIZE; k++)
        {
          int status = gsl_spline_init(fCL_NL[k], ll, Cl_NL[k], limits.LMAX_NOLIMBER);
          if (status) 
          {
            log_fatal(gsl_strerror(status));
            exit(1);
          }
        }
        for (int i=0; i<NSIZE; i++)
        {
          free(Cl_NL[i]);   // gsl_spline_init copies the data
        }
        free(Cl_NL);        // gsl_spline_init copies the data
        free(ll);           // gsl_spline_init copies the data
      }
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wunused-variable"
      {
        double init_static_vars_only = C_gs_tomo_limber(limits.LMIN_tab + 1, ZL(0), ZS(0));
      }
      #pragma GCC diagnostic pop
      #pragma omp parallel for collapse(2)
      for (int k=0; k<NSIZE; k++)
      {
        for (int p=0; p<ntheta; p++)
        {
          const int ZLNZ = ZL(k);
          const int ZSNZ = ZS(k);
          const double l = exp(lnrc + (p - nc)*dlnl);
          if (limber == 1 || (limber != 1 && l > limits.LMAX_NOLIMBER - 1))
          {
            lP[k][p] = (l > limits.LMIN_tab) ?
              l*C_gs_tomo_limber(l, ZLNZ, ZSNZ) :
              l*C_gs_tomo_limber_nointerp(l, ZLNZ, ZSNZ, use_linear_ps_limber, 0);
          }
          else
          {
            double CLNL;
            int status = gsl_spline_eval_e(fCL_NL[k], l, NULL, &CLNL);
            if (status) 
            {
              log_fatal(gsl_strerror(status));
              exit(1);
            }
            lP[k][p] = l*CLNL;
          }
        }
      }
      if (limber != 1)
      {
        for (int i=0; i<NSIZE; i++)
        {
           gsl_spline_free(fCL_NL[i]);
        }
        free(fCL_NL);
      }
      // --------------------------------------------------------------------------------
      // Power spectrum on logarithmic bins (ends)
      // --------------------------------------------------------------------------------
      #pragma omp parallel for
      for (int j=0; j<NSIZE; j++)
      {
        fftw_execute(plan[j]); // Execute FFTW in parallel (thread-safe)
      }
      for (int j=0; j<NSIZE; j++)
      {
        fftw_free(lP[j]);
        fftw_destroy_plan(plan[j]);
      }
      free(lP);
      free(plan);
    }

    double** lP = (double**) malloc(sizeof(double*)*NSIZE);
    fftwZ** kernel = (fftwZ**) fftw_malloc(sizeof(fftwZ*)*NSIZE);
    fftwZ** conv = (fftwZ**) fftw_malloc(sizeof(fftwZ*)*NSIZE);
    fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
    double*** tab = (double***) malloc(sizeof(double**)*NSIZE);
    for (int j=0; j<NSIZE; j++)
    {
      const int COVSZ = (ntheta/2+1);
      lP[j] = (double*) malloc(ntheta*sizeof(double));
      kernel[j] = (fftwZ*) fftw_malloc(sizeof(fftwZ));
      conv[j] = (fftwZ*) fftw_malloc(COVSZ*sizeof(fftwZ));
      plan[j] = fftw_plan_dft_c2r_1d(ntheta, conv[j], lP[j], FFTW_ESTIMATE);
      tab[j] = (double**) malloc(sizeof(double*)*1);
      tab[j][0] = (double*) malloc(sizeof(double)*ntheta);
    }
    #pragma omp parallel for
    for (int j=0; j<NSIZE; j++)
    {
      double arg[2];
      arg[0] = 0; // bias
      arg[1] = 2; // order of Bessel function

      for (int i=0; i<(ntheta/2+1); i++)
      {
        const double kk = 2*M_PI*i/(dlnl*ntheta);
        hankel_kernel_FT(kk, kernel[j], arg, 2);
        conv[j][i][0] = flP[j][i][0]*kernel[j][0][0] - flP[j][i][1]*kernel[j][0][1];
        conv[j][i][1] = flP[j][i][1]*kernel[j][0][0] + flP[j][i][0]*kernel[j][0][1];
      }

      // force Nyquist- and 0-frequency-components to be double
      conv[j][0][1] = 0;
      conv[j][ntheta/2][1] = 0;

      fftw_execute(plan[j]); // Execute FFTW in parallel (thread-safe)

      for (int k=0; k<ntheta; k++)
      {
        const double t = exp((nc - k)*dlnl - lnrc); // theta=1/l
        tab[j][0][ntheta-k-1] = lP[j][k]/(t*2*M_PI*ntheta);
      }
      
      for (int k=0; k<ntheta; k++)
      {
        table[j][k] = tab[j][0][k];
      }
    }
    for (int j=0; j<NSIZE; j++)
    {
      fftw_free(flP[j]);
      fftw_free(lP[j]);
      fftw_free(conv[j]);
      fftw_free(kernel[j]);
      fftw_destroy_plan(plan[j]);
      free(tab[j][0]);
      free(tab[j]);
    }
    free(flP);
    free(lP);
    free(conv);
    free(kernel);
    free(plan);
    free(tab);

    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  if (ni < -1 || ni > tomo.clustering_Nbin - 1 || nj < -1 || nj > tomo.shear_Nbin - 1)
  {
    log_fatal("error in selecting bin number (ni, nj) = [%d,%d]", ni, nj);
    exit(1);
  }
  const double lntheta = log(theta);
  if (lntheta < lnthetamin || lntheta > lnthetamax)
  {
    const double theta = exp(lntheta);
    const double theta_min = exp(lnthetamin);
    const double theta_max = exp(lnthetamax);
    log_fatal("theta = %e outside look-up table range [%e, %e]", theta, theta_min, theta_max);
    exit(1);
  }
  if (test_zoverlap(ni, nj))
  {
    const int q = N_ggl(ni, nj);
    if (q  < 0 || q > NSIZE - 1)
    {
      log_fatal("internal logic error in selecting bin number");
      exit(1);
    }
    return interpol(table[q], ntheta, lnthetamin, lnthetamax, dlntheta, lntheta, 1, 1);
  } 
  else
  {
    return 0.0;
  }
}

// ---------------------------------------------------------------------------

double w_gg_tomo_flatsky(double theta, int ni, int nj, int limber)
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double** table;

  const int NSIZE = tomo.clustering_Npowerspectra;
  const int ntheta = Ntable.N_thetaH;
  const double l_min = limits.LMIN_hankel;
  const double l_max = limits.LMAX_hankel;
  const double lnlmax = log(l_max);
  const double lnlmin = log(l_min);
  const double dlnl = (lnlmax - lnlmin)/((double) ntheta - 1.0);
  const double lnrc = 0.5*(lnlmax + lnlmin);
  const double nc = ntheta/2 + 1;
  const double lnthetamin = (nc - ntheta + 1)*dlnl - lnrc;
  const double lnthetamax = nc*dlnl - lnrc;
  const double dlntheta = (lnthetamax - lnthetamin)/((double) ntheta);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for (int i=0; i<NSIZE; i++) 
    {
      table[i] = (double*) malloc(sizeof(double)*ntheta);
    }
  }  
  if (recompute_gg(C, G, N)) 
  {
    typedef fftw_complex fftwZ;

    fftwZ** flP = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE); // go to log-Fourier-space
    for (int j=0; j<NSIZE; j++)
    {
      flP[j] = fftw_malloc((ntheta/2 + 1)*sizeof(fftwZ));
    }
    {
      double** lP = (double**) malloc(sizeof(double*)*NSIZE);
      fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
      for (int j=0; j<NSIZE; j++)
      {
        int ARRAYSZ = ntheta;
        lP[j] = (double*) malloc(ARRAYSZ*sizeof(double));
        plan[j] = fftw_plan_dft_r2c_1d(ARRAYSZ, lP[j], flP[j], FFTW_ESTIMATE);
      }
      // --------------------------------------------------------------------------------
      // Power spectrum on logarithmic bins (begins)
      // --------------------------------------------------------------------------------
      { 
        gsl_spline** fCL_NL = NULL; 
        if (limber != 1)
        { 
          fCL_NL = (gsl_spline**) malloc(sizeof(gsl_spline*)*NSIZE);
          double** Cl_NL = (double**) malloc(sizeof(double*)*NSIZE);
          for (int i=0; i<NSIZE; i++)
          {
            Cl_NL[i] = calloc(limits.LMAX_NOLIMBER, sizeof(double));
          }
          double* ll = calloc(limits.LMAX_NOLIMBER, sizeof(double));
          for (int i=0; i<limits.LMAX_NOLIMBER; i++)
          {
            ll[i] = i;
          }

          const int L = 1;
          const double tolerance = 0.0075;    // required fractional accuracy in C(l)
          const double dev = 10. * tolerance; // will be diff exact vs Limber init to large
                                              // value in order to start while loop
          for (int k=0; k<NSIZE; k++)
          { // Cocoa: no threading allowed here - (fftw allocation @C_gl_tomo)
            const int ZCL1 = k; // cross redshift bin not supported so not using ZCL1(k)
            const int ZCL2 = k; // cross redshift bin not supported so not using ZCL2(k) 
            C_cl_tomo(L, ZCL1, ZCL2, Cl_NL[k], dev, tolerance);

            const gsl_interp_type* T = gsl_interp_linear;
            fCL_NL[k] = gsl_spline_alloc(T, limits.LMAX_NOLIMBER);
            if (fCL_NL[k] == NULL)
            {
              log_fatal("fail allocation");
              exit(1);
            }
          }
          #pragma omp parallel for
          for (int k=0; k<NSIZE; k++)
          {
            int status = gsl_spline_init(fCL_NL[k], ll, Cl_NL[k], limits.LMAX_NOLIMBER);
            if (status) 
            {
              log_fatal(gsl_strerror(status));
              exit(1);
            }
          } 
          for (int i=0; i<NSIZE; i++)
          {
            free(Cl_NL[i]);
          }
          free(Cl_NL);
          free(ll);
        } 
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-variable"
        {
          double init_static_vars_only = C_gg_tomo_limber(limits.LMIN_tab + 1, 0, 0);
        }
        #pragma GCC diagnostic pop
        #pragma omp parallel for collapse(2)
        for (int k=0; k<NSIZE; k++) 
        { 
          for (int p=0; p<ntheta; p++)
          {
            const int ZCL1 = k; // cross redshift bin not supported so not using ZCL1(k)
            const int ZCL2 = k; // cross redshift bin not supported so not using ZCL2(k) 
            const double l = exp(lnrc + (p - nc)*dlnl);
            if (limber == 1 || (limber != 1 && l > limits.LMAX_NOLIMBER - 1))
            {
              lP[k][p] = (l > limits.LMIN_tab) ?
                l*C_gg_tomo_limber(l, ZCL1, ZCL2) :
                l*C_gg_tomo_limber_nointerp(l, ZCL1, ZCL2, use_linear_ps_limber, 0);
            }
            else
            {
              double CLNL;
              int status = gsl_spline_eval_e(fCL_NL[k], l, NULL, &CLNL);
              if (status) 
              {
                log_fatal(gsl_strerror(status));
                exit(1);
              }
              lP[k][p] = l*CLNL;
            }
          }
        }
        if (limber != 1)
        {
          for (int i=0; i<NSIZE; i++)
          {
             gsl_spline_free(fCL_NL[i]);
          }
          free(fCL_NL);
        }
      } 
      // --------------------------------------------------------------------------------
      // Power spectrum on logarithmic bins (ends)
      // --------------------------------------------------------------------------------
      #pragma omp parallel for
      for (int j=0; j<NSIZE; j++)
      {
        fftw_execute(plan[j]); // Execute FFTW in parallel (thread-safe)
      }
      for (int j=0; j<NSIZE; j++)
      {
        fftw_free(lP[j]);
        fftw_destroy_plan(plan[j]);
      }
      free(lP);
      free(plan);
    }

    double** lP = (double**) malloc(sizeof(double*)*NSIZE);
    fftwZ** kernel = (fftwZ**) fftw_malloc(sizeof(fftwZ*)*NSIZE);
    fftwZ** conv = (fftwZ**) fftw_malloc(sizeof(fftwZ*)*NSIZE);
    fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
    double*** tab = (double***) malloc(sizeof(double**)*NSIZE);
    for (int j=0; j<NSIZE; j++)
    {
      const int ARRAYSZ = ntheta;
      const int COVSZ = (ntheta/2+1);
      lP[j] = (double*) malloc(ARRAYSZ*sizeof(double));
      kernel[j] = (fftwZ*) fftw_malloc(sizeof(fftwZ));
      conv[j] = (fftwZ*) fftw_malloc(COVSZ*sizeof(fftwZ));
      plan[j] = fftw_plan_dft_c2r_1d(ARRAYSZ, conv[j], lP[j], FFTW_ESTIMATE);
      tab[j] = (double**) malloc(sizeof(double*)*1);
      tab[j][0] = (double*) malloc(sizeof(double)*ntheta);
    }
    #pragma omp parallel for
    for (int j=0; j<NSIZE; j++)
    {
      double arg[2];
      arg[0] = 0; // bias
      arg[1] = 0; // order of Bessel function

      for (int i=0; i<(ntheta/2+1); i++)
      {
        const double kk = 2*M_PI*i/(dlnl*ntheta);
        hankel_kernel_FT(kk, kernel[j], arg, 2);
        conv[j][i][0] = flP[j][i][0]*kernel[j][0][0] - flP[j][i][1]*kernel[j][0][1];
        conv[j][i][1] = flP[j][i][1]*kernel[j][0][0] + flP[j][i][0]*kernel[j][0][1];
      }

      // force Nyquist- and 0-frequency-components to be double
      conv[j][0][1] = 0;
      conv[j][ntheta/2][1] = 0;

      fftw_execute(plan[j]); // Execute FFTW in parallel (thread-safe)

      for (int k=0; k<ntheta; k++)
      {
        const double t = exp((nc - k)*dlnl - lnrc); // theta=1/l
        tab[j][0][ntheta-k-1] = lP[j][k]/(t*2*M_PI*ntheta);
      }
      for (int k=0; k<ntheta; k++)
      {
        table[j][k] = tab[j][0][k];
      }
    }
    for (int j=0; j<NSIZE; j++)
    {
      fftw_free(flP[j]);
      fftw_free(lP[j]);
      fftw_free(conv[j]);
      fftw_free(kernel[j]);
      fftw_destroy_plan(plan[j]);
      free(tab[j][0]);
      free(tab[j]);

    }
    free(flP);
    free(lP);
    free(conv);
    free(kernel);
    free(plan);
    free(tab);
    
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  if (ni != nj) 
  {
    log_fatal("cross-tomography not supported");
    exit(1);
  }
  if (ni < -1 || ni > tomo.clustering_Nbin - 1 || nj < -1 || nj > tomo.clustering_Nbin - 1)
  {
    log_fatal("error in selecting bin number (ni, nj) = [%d,%d]", ni, nj);
    exit(1);
  }
  const double lntheta = log(theta);
  if (lntheta < lnthetamin || lntheta > lnthetamax)
  {
    const double theta = exp(lntheta);
    const double theta_min = exp(lnthetamin);
    const double theta_max = exp(lnthetamax);
    log_fatal("theta = %e outside look-up table range [%e, %e]", theta, theta_min, theta_max);
    exit(1);
  }
  const int q = ni; // cross redshift bin not supported so not using N_CL(ni, nj) instead of ni
  if (q < 0 || q > NSIZE - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  }
  return interpol(table[q], ntheta, lnthetamin, lnthetamax, dlntheta, lntheta, 1, 1);
}

// ---------------------------------------------------------------------------

double w_gk_tomo_flatsky(const double theta, const int ni, const int limber)
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double** table;

  const int nisize = tomo.clustering_Nbin;
  const int ntheta = Ntable.N_thetaH;
  const int NSIZE = nisize;
  const double l_min = limits.LMIN_hankel;
  const double l_max = limits.LMAX_hankel;
  const double lnlmax = log(l_max);
  const double lnlmin = log(l_min);
  const double dlnl = (lnlmax - lnlmin)/((double) ntheta - 1.0);
  const double lnrc = 0.5*(lnlmax + lnlmin);
  const double nc = ntheta/2 + 1;
  const double lnthetamin = (nc - ntheta + 1)*dlnl - lnrc;
  const double lnthetamax = nc*dlnl - lnrc;
  const double dlntheta = (lnthetamax-lnthetamin)/((double)ntheta);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for (int i=0; i<NSIZE; i++) 
    {
      table[i] = (double*) malloc(sizeof(double)*ntheta);
    }
  }
  if (recompute_gk(C, G, N))
  {
    if (limber == 1)
    {
      typedef fftw_complex fftwZ;

      // go to log-Fourier-space
      fftwZ** flP = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
      for (int i=0; i<NSIZE; i++)
      {
        flP[i] = fftw_malloc((ntheta/2 + 1)*sizeof(fftwZ));
      }
      {
        double** lP = (double**) malloc(sizeof(double*)*NSIZE);
        fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
        for (int i=0; i<NSIZE; i++)
        {
          lP[i] = (double*) malloc(ntheta*sizeof(double));
          plan[i] = fftw_plan_dft_r2c_1d(ntheta, lP[i], flP[i], FFTW_ESTIMATE);
        }
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-variable"
        {
          double init_static_vars_only = C_gk_tomo_limber_wrapper(limits.LMIN_tab + 1, 0);
        }
        #pragma GCC diagnostic pop 
        // --------------------------------------------------------------------------------
        // Power spectrum on logarithmic bins (begins)
        // --------------------------------------------------------------------------------
        #pragma omp parallel for collapse(2)
        for (int i=0; i<NSIZE; i++)
        {
          for (int p=0; p<ntheta; p++)
          {
            const double l = exp(lnrc + (p - nc)*dlnl);
            lP[i][p] = (l > limits.LMIN_tab) ?
              l*C_gk_tomo_limber_wrapper(l, i) :
              l*C_gk_tomo_limber_nointerp_wrapper(l, i, use_linear_ps_limber, 0);
          }
        }
        // --------------------------------------------------------------------------------
        // Power spectrum on logarithmic bins (ends)
        // --------------------------------------------------------------------------------
        #pragma omp parallel for
        for (int i=0; i<NSIZE; i++)
        {
          fftw_execute(plan[i]); // Execute FFTW in parallel (thread-safe)
        }
        for (int i=0; i<NSIZE; i++)
        {
          fftw_free(lP[i]);
          fftw_destroy_plan(plan[i]);
        } 
        free(lP);
        free(plan);
      }

      double** lP = (double**) fftw_malloc(sizeof(double**)*NSIZE);
      fftwZ** kernel = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
      fftwZ** conv = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
      fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
      double*** tab = (double***) malloc(sizeof(double**)*NSIZE);
      for (int i=0; i<NSIZE; i++)
      {
        lP[i] = (double*) fftw_malloc(ntheta*sizeof(double));
        kernel[i] = (fftwZ*) fftw_malloc(sizeof(fftwZ));
        conv[i] = (fftwZ*) fftw_malloc((ntheta/2 + 1)*sizeof(fftwZ));
        plan[i] = fftw_plan_dft_c2r_1d(ntheta, conv[i], lP[i], FFTW_ESTIMATE);
        tab[i] = (double**) malloc(sizeof(double*)*1);
        tab[i][0] = (double*) malloc(sizeof(double)*ntheta);
      }
      #pragma omp parallel for
      for (int i=0; i<NSIZE; i++)
      {
        double arg[2];
        arg[0] = 0; // bias
        arg[1] = 0; // order of Bessel function

        for (int k=0; k<(ntheta/2+1); k++)
        {
          const double kk = 2*M_PI*i/(dlnl*ntheta);
          hankel_kernel_FT(kk, kernel[i], arg, 2);
          
          conv[i][k][0] = flP[i][k][0]*kernel[i][0][0] - flP[i][k][1]*kernel[i][0][1];
          conv[i][k][1] = flP[i][k][1]*kernel[i][0][0] + flP[i][k][0]*kernel[i][0][1];
        }

        // force Nyquist- and 0-frequency-components to be double
        conv[i][0][1] = 0;
        conv[i][ntheta/2][1] = 0;

        fftw_execute(plan[i]); // Execute FFTW in parallel (thread-safe)

        for (int k=0; k<ntheta; k++)
        {
          const double t = exp((nc - k)*dlnl - lnrc); // theta=1/l
          tab[i][0][ntheta-k-1] = lP[i][k]/(t*2*M_PI*ntheta);
        }
        for (int k=0; k<ntheta; k++)
        {
          table[i][k] = tab[i][0][k];
        }
      }
      for (int i=0; i<NSIZE; i++)
      {
        fftw_free(flP[i]);
        fftw_free(lP[i]);
        fftw_free(kernel[i]);
        fftw_free(conv[i]);
        fftw_destroy_plan(plan[i]);
        free(tab[i][0]);
        free(tab[i]);
      }
      free(flP);
      free(lP);
      free(kernel);
      free(conv);
      free(plan);
      free(tab);
    }
    else
    {
      log_fatal("NonLimber not implemented");
      exit(1);
    }
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  if (ni < -1 || ni > tomo.clustering_Nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }
  const double lntheta = log(theta);
  if (lntheta < lnthetamin || lntheta > lnthetamax)
  {
    const double theta = exp(lntheta);
    const double theta_min = exp(lnthetamin);
    const double theta_max = exp(lnthetamax);
    log_fatal("theta = %e outside look-up table range [%e, %e]", theta, theta_min, theta_max);
    exit(1);
  }
  const int q = ni;
  if (q < 0 || q > NSIZE - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  }
  return interpol(table[q], ntheta, lnthetamin, lnthetamax, dlntheta, lntheta, 1, 1);
}

// ---------------------------------------------------------------------------

double w_ks_tomo_flatsky(const double theta, const int ni, const int limber)
{
  static cosmopara C;
  static nuisancepara N;
  static double** table;

  const int nisize = tomo.shear_Nbin;
  const int ntheta = Ntable.N_thetaH;
  const int NSIZE = nisize;
  const double l_min = limits.LMIN_hankel;
  const double l_max = limits.LMAX_hankel;
  const double lnlmax = log(l_max);
  const double lnlmin = log(l_min);
  const double dlnl = (lnlmax - lnlmin)/((double) ntheta - 1.0);
  const double lnrc = 0.5*(lnlmax + lnlmin);
  const double nc = ntheta/2 + 1;
  const double lnthetamin = (nc - ntheta + 1)*dlnl - lnrc;
  const double lnthetamax = nc*dlnl - lnrc;
  const double dlntheta = (lnthetamax - lnthetamin)/((double) ntheta);

  if (table == 0) 
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for (int i=0; i<NSIZE; i++) 
    {
      table[i] = (double*) malloc(sizeof(double)*ntheta);
    }
  }
  if (recompute_shear(C, N))
  {
    if (limber == 1)
    {
      typedef fftw_complex fftwZ;

      // go to log-Fourier-space
      fftwZ** flP = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
      for (int i=0; i<NSIZE; i++)
      {
        flP[i] = fftw_malloc((ntheta/2 + 1)*sizeof(fftwZ));
      }
      {
        double** lP = (double**) malloc(sizeof(double*)*NSIZE);
        fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
        for (int i=0; i<NSIZE; i++)
        {
          lP[i] = (double*) malloc(ntheta*sizeof(double));
          plan[i] = fftw_plan_dft_r2c_1d(ntheta, lP[i], flP[i], FFTW_ESTIMATE);
        }
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-variable"
        {
          double init_static_vars_only = C_ks_tomo_limber_wrapper(limits.LMIN_tab + 1, 0);
        }
        #pragma GCC diagnostic pop        
        // --------------------------------------------------------------------------------
        // Power spectrum on logarithmic bins (begins)
        // --------------------------------------------------------------------------------
        #pragma omp parallel for collapse(2)
        for (int i=0; i<NSIZE; i++)
        {
          for (int p=0; p<ntheta; p++)
          {
            const double l = exp(lnrc + (p - nc)*dlnl);
            lP[i][p] = (l > limits.LMIN_tab) ?
              l*C_ks_tomo_limber_wrapper(l, i) :
              l*C_ks_tomo_limber_nointerp_wrapper(l, i, use_linear_ps_limber, 0);
          }
        }
        // --------------------------------------------------------------------------------
        // Power spectrum on logarithmic bins (ends)
        // --------------------------------------------------------------------------------
        #pragma omp parallel for
        for (int i=0; i<NSIZE; i++)
        {
          fftw_execute(plan[i]); // Execute FFTW in parallel (thread-safe)
        }
        for (int i=0; i<NSIZE; i++)
        {
          fftw_free(lP[i]);
          fftw_destroy_plan(plan[i]);
        }
        free(lP);
        free(plan);
      }

      double** lP = (double**) fftw_malloc(sizeof(double**)*NSIZE);
      fftwZ** kernel = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
      fftwZ** conv = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
      fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
      double*** tab = (double***) malloc(sizeof(double**)*NSIZE);
      for (int i=0; i<NSIZE; i++)
      {
        lP[i] = (double*) fftw_malloc(ntheta*sizeof(double));
        kernel[i] = (fftwZ*) fftw_malloc(sizeof(fftwZ));
        conv[i] = (fftwZ*) fftw_malloc((ntheta/2 + 1)*sizeof(fftwZ));
        plan[i] = fftw_plan_dft_c2r_1d(ntheta, conv[i], lP[i], FFTW_ESTIMATE);
        tab[i] = (double**) malloc(sizeof(double*)*1);
        tab[i][0] = (double*) malloc(sizeof(double)*ntheta);
      }
      
      #pragma omp parallel for
      for (int i=0; i<NSIZE; i++)
      {
        double arg[2];
        arg[0] = 0; // bias
        arg[1] = 2; // order of Bessel function

        for (int k=0; k<(ntheta/2+1); k++)
        {
          const double kk = 2*M_PI*i/(dlnl*ntheta);
          hankel_kernel_FT(kk, kernel[i], arg, 2);
          conv[i][k][0] = flP[i][k][0]*kernel[i][0][0] - flP[i][k][1]*kernel[i][0][1];
          conv[i][k][1] = flP[i][k][1]*kernel[i][0][0] + flP[i][k][0]*kernel[i][0][1];
        }

        // force Nyquist- and 0-frequency-components to be double
        conv[i][0][1] = 0;
        conv[i][ntheta/2][1] = 0;

        fftw_execute(plan[i]); // Execute FFTW in parallel (thread-safe)

        for (int k=0; k<ntheta; k++)
        {
          const double t = exp((nc - k)*dlnl - lnrc); // theta=1/l
          tab[i][0][ntheta-k-1] = lP[i][k]/(t*2*M_PI*ntheta);
        }
        for (int k=0; k<ntheta; k++)
        {
          table[i][k] = tab[i][0][k];
        }
      }
      for (int i=0; i<NSIZE; i++)
      {
        fftw_free(flP[i]);
        fftw_free(lP[i]);
        fftw_free(kernel[i]);
        fftw_free(conv[i]);
        fftw_destroy_plan(plan[i]);
        free(tab[i][0]);
        free(tab[i]);
      }  
      free(flP);
      free(lP);
      free(kernel);
      free(conv);
      free(plan);
      free(tab);
    }
    else
    {
      log_fatal("NonLimber not implemented");
      exit(1);
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  if (ni < 0 || ni > tomo.shear_Nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }
  const double lntheta = log(theta);    
  if (lntheta < lnthetamin || lntheta > lnthetamax)
  {
    const double theta = exp(lntheta);
    const double theta_min = exp(lnthetamin);
    const double theta_max = exp(lnthetamax);
    log_fatal("theta = %e outside look-up table range [%e, %e]", theta, theta_min, theta_max);
    exit(1);
  }
  const int q = ni;
  if (q < 0 || q > NSIZE - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  }
  return interpol(table[q], ntheta, lnthetamin, lnthetamax, dlntheta, lntheta, 1, 1);
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Limber Approximation (Angular Power Spectrum)
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// SS ANGULAR CORRELATION FUNCTION - TATT
// -----------------------------------------------------------------------------

// NLA/TA amplitude C1, nz argument only need if per-bin amplitude
double C1_TA(double a, double nz, double growfac_a)
{
  // per-bin IA parameters
  if (like.IA == 3 || like.IA == 5)
  {
    return -nuisance.A_z[(int)nz]*cosmology.Omega_m*nuisance.c1rhocrit_ia/ growfac_a;
  }
  // power law evolution
  return -cosmology.Omega_m * nuisance.c1rhocrit_ia /
         growfac_a * nuisance.A_ia *
         pow(1. / (a * nuisance.oneplusz0_ia), nuisance.eta_ia);
}

// TA source bias parameter, nz argument only need if per-bin amplitude
double b_TA(double a __attribute__((unused)), double nz)
{
  // per-bin IA parameters
  if (like.IA == 5) {
    return nuisance.b_ta_z[(int)nz];
  }
  // power law evolution
  return nuisance.b_ta_z[0];
}

// TT amplitude C2, nz argument only need if per-bin amplitude
double C2_TT(double a, double nz, double growfac_a)
{
  // per-bin IA parameters
  if (like.IA == 5)
  {
    return 5. * nuisance.A2_z[(int)nz] * cosmology.Omega_m *
           nuisance.c1rhocrit_ia * pow(1.0/growfac_a, 2.0);
  }
  // power law evolution
  return 5. * nuisance.A2_ia * cosmology.Omega_m * nuisance.c1rhocrit_ia *
         (1.0 /(growfac_a*growfac_a)) *
         pow(1. / (a * nuisance.oneplusz0_ia), nuisance.eta_ia_tt);
}

double int_for_C_ss_tomo_TATT_EE_limber(double a, void* params)
{
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;
  const int n1 = (int) ar[0]; // first source bin
  const int n2 = (int) ar[1]; // second source bin
  if (n1 < 0 || n1 > tomo.shear_Nbin - 1 || n2 < 0 || n2 > tomo.shear_Nbin - 1)
  {
    log_fatal("error in selecting bin number (ni, nj) = [%d,%d]", n1, n2);
    exit(1);
  }
  const double ell = ar[2] + 0.5;

  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double fK = f_K(chidchi.chi);
  const double k = ell / fK;
  const double PK = Pdelta(k, a);

  const double ws1 = W_source(a, n1, hoverh0); // radial n_z weight for first source bin
  const double ws2 = W_source(a, n2, hoverh0); // radial n_z weight for second source bin
  const double wk1 = W_kappa(a, fK, n1); // radial lens efficiency for first source bin
  const double wk2 = W_kappa(a, fK, n2); // radial lens efficiency for second source bin

  const double C1 = C1_TA(a, n1, growfac_a);   // IA parameters for first source bin
  const double b_ta = b_TA(a, n1);             // IA parameters for first source bin
  const double C2 = C2_TT(a, n1, growfac_a);   // IA parameters for first source bin

  const double C1_2 = C1_TA(a, n2, growfac_a);  // IA parameters for second source bin
  const double b_ta_2 = b_TA(a, n2);            // IA parameters for second source bin
  const double C2_2 = C2_TT(a, n2, growfac_a);  // IA parameters for second source bin
  
  const double tmp1 = TATT_II_EE(k, a, C1, C2, b_ta, C1_2, C2_2, b_ta_2, growfac_a, PK);
  const double tmp2 = TATT_GI_E(k, a, C1, C2, b_ta, growfac_a, PK);
  const double tmp3 = TATT_GI_E(k, a, C1_2, C2_2, b_ta_2, growfac_a, PK);
  
  const double res = wk1 * wk2 * PK + ws1 * ws2 * tmp1 + ws1 * wk2 * tmp2 + ws2 * wk1 * tmp3;
  return res * chidchi.dchida / (fK * fK);
}

double int_for_C_ss_tomo_TATT_BB_limber(double a, void* params)
{
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;
  const int n1 = (int) ar[0]; // first source bin 
  const int n2 = (int) ar[1]; // second source bin
  if (n1 < 0 || n1 > tomo.shear_Nbin - 1 || n2 < 0 || n2 > tomo.shear_Nbin - 1)
  {
    log_fatal("error in selecting bin number (ni, nj) = [%d,%d]", n1, n2);
    exit(1);
  }
  const double ell = ar[2] + 0.5;

  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double fK = f_K(chidchi.chi);
  const double k = ell / fK;

  const double ws_1 = W_source(a, n1, hoverh0);  // radial n_z weight for first source bin 
  const double C1_1 = C1_TA(a, n1, growfac_a);   // IA parameters for first source bin
  const double b_ta_1 = b_TA(a, n1);             // IA parameters for first source bin
  const double C2_1 = C2_TT(a, n1, growfac_a);   // IA parameters for first source bin

  const double ws_2 = W_source(a, n2, hoverh0); // radial n_z weight for second source bin
  const double C1_2 = C1_TA(a, n2, growfac_a);  // IA parameters for second source bin
  const double b_ta_2 = b_TA(a, n2);            // IA parameters for second source bin
  const double C2_2 = C2_TT(a, n2, growfac_a);  // IA parameters for second source bin

  const double tmp1 = TATT_II_BB(k, a, C1_1, C2_1, b_ta_1, C1_2, C2_2, b_ta_2, growfac_a);
  return (ws_1 * ws_2 * tmp1) * chidchi.dchida / (fK * fK);
}

double C_ss_tomo_TATT_EE_limber_nointerp(double l, int ni, int nj, const int init_static_vars_only)
{
  if(ni < -1 || ni > tomo.shear_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }
  double ar[3] = {(double) ni, (double) nj, l};
  const double amin = fmax(amin_source(ni), amin_source(nj));
  const double amax = fmin(amax_source(ni), amax_source(nj));

  return (init_static_vars_only == 1) ? int_for_C_ss_tomo_TATT_EE_limber(amin, (void*) ar) : 
    like.high_def_integration == 2 ?
      int_gsl_integrate_high_precision(int_for_C_ss_tomo_TATT_EE_limber, (void*) ar, 
        amin, amax, NULL, GSL_WORKSPACE_SIZE) :
    like.high_def_integration == 1 ?
      int_gsl_integrate_medium_precision(int_for_C_ss_tomo_TATT_EE_limber, (void*) ar, 
        amin, amax, NULL, GSL_WORKSPACE_SIZE) :
      int_gsl_integrate_low_precision(int_for_C_ss_tomo_TATT_EE_limber, (void*) ar, 
        amin, amax, NULL, GSL_WORKSPACE_SIZE);
}

double C_ss_tomo_TATT_BB_limber_nointerp(double l, int ni, int nj, const int init_static_vars_only)
{
  if(ni < -1 || ni > tomo.shear_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  double ar[3] = {(double) ni, (double) nj, l};
  const double amin = fmax(amin_source(ni), amin_source(nj));
  const double amax = fmin(amax_source(ni), amax_source(nj));
  if (!(amin>0) || !(amin<1) || !(amax>0) || !(amax<1)) 
  {
    log_fatal("0 < amin/amax < 1 not true");
    exit(1);
  }

  return (init_static_vars_only == 1) ? int_for_C_ss_tomo_TATT_BB_limber(amin, (void*) ar) : 
    like.high_def_integration == 2 ?
    int_gsl_integrate_high_precision(int_for_C_ss_tomo_TATT_BB_limber, (void*) ar,
      amin, amax, NULL, GSL_WORKSPACE_SIZE) :
    like.high_def_integration == 1 ?
    int_gsl_integrate_medium_precision(int_for_C_ss_tomo_TATT_BB_limber, (void*) ar,
      amin, amax, NULL, GSL_WORKSPACE_SIZE) :
    int_gsl_integrate_low_precision(int_for_C_ss_tomo_TATT_BB_limber, (void*) ar,
      amin, amax, NULL, GSL_WORKSPACE_SIZE);
}

double C_ss_tomo_TATT_EE_limber(const double l, const int ni, const int nj)
{
  static cosmopara C;
  static nuisancepara N;
  static double** table;
  static double* sig;
  static int osc[100];

  const int NSIZE = tomo.shear_Npowerspectra;
  const int nell = Ntable.N_ell_TATT;
  const double lnlmin = log(fmax(limits.LMIN_tab - 1., 1.0));
  const double lnlmax = log(fmax(limits.LMAX, limits.LMAX_hankel) + 1);
  const double dlnl = (lnlmax - lnlmin) / ((double) nell - 1.);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for (int i=0; i<NSIZE; i++) 
    {
      table[i] = (double*) malloc(sizeof(double)*nell);
    }
    sig = (double*) malloc(sizeof(double)*NSIZE);
  }
  if (recompute_shear(C, N))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      const int k = 0;
      const int Z1NZ = Z1(k);
      const int Z2NZ = Z2(k);
      const double lnl = lnlmin;
      double init = C_ss_tomo_TATT_EE_limber_nointerp(exp(lnl), Z1NZ, Z2NZ, 1);
    }
    #pragma GCC diagnostic pop
    #pragma omp parallel for collapse(2)
    for (int k=0; k<NSIZE; k++)
    {
      for (int i=0; i<nell; i++)
      {
        const int Z1NZ = Z1(k);
        const int Z2NZ = Z2(k);
        const double lnl = lnlmin + i*dlnl;
        table[k][i] = C_ss_tomo_TATT_EE_limber_nointerp(exp(lnl), Z1NZ, Z2NZ, 0);
      }
    }
    #pragma omp parallel for
    for (int k=0; k<NSIZE; k++)
    {
      const int Z1NZ = Z1(k);
      const int Z2NZ = Z2(k);
      sig[k] = 1.;
      osc[k] = 0;
      if (C_ss_tomo_TATT_EE_limber_nointerp(500., Z1NZ, Z2NZ, 0) < 0)
      {
        sig[k] = -1.;
      }
      for (int i=0; i<nell; i++)
      {
        if (table[k][i] * sig[k] < 0.)
        {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0)
      {
        for (int i = 0; i<nell; i++)
        {
          table[k][i] = log(sig[k] * table[k][i]);
        }
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  if (ni < 0 || ni > tomo.shear_Nbin - 1 || nj < 0 || nj > tomo.shear_Nbin - 1)
  {
    log_fatal("error in selecting bin number (ni, nj) = [%d,%d]", ni, nj);
    exit(1);
  }
  const double lnl = log(l);
  if (lnl < lnlmin)
  {
    log_warn("l = %e < l_min = %e. Extrapolation adopted", l, exp(lnlmin));
  }
  if (lnl > lnlmax)
  {
    log_warn("l = %e > l_max = %e. Extrapolation adopted", l, exp(lnlmax));
  }
  const int q = N_shear(ni, nj);
  if (q < 0 || q > NSIZE - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  }
  double f1;
  if (osc[q] == 0)
  {
    f1 = sig[q] * exp(interpol(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 1, 1));
  }
  else
  {
    f1 = interpol(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 1, 1);
  }
  return isnan(f1) ? 0.0 : f1;
}

double C_ss_tomo_TATT_BB_limber(const double l, const int ni, const int nj)
{
  static cosmopara C;
  static nuisancepara N;
  static double** table;
  static double* sig;
  static int osc[100];

  const int nell = Ntable.N_ell_TATT;
  const int NSIZE = tomo.shear_Npowerspectra;
  const double lnlmin = log(fmax(limits.LMIN_tab, 1.0));
  const double lnlmax = log(fmax(limits.LMAX, limits.LMAX_hankel) + 1);
  const double dlnl = (lnlmax - lnlmin) / ((double) nell - 1.0);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for (int i=0; i<NSIZE; i++) 
    {
      table[i] = (double*) malloc(sizeof(double)*nell);
    }
    sig = (double*) malloc(sizeof(double)*NSIZE);
  }
  if (recompute_shear(C, N))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      const int k = 0;
      const int Z1NZ = Z1(k);
      const int Z2NZ = Z2(k);
      const double lnl = lnlmin;
      double init = C_ss_tomo_TATT_BB_limber_nointerp(exp(lnl), Z1NZ, Z2NZ, 1);
    }
    #pragma GCC diagnostic pop
    #pragma omp parallel for collapse(2)
    for (int k=0; k<NSIZE; k++)
    {
      for (int i=0; i<nell; i++)
      {
        const int Z1NZ = Z1(k);
        const int Z2NZ = Z2(k);
        const double lnl = lnlmin + i*dlnl;
        table[k][i] = C_ss_tomo_TATT_BB_limber_nointerp(exp(lnl), Z1NZ, Z2NZ, 0);
      }
    }
    #pragma omp parallel for
    for (int k=0; k<NSIZE; k++)
    {
      const int Z1NZ = Z1(k);
      const int Z2NZ = Z2(k);
      sig[k] = 1.;
      osc[k] = 0;
      if (C_ss_tomo_TATT_BB_limber_nointerp(500., Z1NZ, Z2NZ, 0) < 0)
      {
        sig[k] = -1.;
      }
      for (int i=0; i<nell; i++)
      {
        if (table[k][i] * sig[k] < 0.)
        {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0)
      {
        #pragma omp parallel for
        for (int i=0; i<nell; i++)
        {
          table[k][i] = log(sig[k] * table[k][i]);
        }
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  if (ni < 0 || ni > tomo.shear_Nbin - 1 || nj < 0 || nj > tomo.shear_Nbin - 1)
  {
    log_fatal("error in selecting bin number (ni, nj) = [%d,%d]", ni, nj);
    exit(1);
  }
  const double lnl = log(l);
  if (lnl < lnlmin)
  {
    log_warn("l = %e < l_min = %e. Extrapolation adopted", l, exp(lnlmin));
  }
  if (lnl > lnlmax)
  {
    log_warn("l = %e > l_max = %e. Extrapolation adopted", l, exp(lnlmax));
  }
  const int q = N_shear(ni, nj);
  if (q < 0 || q > NSIZE - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  }
  double f1;
  if (osc[q] == 0)
  {
    f1 = sig[q] * exp(interpol(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 1, 1));
  }
  else
  {
    f1 = interpol(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 1, 1);
  }
  return isnan(f1) ? 0.0 : f1;
}

// -----------------------------------------------------------------------------
// SS ANGULAR CORRELATION FUNCTION - NON TATT
// -----------------------------------------------------------------------------

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
  if (n1 < 0 || n1 > tomo.shear_Nbin - 1 || n2 < 0 || n2 > tomo.shear_Nbin - 1)
  {
    log_fatal("error in selecting bin number (ni, nj) = [%d,%d]", n1, n2);
    exit(1);
  }
  const double l = ar[2];
  const int use_linear_ps = (int) ar[3];

  const double ell = l + 0.5;
  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
 
  const double ws1 = W_source(a, n1, hoverh0);
  const double ws2 = W_source(a, n2, hoverh0);
  const double wk1 = W_kappa(a, fK, n1);
  const double wk2 = W_kappa(a, fK, n2);
  const double PK = (use_linear_ps == 1) ? p_lin(k,a) : Pdelta(k,a);
  
  const double ell4 = ell*ell*ell*ell; // correction (1812.05995 eqs 74-79)
  const double ell_prefactor = l*(l - 1.)*(l + 1.)*(l + 2.)/ell4; 

  double res = 0.0;
  switch(like.IA)
  {
    case 0:
    {
      res = wk1*wk2;

      break;
    }
    case 1:
    {
      const double norm = A_IA_Joachimi(a)*cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a;
      res = ws1*ws2*norm*norm - (ws1*wk2+ws2*wk1)*norm + wk1*wk2;

      break;
    }
    case 3:
    {
      const double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a;
      res = ws1*ws2*norm*norm - (ws1*wk2+ws2*wk1)*norm + wk1*wk2;

      break;
    }
    case 4:
    {
      const double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a*
        nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia), nuisance.eta_ia);
      res = ws1*ws2*norm*norm - (ws1*wk2+ws2*wk1)*norm + wk1*wk2;

      break;
    }
    default:
    {
      log_fatal("like.IA = %d not supported", like.IA);
      exit(1);
    }
  }
  return res*PK*(chidchi.dchida/(fK*fK))*ell_prefactor;
}

double C_ss_tomo_limber_nointerp(double l, int ni, int nj, int use_linear_ps, 
const int init_static_vars_only)
{
  if(ni < -1 || ni > tomo.shear_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  double ar[4] = {(double) ni, (double) nj, l, (double) use_linear_ps};
  const double amin = fmax(amin_source(ni), amin_source(nj));
  const double amax = fmin(amax_source(ni), amax_source(nj));;
  if (!(amin>0) || !(amin<1) || !(amax>0) || !(amax<1)) 
  {
    log_fatal("0 < amin/amax < 1 not true");
    exit(1);
  }

  return (init_static_vars_only == 1) ? int_for_C_ss_tomo_limber(amin, (void*) ar) :
    like.high_def_integration == 2 ?
    int_gsl_integrate_high_precision(int_for_C_ss_tomo_limber, (void*) ar, amin, amax, NULL, 
      GSL_WORKSPACE_SIZE) :
    like.high_def_integration == 1 ?
    int_gsl_integrate_medium_precision(int_for_C_ss_tomo_limber, (void*) ar, amin, amax, NULL, 
      GSL_WORKSPACE_SIZE) :
    int_gsl_integrate_low_precision(int_for_C_ss_tomo_limber, (void*) ar, amin, amax, NULL, 
      GSL_WORKSPACE_SIZE);       
}

double C_ss_tomo_limber(double l, int ni, int nj)
{
  if(ni < -1 || ni > tomo.shear_Nbin -1 || nj < -1 || nj > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", ni, nj);
    exit(1);
  }

  static cosmopara C;
  static nuisancepara N;
  static double **table;  

  const int nell = Ntable.N_ell;
  const int NSIZE = tomo.shear_Npowerspectra;
  const double lnlmin = log(fmax(limits.LMIN_tab - 1., 1.0));
  const double lnlmax = log(fmax(limits.LMAX, limits.LMAX_hankel) + 1);
  const double dlnl = (lnlmax - lnlmin)/(nell - 1.);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for (int i=0; i<NSIZE; i++) 
    {
      table[i] = (double*) malloc(sizeof(double)*nell);
    }
  }

  if (recompute_shear(C, N))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      const int k = 0;
      double init = C_ss_tomo_limber_nointerp(exp(lnlmin), Z1(k), Z2(k), use_linear_ps_limber, 1);
    }
    #pragma GCC diagnostic pop
    #pragma omp parallel for collapse(2)
    for (int k=0; k<NSIZE; k++)
    {
      for (int i=0; i<nell; i++)
      {
        const int Z1NZ = Z1(k);
        const int Z2NZ = Z2(k);
        const double lnl = lnlmin + i*dlnl;
        const double l = exp(lnl);
        table[k][i]= log(C_ss_tomo_limber_nointerp(l, Z1NZ, Z2NZ, use_linear_ps_limber, 0));
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }

  const double lnl = log(l);
  if (lnl < lnlmin || lnl > lnlmax)
  {
    log_fatal("l = %e outside look-up table range [%e,%e]", l, exp(lnlmin), exp(lnlmax));
    exit(1);
  }

  const int q = N_shear(ni, nj);
  if(q > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }
  
  const double f1 = exp(interpol_fitslope(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 1.));
  if (isnan(f1)) 
  {
    return 0.0;
  }
  return f1;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// GS ANGULAR CORRELATION FUNCTION - TATT
// -----------------------------------------------------------------------------

double int_for_C_gs_tomo_limber_TATT(double a, void* params)
{
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;
  const int nl = (int) ar[0];
  const int ns = (int) ar[1];
  if (nl < 0 || nl > tomo.clustering_Nbin - 1 || ns < 0 || ns > tomo.shear_Nbin - 1)
  {
    log_fatal("error in selecting bin number (nl, ns) = [%d,%d]", nl, ns);
    exit(1);
  }
  const double l = ar[2];

  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double g4 = growfac_a*growfac_a*growfac_a*growfac_a;
  const double ell = l + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell / fK;
  const double z = 1. / a - 1.;

  const double WGAL = W_gal(a, nl, hoverh0);
  const double WMAG = W_mag(a, fK, nl); 
  const double WS = W_source(a, ns, hoverh0);
  const double WK = W_kappa(a, fK, ns);
  
  const double b1 = gbias.b1_function(z, nl);
  const double bmag = gbias.b_mag[nl];
  const double b2 = gbias.b2[nl];
  const double bs2 = gbias.bs2[nl];
  const double b3 = b3nl_from_b1(b1);

  const double PK = Pdelta(k, a);
  const double tmp1 = 0.5*g4*(b2 * PT_d1d2(k) + bs2 * PT_d1s2(k) + b3 * PT_d1d3(k));
  const double PK1loop = (WGAL*b2 != 0) ? tmp1 : 0.0; 

  const double C1 = C1_TA(a, ns, growfac_a);
  const double b_ta = b_TA(a, ns);
  const double C2 = C2_TT(a, ns, growfac_a);
  const double GIE = TATT_GI_E(k, a, C1, C2, b_ta, growfac_a, PK);

  const double res = WK*(WMAG*bmag*PK + WGAL*(b1*PK + PK1loop)) + (WGAL*b1 + WMAG*bmag)*WS*GIE;
  return res*chidchi.dchida/(fK*fK);
}

// -----------------------------------------------------------------------------
// GS ANGULAR CORRELATION FUNCTION - NON TATT
// -----------------------------------------------------------------------------

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
  if (nl < 0 || nl > tomo.clustering_Nbin - 1 || ns < 0 || ns > tomo.shear_Nbin - 1)
  {
    log_fatal("error in selecting bin number (nl, ns) = [%d,%d]", nl, ns);
    exit(1);
  }
  const double l = ar[2];
  const int use_linear_ps = (int) ar[3];

  const double ell = l + 0.5;
  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double z = 1.0/a - 1.0;

  const double b1 = gbias.b1_function(z, nl);
  const double bmag = gbias.b_mag[nl];

  const double WS = W_source(a, ns, hoverh0);
  const double WK = W_kappa(a, fK, ns);
  const double WGAL = W_gal(a, nl, hoverh0);
  const double WMAG = W_mag(a, fK, nl);

  const double ell_prefactor = l*(l + 1.)/(ell*ell);          // correction (1812.05995 eqs 74-79)
  const double tmp = (l - 1.)*l*(l + 1.)*(l + 2.);            // correction (1812.05995 eqs 74-79)
  const double ell_prefactor2 = (tmp > 0) ? sqrt(tmp)/(ell*ell) : 0.0;

  double res = 0.0;
  switch(like.IA)
  {
    case 0:
    {
      res = WK;

      break;
    }
    case 1:
    {
      const double norm = A_IA_Joachimi(a)*cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a;

      res = (WK - WS*norm);

      break;
    }
    case 3:
    {
      const double norm = nuisance.A_z[ns]*cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a;

      res = (WK - WS*norm);

      break;
    }
    case 4:
    {
      const double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a*
        nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);

      res = (WK - WS*norm);

      break;
    }
    default:
    {
      log_fatal("like.IA = %d not supported", like.IA);
      exit(1);
    }
  }

  if (include_HOD_GX == 1)
  {
    if (include_RSD_GS == 1)
    {
      log_fatal("RSD not implemented with (HOD = TRUE)");
      exit(1);
    }
    else
    {
      res *= WGAL;
    }

    const double PK = p_gm(k, a, nl, use_linear_ps);
    res *= PK;
  }
  else
  {
    if (include_RSD_GS == 1)
    {
      const double chi_0 = f_K(ell/k);
      const double chi_1 = f_K((ell+1.)/k);
      const double a_0 = a_chi(chi_0);
      const double a_1 = a_chi(chi_1);
      const double WRSD = W_RSD(ell, a_0, a_1, nl);

      res *= (WGAL*b1 + WMAG*ell_prefactor*bmag) + WRSD;
    }
    else
    {
      res *= (WGAL*b1 + WMAG*ell_prefactor*bmag);
    }

    const double PK = use_linear_ps == 1 ? p_lin(k,a) : Pdelta(k,a);
    res *= PK;
  }
  
  return res*(chidchi.dchida/(fK*fK))*ell_prefactor2;
}

double int_for_C_gs_tomo_limber_withb2(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;
  const int nl = (int) ar[0];
  const int ns = (int) ar[1];
  if (nl < 0 || nl > tomo.clustering_Nbin - 1 || ns < 0 || ns > tomo.shear_Nbin - 1)
  {
    log_fatal("error in selecting bin number (nl, ns) = [%d,%d]", nl, ns);
    exit(1);
  }
  const double l = ar[2];
  //const int use_linear_ps = (int) ar[3];

  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double g4 = growfac_a*growfac_a*growfac_a*growfac_a;
  const double ell = l + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double z = 1.0/a - 1.0;

  const double b1 = gbias.b1_function(z, nl);
  const double b2 = gbias.b2[nl];
  const double bs2 = gbias.bs2[nl];
  const double b3 = b3nl_from_b1(b1);
  const double bmag = gbias.b_mag[nl];

  const double WS = W_source(a, ns, hoverh0);
  const double WK = W_kappa(a, fK, ns);
  const double WGAL = W_gal(a, nl, hoverh0);
  const double WMAG = W_mag(a, fK, nl);

  const double ell_prefactor = l*(l + 1.)/(ell*ell);       // correction (1812.05995 eqs 74-79)
  const double tmp = (l - 1.)*l*(l + 1.)*(l + 2.);         // correction (1812.05995 eqs 74-79)
  const double ell_prefactor2 = (tmp > 0) ? sqrt(tmp)/(ell*ell) : 0.0;
  
  double linear_part;
  double non_linear_part;
  switch(like.IA)
  {
    case 0:
    {
      linear_part = WK;
      non_linear_part = WK;

      break;
    }
    case 1:
    {
      const double norm = A_IA_Joachimi(a)*cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a;

      linear_part = (WK - WS*norm);
      non_linear_part = (WK - WS*norm);

      break;
    }
    case 3:
    {
      const double norm = nuisance.A_z[(int)ar[1]]*cosmology.Omega_m*
        nuisance.c1rhocrit_ia/growfac_a;

      linear_part = (WK - WS*norm);
      non_linear_part = (WK - WS*norm);

      break;
    }
    case 4:
    {
      const double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a*
        nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);

      linear_part = (WK - WS*norm);
      non_linear_part = (WK - WS*norm);

      break;
    }
    default:
    {
      log_fatal("like.IA = %d not supported", like.IA);
      exit(1);
    }
  }

  if(include_RSD_GS == 1)
  {
    static double chi_a_min = 0;
    if (chi_a_min == 0)
    {
      chi_a_min = chi(limits.a_min);
    }
    const double chi_0 = f_K(ell/k);
    const double chi_1 = f_K((ell + 1.0)/k);
    if (chi_1 > chi_a_min)
    {
      return 0;
    }
    const double a_0 = a_chi(chi_0);
    const double a_1 = a_chi(chi_1);
    const double WRSD = W_RSD(ell, a_0, a_1, nl);

    linear_part *= (WGAL*b1 + WMAG*ell_prefactor*bmag) + WRSD;
    non_linear_part *= WGAL;
  }
  else
  {
    linear_part *= (WGAL*b1 + WMAG*ell_prefactor*bmag);
    non_linear_part *= WGAL;
  }

  linear_part *= Pdelta(k,a);
  non_linear_part *= g4*0.5*(b2*PT_d1d2(k) + bs2*PT_d1s2(k) + b3*PT_d1d3(k));

  return (linear_part + non_linear_part)*(chidchi.dchida/(fK*fK))*ell_prefactor2;
}

double C_gs_tomo_limber_nointerp(double l, int nl, int ns, int use_linear_ps,
const int init_static_vars_only)
{
  if(nl < -1 || nl > tomo.clustering_Nbin -1 || ns < -1 || ns > tomo.shear_Nbin -1)
  {
    log_fatal("invalid bin input (ni, nj) = (%d, %d)", nl, ns);
    exit(1);
  }

  double ar[4] = {(double) nl, (double) ns, l, use_linear_ps};
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
  if (like.IA == 0 || like.IA == 1 || like.IA == 3 || like.IA == 4)
  {
    if (has_b2_galaxies() && use_linear_ps == 0)
    {
      res = (init_static_vars_only == 1) ? int_for_C_gs_tomo_limber_withb2(amin, (void*) ar) :
        like.high_def_integration == 2 ?
        int_gsl_integrate_high_precision(int_for_C_gs_tomo_limber_withb2, (void*) ar, amin, amax, 
          NULL, GSL_WORKSPACE_SIZE) :
        like.high_def_integration == 1 ?
        int_gsl_integrate_medium_precision(int_for_C_gs_tomo_limber_withb2, (void*) ar, amin, amax, 
          NULL, GSL_WORKSPACE_SIZE) :
        int_gsl_integrate_low_precision(int_for_C_gs_tomo_limber_withb2, (void*) ar, amin, amax, 
          NULL, GSL_WORKSPACE_SIZE);
    }
    else 
    {
      res = (init_static_vars_only == 1) ? int_for_C_gs_tomo_limber(amin, (void*) ar) :
        like.high_def_integration == 2 ?
        int_gsl_integrate_high_precision(int_for_C_gs_tomo_limber, (void*) ar, amin, amax, NULL, 
          GSL_WORKSPACE_SIZE) :
        like.high_def_integration == 1 ?
        int_gsl_integrate_medium_precision(int_for_C_gs_tomo_limber, (void*) ar, amin, amax, NULL, 
          GSL_WORKSPACE_SIZE) :
        int_gsl_integrate_low_precision(int_for_C_gs_tomo_limber, (void*) ar, amin, amax, NULL, 
          GSL_WORKSPACE_SIZE);
    }
  }
  else if (like.IA == 5 || like.IA == 6)
  {
    if (use_linear_ps == 1)
    {
      log_fatal("use linear power spectrum option not implemented with TATT");
      exit(1);
      return 0;
    }
    else
    {
      res = (init_static_vars_only == 1) ? int_for_C_gs_tomo_limber_TATT(amin, (void*) ar) :
        like.high_def_integration == 2 ?
        int_gsl_integrate_high_precision(int_for_C_gs_tomo_limber_TATT, (void*) ar, amin, amax, 
          NULL, GSL_WORKSPACE_SIZE) :
        like.high_def_integration == 1 ?
        int_gsl_integrate_medium_precision(int_for_C_gs_tomo_limber_TATT, (void*) ar, amin, amax, 
          NULL, GSL_WORKSPACE_SIZE) :
        int_gsl_integrate_low_precision(int_for_C_gs_tomo_limber_TATT, (void*) ar, amin, amax, 
          NULL, GSL_WORKSPACE_SIZE);
    }
  }
  else
  {
    log_fatal("like.IA = %d not supported", like.IA);
    exit(1);
    return 0;
  }
  return res;
}

double C_gs_tomo_limber(double l, int ni, int nj)
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double** table;
  static double* sig;
  static int osc[100];

  const int NSIZE = tomo.ggl_Npowerspectra;
  const int nell = (like.IA == 5 || like.IA == 6) ? Ntable.N_ell_TATT :  Ntable.N_ell;
  const double lnlmin = log(fmax(limits.LMIN_tab, 1.0));
  const double lnlmax = log(fmax(limits.LMAX, limits.LMAX_hankel) + 1);
  const double dlnl = (lnlmax - lnlmin) / ((double) nell - 1.0);

  if (table == 0) 
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for (int i=0; i<NSIZE; i++) 
    {
      table[i] = (double*) malloc(sizeof(double)*nell);
    }
    sig = (double*) malloc(sizeof(double)*NSIZE);
  }
  if (recompute_gs(C, G, N))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init = C_gs_tomo_limber_nointerp(exp(lnlmin), ZL(0), ZS(0), use_linear_ps_limber, 1);
    }
    #pragma GCC diagnostic pop

    if (like.IA == 5 || like.IA == 6) // TATT MODELING
    { 
      #pragma omp parallel for collapse(2)
      for (int k=0; k<NSIZE; k++)
      {
        for (int i=0; i<nell; i++)
        {
          const int ZLNZ = ZL(k);
          const int ZSNZ = ZS(k);
          const double lnl = lnlmin + i*dlnl;
          table[k][i] = C_gs_tomo_limber_nointerp(exp(lnl), ZLNZ, ZSNZ, use_linear_ps_limber, 0);
        }
      }
    }
    else
    {
      #pragma omp parallel for collapse(2)
      for (int k=0; k<NSIZE; k++)
      {
        for (int i=0; i<nell; i++)
        {
          const int ZLNZ = ZL(k);
          const int ZSNZ = ZS(k);
          const double lnl = lnlmin + i*dlnl;
          const double ll = exp(lnl);
          table[k][i] = log(C_gs_tomo_limber_nointerp(ll, ZLNZ, ZSNZ, use_linear_ps_limber, 0));
        }
      }
    }
    #pragma omp parallel for
    for (int k=0; k<NSIZE; k++)
    {
      const int ZLNZ = ZL(k);
      const int ZSNZ = ZS(k);
      sig[k] = 1.;
      osc[k] = 0;
      if (C_gs_tomo_limber_nointerp(500., ZLNZ, ZSNZ, use_linear_ps_limber, 0) < 0)
      {
        sig[k] = -1.;
      }
      for (int i=0; i<nell; i++)
      {
        if (table[k][i] * sig[k] < 0.)
        {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0)
      {
        for (int i=0; i<nell; i++)
        {
          table[k][i] = log(sig[k] * table[k][i]);
        }
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
    update_galpara(&G);
  }
    
  const int q = N_ggl(ni, nj);
  if (q < 0 || q > NSIZE - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  }

  const double lnl = log(l);
  if (lnl < lnlmin)
  {
    log_warn("l = %e < l_min = %e. Extrapolation adopted", l, exp(lnlmin));
  }
  if (lnl > lnlmax)
  {
    log_warn("l = %e > l_max = %e. Extrapolation adopted", l, exp(lnlmax));
  }

  double f1 = 0.;
  if (test_zoverlap(ni, nj) && osc[q] == 0)
  {
    f1 = sig[q] * exp(interpol(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 0, 0));
  }
  else if (test_zoverlap(ni, nj) && osc[q] == 1)
  {
    f1 = interpol(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 0, 0);
  }
  if (isnan(f1))
  {
    f1 = 0;
  }
  return f1;
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

  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double ell = l + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell / fK;
  const double z = 1.0/a - 1.0;

  const double b1i = gbias.b1_function(z, ni);
  const double bmagi = gbias.b_mag[ni];
  const double WGALi = W_gal(a, ni, hoverh0);
  const double WMAGi = W_mag(a, fK, ni);
  
  const double b1j = gbias.b1_function(z, nj);
  const double bmagj = gbias.b_mag[nj];
  const double WGALj = W_gal(a, nj, hoverh0);
  const double WMAGj = W_mag(a, fK, nj);

  const double ell_prefactor = l*(l + 1.)/(ell*ell); // prefactor correction (1812.05995 eqs 74-79)
  
  double res = 1.0;
  if (include_HOD_GX == 1)
  {
    if (include_RSD_GG == 1)
    {
      log_fatal("RSD not implemented with (HOD = TRUE)");
      exit(1);
    }
    else
    { 
      res *= WGALi*WGALj;
    }
    const double PK = p_gg(k, a, ni, nj, use_linear_ps);
    res *= PK;
  }
  else
  {
    if(include_RSD_GG == 1)
    {
      static double chi_a_min = 0;
      if (chi_a_min == 0)
      {
        chi_a_min = chi(limits.a_min);
      }
      const double chi_0 = f_K(ell/k);
      const double chi_1 = f_K((ell + 1.0)/k);
      if (chi_1 > chi_a_min)
      {
        return 0;
      }
      const double a_0 = a_chi(chi_0);
      const double a_1 = a_chi(chi_1);
      
      const double WRSDi =  W_RSD(ell, a_0, a_1, ni);
      const double WRSDj =  W_RSD(ell, a_0, a_1, nj);

      res *= (WGALi*b1i + WMAGi*ell_prefactor*bmagi + WRSDi);
      res *= (WGALj*b1j + WMAGj*ell_prefactor*bmagj + WRSDj);
    }
    else
    {
      res *= (WGALi*b1i + WMAGi*ell_prefactor*bmagi);
      res *= (WGALj*b1j + WMAGj*ell_prefactor*bmagj);
    }

    const double PK = (use_linear_ps == 1) ? p_lin(k,a) : Pdelta(k,a);
    res *= PK;
  }

  return res*chidchi.dchida/(fK*fK);
}

double int_for_C_gg_tomo_limber_withb2(double a, void* params)
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
  //const int use_linear_ps = (int) ar[3];

  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double g4 = growfac_a*growfac_a*growfac_a*growfac_a;
  const double ell = l + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double z = 1./a - 1.;

  const double b1i = gbias.b1_function(z, ni);
  const double bmagi = gbias.b_mag[ni];
  const double WGALi = W_gal(a, ni, hoverh0);
  const double WMAGi = W_mag(a, fK, ni);
  
  const double b1j = gbias.b1_function(z, nj);
  const double bmagj = gbias.b_mag[nj];
  const double WGALj = W_gal(a, nj, hoverh0);
  const double WMAGj = W_mag(a, fK, nj);
  
  // only i = j has been implemented for higher order bias
  const double b2 = gbias.b2[ni];
  const double bs2 = gbias.bs2[ni];
  const double b3 = b3nl_from_b1(b1i);

  const double s4 = 0.; // PT_sigma4(k);

  const double ell_prefactor = l*(l + 1.)/(ell*ell); // prefactor correction (1812.05995 eqs 74-79)

  double linear_part = 1.0;
  double non_linear_part = 1.0;
  if(include_RSD_GG == 1)
  {
    static double chi_a_min = 0;
    if (chi_a_min == 0)
    {
      chi_a_min = chi(limits.a_min);
    }
    const double chi_0 = f_K(ell/k);
    const double chi_1 = f_K((ell+1.)/k);
    if (chi_1 > chi_a_min)
    {
      return 0;
    }
    const double a_0 = a_chi(chi_0);
    const double a_1 = a_chi(chi_1);
    const double WRSDi = W_RSD(ell, a_0, a_1, ni);
    const double WRSDj = W_RSD(ell, a_0, a_1, nj);

    linear_part *= WGALi*b1i + WMAGi*ell_prefactor*bmagi + WRSDi;
    linear_part *= WGALj*b1j + WMAGj*ell_prefactor*bmagj + WRSDj;

    non_linear_part *= WGALi*WGALj;
  }
  else
  {
    linear_part *= WGALi*b1i + WMAGi*ell_prefactor*bmagi;
    linear_part *= WGALj*b1j + WMAGj*ell_prefactor*bmagj;
    
    non_linear_part *= WGALi*WGALj;
  }

  linear_part *= Pdelta(k,a);

  non_linear_part *= g4*(b1i * b2 * PT_d1d2(k) + 0.25 * b2 * b2 * (PT_d2d2(k) - 2. * s4) +
    b1i * bs2 * PT_d1s2(k) + 0.5 * b2 * bs2 * (PT_d2s2(k) - 4. / 3. * s4) +
    0.25 * bs2 * bs2 * (PT_s2s2(k) - 8. / 9. * s4) + b1i * b3 * PT_d1d3(k));

  return (linear_part + non_linear_part)*(chidchi.dchida/(fK*fK));
}

double C_gg_tomo_limber_nointerp(double l, int ni, int nj, int use_linear_ps,
const int init_static_vars_only)
{
  if (ni < 0 || ni > tomo.clustering_Nbin - 1 || nj < 0 || nj > tomo.clustering_Nbin - 1)
  {
    log_fatal("error in selecting bin number (ni, nj) = [%d,%d]", ni, nj);
    exit(1);
  }
  if (ni != nj)
  {
    log_fatal("cross-tomography (ni,nj) = (%d,%d) bins not supported", ni, nj);
    exit(1);
  }
  double ar[4] = {(double) ni, (double) nj, l, (double) use_linear_ps};
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

  double res;
  if (has_b2_galaxies() && use_linear_ps == 0)
  {
    res = (init_static_vars_only == 1) ? int_for_C_gg_tomo_limber_withb2(amin, (void*) ar) :
      like.high_def_integration == 2 ?
      int_gsl_integrate_high_precision(int_for_C_gg_tomo_limber_withb2, (void*) ar, amin, amax, 
        NULL, GSL_WORKSPACE_SIZE) :
      like.high_def_integration == 1 ?
      int_gsl_integrate_medium_precision(int_for_C_gg_tomo_limber_withb2, (void*) ar, amin, amax, 
        NULL, GSL_WORKSPACE_SIZE) :
      int_gsl_integrate_low_precision(int_for_C_gg_tomo_limber_withb2, (void*) ar, amin, amax, NULL, 
        GSL_WORKSPACE_SIZE);
  }
  else
  {
    res = (init_static_vars_only == 1) ? int_for_C_gg_tomo_limber(amin, (void*) ar) :
      like.high_def_integration > 0 ?
      int_gsl_integrate_high_precision(int_for_C_gg_tomo_limber, (void*) ar, amin, amax, NULL,
        GSL_WORKSPACE_SIZE) :
      int_gsl_integrate_medium_precision(int_for_C_gg_tomo_limber, (void*) ar, amin, amax, NULL,
        GSL_WORKSPACE_SIZE);
  }
  return res;
}

double C_gg_tomo_limber(double l, int ni, int nj)
{ // cross redshift bin not supported
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double** table;

  const int nell = Ntable.N_ell;
  const double lnlmin = log(fmax(limits.LMIN_tab, 1.0));
  const double lnlmax = log(fmax(limits.LMAX, limits.LMAX_hankel) + 1);
  const double dlnl = (lnlmax - lnlmin) / ((double) nell - 1.0);
  const int NSIZE = tomo.clustering_Nbin;

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for (int i=0; i<NSIZE; i++) 
    {
      table[i] = (double*) malloc(sizeof(double)*nell);
    }
  }
  if (recompute_gg(C, G, N))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init = C_gg_tomo_limber_nointerp(exp(lnlmin), 0, 0, use_linear_ps_limber, 1);
    }
    #pragma GCC diagnostic pop
    
    #pragma omp parallel for collapse(2)
    for (int k=0; k<NSIZE; k++)  
    {
      for (int p=0; p<nell; p++)
      {
        const int ZCL1 = k; // cross redshift bin not supported so not using ZCL1(k)
        const int ZCL2 = k; // cross redshift bin not supported so not using ZCL2(k)
        const double lnl = lnlmin + p*dlnl;
        const double res = C_gg_tomo_limber_nointerp(exp(lnl), ZCL1, ZCL2, use_linear_ps_limber, 0);
        table[k][p] = (res <= 0) ? -100 : log(res);
      }
    }
    update_galpara(&G);
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  if (ni < 0 || ni > tomo.clustering_Nbin - 1 || nj < 0 || nj > tomo.clustering_Nbin - 1)
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
  {
    log_warn("l = %e < l_min = %e. Extrapolation adopted", l, exp(lnlmin));
  }
  if (lnl > lnlmax)
  {
    log_warn("l = %e > l_max = %e. Extrapolation adopted", l, exp(lnlmax));
  }
  const int q = ni; // cross redshift bin not supported so not using N_CL(ni, nj) instead of ni
  if (q < 0 || q > NSIZE - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  }
  const double f1 = exp(interpol_fitslope(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 1));
  return isnan(f1) ? 0.0 : f1;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_gk_limber(double a, void* params)
{
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;

  const int nl = (int) ar[0];
  if (nl < 0 || nl > tomo.clustering_Nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", nl);
    exit(1);
  }
  const double l = ar[1];
  const int use_linear_ps = (int) ar[2];

  const double ell = l + 0.5;
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double z = 1./a - 1.;

  const double b1 = gbias.b1_function(z, nl);
  const double bmag = gbias.b_mag[nl];

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
    const double PK = p_gm(k, a, nl, use_linear_ps);
    res *= PK;
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
    const double PK = use_linear_ps == 1 ? p_lin(k,a) : Pdelta(k,a);
    res *= PK;
  }

  return (res*chidchi.dchida/(fK*fK))*ell_prefactor;
}

double int_for_C_gk_limber_withb2(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;

  const int nl = (int) ar[0];
  if (nl < 0 || nl > tomo.clustering_Nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", nl);
    exit(1);
  }
  const double l = ar[1];
  //const int use_linear_ps = (int) ar[2];
  
  const double growfac_a = growfac(a);
  const double g4 = growfac_a*growfac_a*growfac_a*growfac_a;
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double ell = l + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK; 
  const double z = 1./a - 1.;
  
  const double b2 = gbias.b2[nl];
  const double bs2 = gbias.bs2[nl];
  const double b1 = gbias.b1_function(z, nl);
  const double bmag = gbias.b_mag[nl];

  const double WK = W_k(a, fK);
  const double WGAL = W_gal(a, nl, hoverh0);
  const double WMAG = W_mag(a, fK, nl);
  
  const double ell_prefactor = l*(l + 1.)/(ell*ell); // prefactor correction (1812.05995 eqs 74-79)

  double linear_part = WK; 
  double non_linear_part = WK; 

  if(include_RSD_GK == 1)
  {
    static double chi_a_min = 0;
    if (chi_a_min == 0)
    {
      chi_a_min = chi(limits.a_min);
    }
    const double chi_0 = f_K(ell/k);
    const double chi_1 = f_K((ell + 1.0)/k);
    if (chi_1 > chi_a_min)
    {
      return 0;
    }
    const double a_0 = a_chi(chi_0);
    const double a_1 = a_chi(chi_1);
    const double WRSD = W_RSD(ell, a_0, a_1, nl);

    linear_part *= (WGAL*b1 + WMAG*ell_prefactor*bmag) + WRSD;
    non_linear_part *= WGAL;
  }
  else
  {
    linear_part *= (WGAL*b1 + WMAG*ell_prefactor*bmag);
    non_linear_part *= WGAL;
  }

  linear_part *= Pdelta(k,a);

  non_linear_part *= g4*(0.5*b2*PT_d1d2(k) + 0.5*bs2*PT_d1s2(k));

  return (linear_part + non_linear_part)*(chidchi.dchida/(fK*fK))*ell_prefactor;
}

double C_gk_tomo_limber_nointerp(double l, int ni, int use_linear_ps, 
const int init_static_vars_only)
{
  if (ni < 0 || ni > tomo.clustering_Nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }
  double ar[3] = {(double) ni, l, (double) use_linear_ps};
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
  if (has_b2_galaxies() && use_linear_ps == 0)
  {
    return (init_static_vars_only == 1) ? int_for_C_gk_limber_withb2(amin, (void*) ar) :
      like.high_def_integration > 0 ?
      int_gsl_integrate_high_precision(int_for_C_gk_limber_withb2, (void*) ar, amin, amax, NULL, 
        GSL_WORKSPACE_SIZE) :
      int_gsl_integrate_medium_precision(int_for_C_gk_limber_withb2, (void*) ar, amin, amax, NULL, 
        GSL_WORKSPACE_SIZE);
  }
  else
  {
    return (init_static_vars_only == 1) ? int_for_C_gk_limber(amin, (void*) ar) :
      like.high_def_integration == 2 ?
      int_gsl_integrate_high_precision(int_for_C_gk_limber, (void*) ar, amin, amax, NULL, 
        GSL_WORKSPACE_SIZE) :
      like.high_def_integration == 1 ?
      int_gsl_integrate_medium_precision(int_for_C_gk_limber, (void*) ar, amin, amax, NULL, 
        GSL_WORKSPACE_SIZE) :
      int_gsl_integrate_low_precision(int_for_C_gk_limber, (void*) ar, amin, amax, NULL, 
        GSL_WORKSPACE_SIZE);
  }
}

double C_gk_tomo_limber(double l, int ni)
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double** table;

  const int NSIZE = tomo.clustering_Nbin;
  const int nell = Ntable.N_ell;
  const double lnlmin = log(fmax(limits.LMIN_tab, 1.0));
  const double lnlmax = log(fmax(limits.LMAX, limits.LMAX_hankel) + 1);
  const double dlnl = (lnlmax - lnlmin)/((double) nell - 1.0);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for (int i=0; i<NSIZE; i++) 
    {
      table[i] = (double*) malloc(sizeof(double)*nell);
    }
  }
  if (recompute_gk(C, G, N))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      const int k = 0;
      const double lnl = lnlmin;
      double init = C_gk_tomo_limber_nointerp(exp(lnl), k, use_linear_ps_limber, 1);
    }
    #pragma GCC diagnostic pop
    #pragma omp parallel for collapse(2)
    for (int k=0; k<NSIZE; k++)
    {
      for (int i=0; i<nell; i++)
      {
        const double lnl = lnlmin + i*dlnl;
        table[k][i]= log(C_gk_tomo_limber_nointerp(exp(lnl), k, use_linear_ps_limber, 0));
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
    update_galpara(&G);
  }
  if (ni < -1 || ni > tomo.clustering_Nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }
  const double lnl = log(l);
  if (lnl < lnlmin)
  {
    log_warn("l = %e < l_min = %e. Extrapolation adopted", l, exp(lnlmin));
  }
  if (lnl > lnlmax)
  {
    log_warn("l = %e > l_max = %e. Extrapolation adopted", l, exp(lnlmax));
  }
  const int q =  ni; 
  if (q < 0 || q > NSIZE - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  } 
  const double f1 = exp(interpol(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 1, 1));
  return isnan(f1) ? 0.0 : f1;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_ks_limber(double a, void* params)
{
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }

  double* ar = (double*) params;
  const int ni = (int) ar[0];
  if (ni < -1 || ni > tomo.shear_Nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }
  const double l = ar[1];
  const int use_linear_ps = (int) ar[2];
  
  const double ell = l + 0.5;  
  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double PK = use_linear_ps == 1 ? p_lin(k,a) : Pdelta(k,a);

  const double WS1 = W_source(a, ni, hoverh0);
  const double WK1 = W_kappa(a, fK, ni);
  const double WK2 = W_k(a, fK);

  const double ell_prefactor1 = l*(l + 1.)/(ell*ell); // prefactor correction (1812.05995 eqs 74-79)
  const double tmp = (l - 1.)*l*(l + 1.)*(l + 2.);    // prefactor correction (1812.05995 eqs 74-79)
  const double ell_prefactor2 = (tmp > 0) ? sqrt(tmp)/(ell*ell) : 0.0; 

  double res = 0;
  switch(like.IA)
  {
    case 0:
    {
      res = WK1*WK2;

      break;
    }
    case 1:
    {
      const double norm = A_IA_Joachimi(a)*cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a;
      res = (-WS1*WK2*norm + WK1*WK2);

      break;
    }
    case 3:
    {
      const double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a;
      res = (-WS1*WK2*norm + WK1*WK2);

      break;
    }
    case 4:
    {
      const double norm = (cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a)*
        nuisance.A_ia*pow(1.0/(a*nuisance.oneplusz0_ia), nuisance.eta_ia);
      res = -(WS1*WK2*norm + WK1*WK2);

      break;
    }
    default:
    {
      log_fatal("like.IA = %d not supported", like.IA);
      exit(1);
    }
  }

  return (res*PK*chidchi.dchida/(fK*fK))*ell_prefactor1*ell_prefactor2;
}

double C_ks_tomo_limber_nointerp(double l, int ni, int use_linear_ps,
const int init_static_vars_only)
{
  if (ni < -1 || ni > tomo.shear_Nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }
  double ar[3] = {(double) ni, l, (double) use_linear_ps};
  const double amin = amin_source(ni);
  const double amax = amax_source(ni);
  if (!(amin>0) || !(amin<1) || !(amax>0) || !(amax<1)) 
  {
    log_fatal("0 < amin/amax < 1 not true");
    exit(1);
  }
  return (init_static_vars_only == 1) ? int_for_C_ks_limber(amin, (void*) ar) :
    like.high_def_integration == 2 ?
    int_gsl_integrate_high_precision(int_for_C_ks_limber, (void*) ar, amin, amax, NULL,
     GSL_WORKSPACE_SIZE) :
    like.high_def_integration == 1 ?
    int_gsl_integrate_medium_precision(int_for_C_ks_limber, (void*) ar, amin, amax, NULL,
     GSL_WORKSPACE_SIZE) :
    int_gsl_integrate_low_precision(int_for_C_ks_limber, (void*) ar, amin, amax, NULL,
     GSL_WORKSPACE_SIZE);
}

double C_ks_tomo_limber(double l, int ni)
{
  static cosmopara C;
  static nuisancepara N;
  static double** table;
  static double* sig;
  static int osc[100];

  const int NSIZE = tomo.shear_Nbin;
  const int nell = Ntable.N_ell;
  const double lnlmin = log(fmax(limits.LMIN_tab, 1.0));
  const double lnlmax = log(fmax(limits.LMAX, limits.LMAX_hankel) + 1);
  const double dlnl = (lnlmax - lnlmin)/((double) nell - 1.0);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for (int i=0; i<NSIZE; i++) 
    {
      table[i] = (double*) malloc(sizeof(double)*nell);
    }
    sig = (double*) malloc(sizeof(double)*NSIZE);
  }
  if (recompute_ks(C, N))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init = C_ks_tomo_limber_nointerp(exp(lnlmin), 0, use_linear_ps_limber, 1);
    }
    #pragma GCC diagnostic pop
    #pragma omp parallel for collapse(2)
    for (int k=0; k<NSIZE; k++)
    {
      for (int i=0; i<nell; i++)
      {
        const double lnl = lnlmin + i*dlnl;
        table[k][i] = C_ks_tomo_limber_nointerp(exp(lnl), k, use_linear_ps_limber, 0);
      }
    }
    #pragma omp parallel for
    for (int k=0; k<NSIZE; k++)
    {
      sig[k] = 1.;
      osc[k] = 0;
      if (C_ks_tomo_limber_nointerp(500., k, use_linear_ps_limber, 0) < 0)
      {
        sig[k] = -1.;
      }
      for (int i=0; i<nell; i++)
      {
        if (table[k][i]*sig[k] < 0.)
        {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0)
      {
        for (int i=0; i<nell; i++)
        {
          table[k][i] = log(sig[k]*table[k][i]);
        }
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  } 
  if (ni < 0 || ni > tomo.shear_Nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d (max %d)", ni, tomo.shear_Nbin);
    exit(1);
  }
  const double lnl = log(l);
  if (lnl < lnlmin)
  {
    log_warn("l = %e < l_min = %e. Extrapolation adopted", l, exp(lnlmin));
  }
  if (lnl > lnlmax)
  {
    log_warn("l = %e > l_max = %e. Extrapolation adopted", l, exp(lnlmax));
  }
  const int q =  ni; 
  if (q < 0 || q > NSIZE - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  } 
  double f1 = 0.0;
  if (osc[ni] == 0)
  {
    f1 = sig[ni]*exp(interpol_fitslope(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 1));
  }
  else if (osc[ni] == 1)
  {
    f1 = interpol_fitslope(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 1);
  }
  else
  {
    log_fatal("internal logic error in selecting osc[ni]");
    exit(1);
  }
  return isnan(f1) ? 0.0 : f1;
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
  const int use_linear_ps = (int) ar[1];
  
  struct chis chidchi = chi_all(a);
  const double ell = l + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double WK = W_k(a, fK);
  const double PK = use_linear_ps == 1 ? p_lin(k,a) : Pdelta(k,a);
  
  const double ell_prefactor = l*(l + 1.0)/(ell*ell);  // prefac correction (1812.05995 eqs 74-79)

  return WK*WK*PK*(chidchi.dchida/(fK*fK))*ell_prefactor*ell_prefactor;
}

double C_kk_limber_nointerp(double l, int use_linear_ps, const int init_static_vars_only)
{
  double ar[2] = {l, (double) use_linear_ps};
  const double amin = limits.a_min*(1. + 1.e-5);
  const double amax = 0.99999;
  return (init_static_vars_only == 1) ? int_for_C_kk_limber(amin, (void*) ar) :
    like.high_def_integration > 0 ?
    int_gsl_integrate_high_precision(int_for_C_kk_limber, (void*) ar, amin, amax, 
      NULL, GSL_WORKSPACE_SIZE) :
    int_gsl_integrate_medium_precision(int_for_C_kk_limber, (void*) ar, amin, amax, 
      NULL, GSL_WORKSPACE_SIZE);
}

double C_kk_limber(double l)
{
  static cosmopara C;
  static double* table;
  
  const int nell = Ntable.N_ell;
  const double lnlmin = log(fmax(limits.LMIN_tab, 1.0));
  const double lnlmax = log(fmax(limits.LMAX, limits.LMAX_hankel) + 1);
  const double dlnl = (lnlmax - lnlmin)/((double) nell - 1.0);

  if (table == 0)
  {
    table = (double*) malloc(sizeof(double)*nell);
  }
  if (recompute_cosmo3D(C))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init = C_kk_limber_nointerp(exp(lnlmin), use_linear_ps_limber, 1);
    }
    #pragma GCC diagnostic pop
    #pragma omp parallel for
    for (int i=0; i<nell; i++)
    {
      const double lnl = lnlmin + i*dlnl;
      table[i] = log(C_kk_limber_nointerp(exp(lnl), use_linear_ps_limber, 0));
    }
    update_cosmopara(&C);
  }
  const double lnl = log(l);
  if (lnl < lnlmin)
  {
    log_warn("l = %e < l_min = %e. Extrapolation adopted", l, exp(lnlmin));
  }
  if (lnl > lnlmax)
  {
    log_warn("l = %e > l_max = %e. Extrapolation adopted", l, exp(lnlmax));
  }
  const double f1 = exp(interpol(table, nell, lnlmin, lnlmax, dlnl, lnl, 1, 1));
  return isnan(f1) ? 0.0 : f1;
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
  if (nl < 0 || nl > tomo.clustering_Nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", nl);
    exit(1);
  }
  const double l = ar[1];
  const int use_linear_ps = (int) ar[2];
  const int use_2h_only = use_linear_ps;

  const double ell = l + 0.5;
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double z = 1./a - 1.;

  const double b1 = gbias.b1_function(z, nl);
  const double bmag = gbias.b_mag[nl];  

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
    {
      res *= WGAL*b1 + WMAG*ell_prefactor*bmag;
    }
    const double PK = p_my(k, a, use_2h_only);
    res *= PK;
  }

  return res*chidchi.dchida/(fK*fK);
}

double C_gy_tomo_limber_nointerp(double l, int ni, int use_linear_ps, 
const int init_static_vars_only)
{
  double ar[3] = {(double)ni, l, (double) use_linear_ps};
  const double amin = amin_lens(ni);
  const double amax = 0.99999;

  if (gbias.b2[ni] || gbias.b2[ni]) 
  {
    log_fatal("b2 not supported in C_gy_nointerp");
    exit(1);
    return 0.0; // avoid gcc warning
  }
  else
  {
    return  (init_static_vars_only == 1) ? int_for_C_gy_tomo_limber(amin, (void*) ar) :
      like.high_def_integration > 0 ?
      int_gsl_integrate_high_precision(int_for_C_gy_tomo_limber, (void*) ar, amin, amax, 
        NULL, GSL_WORKSPACE_SIZE) :
      int_gsl_integrate_medium_precision(int_for_C_gy_tomo_limber, (void*) ar, amin, amax, 
        NULL, GSL_WORKSPACE_SIZE);
  }
}

double C_gy_tomo_limber(double l, int ni)
{
  static cosmopara C;
  static nuisancepara N;
  static ynuisancepara N2;
  static galpara G;
  static double** table;

  const int NSIZE = tomo.clustering_Nbin;
  const int nell = Ntable.N_ell;
  const double lnlmin = log(fmax(limits.LMIN_tab, 1.0));
  const double lnlmax = log(fmax(limits.LMAX, limits.LMAX_hankel) + 1);
  const double dlnl = (lnlmax - lnlmin)/((double) nell - 1.0);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for (int i=0; i<NSIZE; i++) 
    {
      table[i] = (double*) malloc(sizeof(double)*nell);
    }
  }
  if (recompute_gy(C, G, N, N2))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      const int k = 0;
      const double lnl = lnlmin;
      double init = C_gy_tomo_limber_nointerp(exp(lnl), k, use_linear_ps_limber, 1);
    }
    #pragma GCC diagnostic pop
    #pragma omp parallel for collapse(2)
    for (int k=0; k<NSIZE; k++)
    {
      for (int i=0; i<nell; i++)
      {
        const double lnl = lnlmin + i*dlnl;
        table[k][i]= log(C_gy_tomo_limber_nointerp(exp(lnl), k, use_linear_ps_limber, 0));
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
    update_galpara(&G);
    update_ynuisance(&N2);
  }
  
  const int q =  ni; 
  if (q < 0 || q > NSIZE - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  } 
  
  const double lnl = log(l);
  if (lnl < lnlmin)
  {
    log_warn("l = %e < l_min = %e. Extrapolation adopted", l, exp(lnlmin));
  }
  if (lnl > lnlmax)
  {
    log_warn("l = %e > l_max = %e. Extrapolation adopted", l, exp(lnlmax));
  }
  
  const double f1 = exp(interpol(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 1, 1));
  return isnan(f1) ? 0.0 : f1;
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
  if (ni < 0 || ni > tomo.shear_Nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }
  const double l = ar[1];
  const int use_linear_ps = (int) ar[2];
  const int use_2h_only = use_linear_ps;
  
  const double ell = l + 0.5;
  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double fK = f_K(chidchi.chi);
  const double k  = ell/fK;
  
  const double PK = p_my(k, a, use_2h_only);

  const double WS = W_source(a, ni, hoverh0);
  const double WK = W_kappa(a, fK, ni);
  const double WY = W_y(a);

  const double tmp = (l - 1.0)*l*(l + 1.0)*(l + 2.0); // prefactor correction (1812.05995 eqs 74-79)
  const double ell_prefactor2 = (tmp > 0) ? sqrt(tmp)/(ell*ell) : 0.0;

  double res = 0;
  switch(like.IA)
  {
    case 0:
    {
      res = WS*WY;

      break;
    }
    case 1:
    {
      const double norm = A_IA_Joachimi(a)*cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a;
      
      res = WK*WY - WS*WY*norm;
      
      break;
    }
    case 4:
    {
      const double norm = (cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac_a)*
        nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);
      
      res = WK*WY - WS*WY*norm;

      break;
    }
    default:
    {
      log_fatal("like.IA = %d not supported", like.IA);
      exit(1);
    }
  }

  return res*PK*(chidchi.dchida/(fK*fK))*ell_prefactor2;
}

double C_ys_tomo_limber_nointerp(double l, int ni, int use_linear_ps, 
const int init_static_vars_only)
{
  double ar[3] = {(double) ni, l, (double) use_linear_ps};
  const double amin = amin_source(ni);
  const double amax = 0.99999;
   
  return  (init_static_vars_only == 1) ? int_for_C_ys_tomo_limber(amin, (void*) ar) :
    like.high_def_integration > 0 ?
    int_gsl_integrate_high_precision(int_for_C_ys_tomo_limber, (void*) ar, amin, amax, 
      NULL, GSL_WORKSPACE_SIZE) :
    int_gsl_integrate_medium_precision(int_for_C_ys_tomo_limber, (void*) ar, amin, amax, 
      NULL, GSL_WORKSPACE_SIZE);
}

double C_ys_tomo_limber(double l, int ni)
{
  static cosmopara C;
  static nuisancepara N;
  static ynuisancepara N2;

  static double** table;
  static double* sig;
  static int osc[100];

  const int NSIZE = tomo.shear_Nbin;
  const int nell = Ntable.N_ell;
  const double lnlmin = log(fmax(limits.LMIN_tab, 1.0));
  const double lnlmax = log(fmax(limits.LMAX, limits.LMAX_hankel) + 1);
  const double dlnl = (lnlmax - lnlmin)/((double) nell - 1.0);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for (int i=0; i<NSIZE; i++) 
    {
      table[i] = (double*) malloc(sizeof(double)*nell);
    }
    sig = (double*) malloc(sizeof(double)*NSIZE);
  }
  if (recompute_ys(C, N, N2))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init = C_ys_tomo_limber_nointerp(exp(lnlmin), 0, use_linear_ps_limber, 1);
    }
    #pragma GCC diagnostic pop
    #pragma omp parallel for collapse(2)
    for (int k=0; k<NSIZE; k++)
    {
      for (int i=0; i<nell; i++)
      {
        table[k][i] = C_ys_tomo_limber_nointerp(exp(lnlmin + i*dlnl), k, use_linear_ps_limber, 0);
      }
    }
    #pragma omp parallel for
    for (int k=0; k<NSIZE; k++)
    {
      sig[k] = 1.;
      osc[k] = 0;
      if (C_ys_tomo_limber_nointerp(500., k, use_linear_ps_limber, 0) < 0)
      {
        sig[k] = -1.;
      }
      for (int i=0; i<nell; i++)
      {
        if (table[k][i]*sig[k] < 0.)
        {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0)
      {
        for (int i=0; i<nell; i++)
        {
          table[k][i] = log(sig[k]*table[k][i]);
        }
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
    update_ynuisance(&N2);
  } 
  if (ni < 0 || ni > tomo.shear_Nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d (max %d)", ni, tomo.shear_Nbin);
    exit(1);
  }

  const int q =  ni; 
  if (q < 0 || q > NSIZE - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  } 
  
  const double lnl = log(l);
  if (lnl < lnlmin)
  {
    log_warn("l = %e < l_min = %e. Extrapolation adopted", l, exp(lnlmin));
  }
  if (lnl > lnlmax)
  {
    log_warn("l = %e > l_max = %e. Extrapolation adopted", l, exp(lnlmax));
  }

  double f1 = 0.0;
  if (osc[ni] == 0)
  {
    f1 = sig[ni]*exp(interpol_fitslope(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 1));
  }
  else if (osc[ni] == 1)
  {
    f1 = interpol_fitslope(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 1);
  }
  else
  {
    log_fatal("internal logic error in selecting osc[ni]");
    exit(1);
  }
  return isnan(f1) ? 0.0 : f1;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_ky_limber(double a, void*params)
{
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }

  double *ar = (double*) params;
  const double l = ar[0];
  const int use_linear_ps = (int) ar[1];
  const int use_2h_only = use_linear_ps;

  const double ell = l + 0.5;
  struct chis chidchi = chi_all(a);
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  
  const double PK = p_my(k, a, use_2h_only);
  
  const double WK = W_k(a, fK);
  const double WY = W_y(a);

  const double ell_prefactor = l*(l + 1.0)/(ell*ell); // prefactor correction (1812.05995 eqs 74-79)

  return (WK*WY*PK*chidchi.dchida/(fK*fK))*ell_prefactor;
}

double C_ky_limber_nointerp(double l, int use_linear_ps, const int init_static_vars_only)
{
  double ar[2] = {l, (double) use_linear_ps};
  const double amin = limits.a_min_hm;
  const double amax = 1.0 - 1.e-5;

  return (init_static_vars_only == 1) ? int_for_C_ky_limber(amin, (void*) ar) :
    like.high_def_integration > 0 ?
    int_gsl_integrate_high_precision(int_for_C_ky_limber, (void*) ar, amin, amax, 
      NULL, GSL_WORKSPACE_SIZE) :
    int_gsl_integrate_medium_precision(int_for_C_ky_limber, (void*) ar, amin, amax, 
      NULL, GSL_WORKSPACE_SIZE);
}

double C_ky_limber(double l)
{
  static cosmopara C;
  static ynuisancepara N;
  static double* table;
  
  const int nell = Ntable.N_ell;
  const double lnlmin = log(fmax(limits.LMIN_tab, 1.0));
  const double lnlmax = log(fmax(limits.LMAX, limits.LMAX_hankel) + 1);
  const double dlnl = (lnlmax - lnlmin)/((double) nell - 1.0);

  if (table == 0)
  {
    table = (double*) malloc(sizeof(double)*nell);
  }
  if (recompute_ky(C, N))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init = C_ky_limber_nointerp(exp(lnlmin), use_linear_ps_limber, 1);
    }
    #pragma GCC diagnostic pop
    #pragma omp parallel for
    for (int i=0; i<nell; i++)
    {
      table[i] = log(C_ky_limber_nointerp(exp(lnlmin + i*dlnl), use_linear_ps_limber, 0));
    }
    update_ynuisance(&N);
    update_cosmopara(&C);
  }
  const double lnl = log(l);
  if (lnl < lnlmin)
  {
    log_warn("l = %e < l_min = %e. Extrapolation adopted", l, exp(lnlmin));
  }
  if (lnl > lnlmax)
  {
    log_warn("l = %e > l_max = %e. Extrapolation adopted", l, exp(lnlmax));
  }
  const double f1 = exp(interpol(table, nell, lnlmin, lnlmax, dlnl, lnl, 1, 1));
  return isnan(f1) ? 0.0 : f1;
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
  const int use_linear_ps = (int) ar[1];
  const int use_2h_only = use_linear_ps;

  const double ell = l + 0.5;
  struct chis chidchi = chi_all(a);
  const double fK = f_K(chidchi.chi);
  const double k  = ell/fK;

  const double PK = p_yy(k, a, use_2h_only);
  
  const double WY = W_y(a);

  return WY*WY*PK*chidchi.dchida/(fK*fK);
}

double C_yy_limber_nointerp(double l, int use_linear_ps, const int init_static_vars_only)
{
  double ar[2] = {l, (double) use_linear_ps};
  const double amin = limits.a_min;
  const double amax = 1.0 - 1.e-5;

  return (init_static_vars_only == 1) ? int_for_C_yy_limber(amin, (void*) ar) :
    like.high_def_integration > 0 ?
    int_gsl_integrate_high_precision(int_for_C_yy_limber, (void*) ar, amin, amax, 
      NULL, GSL_WORKSPACE_SIZE) :
    int_gsl_integrate_medium_precision(int_for_C_yy_limber, (void*) ar, amin, amax, 
      NULL, GSL_WORKSPACE_SIZE);
}

double C_yy_limber(double l)
{
  static cosmopara C;
  static ynuisancepara N;
  static double* table;
  
  const int nell = Ntable.N_ell;
  const double lnlmin = log(fmax(limits.LMIN_tab, 1.0));
  const double lnlmax = log(fmax(limits.LMAX, limits.LMAX_hankel) + 1.0);
  const double dlnl = (lnlmax - lnlmin)/((double) nell - 1.0);

  if (table == 0)
  {
    table = (double*) malloc(sizeof(double)*nell);
  }
  if (recompute_yy(C, N))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init = C_yy_limber_nointerp(exp(lnlmin), use_linear_ps_limber, 1);
    }
    #pragma GCC diagnostic pop
    #pragma omp parallel for
    for (int i=0; i<nell; i++)
    {
      table[i] = log(C_yy_limber_nointerp(exp(lnlmin + i*dlnl), use_linear_ps_limber, 0));
    }
    update_ynuisance(&N);
    update_cosmopara(&C);
  }
  const double lnl = log(l);
  if (lnl < lnlmin)
  {
    log_warn("l = %e < l_min = %e. Extrapolation adopted", l, exp(lnlmin));
  }
  if (lnl > lnlmax)
  {
    log_warn("l = %e > l_max = %e. Extrapolation adopted", l, exp(lnlmax));
  }
  const double f1 = exp(interpol(table, nell, lnlmin, lnlmax, dlnl, lnl, 1, 1));
  return isnan(f1) ? 0.0 : f1;
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Non Limber (Angular Power Spectrum)
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// --------------------------------------------------------------------------------------
// Galaxy Clustering
// --------------------------------------------------------------------------------------

void f_chi_for_Psi_cl(double *const chi, int Nchi, double *const f_chi, const int ni,
const double zmin, const double zmax)
{ // Integrand for galaxy density
  if (ni < -1 || ni > tomo.clustering_Nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }
  const double real_coverH0 = cosmology.coverH0/cosmology.h0; // unit Mpc
  for (int i=0; i<Nchi; i++)
  {
    const double a = a_chi(chi[i]/real_coverH0 /* convert unit to c/H0 */);
    const double z = 1.0/a - 1.;
    const double tmp1 = pf_photoz(z, ni);
    const double pf = (tmp1 < 0.) ? 0 : tmp1; // get rid of unphysical negatives
    f_chi[i] = gbias.b1_function(z, ni)*chi[i]*pf*growfac(a)*hoverh0(a)/real_coverH0;
    if ((z < zmin) || (z > zmax))
    {
      f_chi[i] = 0.;
    }
  }
}

void f_chi_for_Psi_cl_RSD(double *const chi, int Nchi, double *const f_chi, const int ni,
const double zmin, const double zmax)
{ // Integrand for galaxy density RSD
  if (ni < -1 || ni > tomo.clustering_Nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }
  const double real_coverH0 = cosmology.coverH0 / cosmology.h0;
  for (int i=0; i<Nchi; i++) 
  {
    const double a = a_chi(chi[i]/real_coverH0 /* convert unit to c/H0 */);
    const double z = 1.0 / a - 1.0;
    const double tmp1 = pf_photoz(z, ni);
    const double pf = (tmp1 < 0.) ? 0 : tmp1; // get rid of unphysical negatives
    struct growths tmp2 = growfac_all(a);
    f_chi[i] = -chi[i]*pf*tmp2.D*tmp2.f*hoverh0(a)/real_coverH0;
    if ((z < zmin) || (z > zmax))
    {
      f_chi[i] = 0.;
    } 
  }
}

void f_chi_for_Psi_cl_Mag(double *const chi, int Nchi, double *const f_chi, const int ni,
const double zmax)
{ // Integrand for lensing magnification of galaxy density
  if (ni < -1 || ni > tomo.clustering_Nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }
  const double real_coverH0 = cosmology.coverH0 / cosmology.h0;
  for (int i=0; i<Nchi; i++)
  {
    const double a = a_chi(chi[i]/real_coverH0 /* convert unit to c/H0 */);
    const double z = 1. / a - 1.;
    const double fK = f_K(chi[i]/real_coverH0 /* convert unit to c/H0 */);
    const double wmag = W_mag(a, fK, ni);
    const double window_M = wmag/fK/(real_coverH0*real_coverH0);
    f_chi[i] = window_M * growfac(a); // unit [Mpc^-2]
    if (z > zmax)
    {
      f_chi[i] = 0.;
    } 
  }
}

void C_cl_tomo(int L, const int ni, const int nj, double *const Cl, double dev, double tol)
{
  if (ni < -1 || ni > tomo.clustering_Nbin - 1 || nj < -1 || nj > tomo.clustering_Nbin - 1)
  {
    log_fatal("error in selecting bin number (ni, nj) = [%d,%d]", ni, nj);
    exit(1);
  }
  if (ni != nj)
  {
    log_fatal("Cocoa disabled cross-spectrum w_gg");
    exit(1);
  }

  static double** k1;
  static double** k2;
  static double** Fk1;
  static double** Fk2;
  static double** Fk1_Mag;
  static double** Fk2_Mag;
  static double* chi_ar;

  const int Nell_block = Ntable.NL_Nell_block;
  const int Nchi = Ntable.NL_Nchi;  
  int ell_ar[Nell_block];
  double f1_chi[Nchi];
  double f2_chi[Nchi];
  double f1_chi_RSD[Nchi];
  double f2_chi_RSD[Nchi];
  double f1_chi_Mag[Nchi];
  double f2_chi_Mag[Nchi];
  const double real_coverH0 = cosmology.coverH0/cosmology.h0;
  const double chi_min = chi(1./(1.0 + 0.002))*real_coverH0; // DIMENSIONELESS
  const double chi_max = chi(1./(1.0 + 4.0))*real_coverH0;   // DIMENSIONELESS
  const double dlnchi = log(chi_max/chi_min) / ((double) Nchi - 1.0);
  const double dlnk = dlnchi;

  if (k1 == 0) 
  {
    k1 = (double**) malloc(Nell_block * sizeof(double*));
    k2 = (double**) malloc(Nell_block * sizeof(double*));
    Fk1 = (double**) malloc(Nell_block * sizeof(double*));
    Fk2 = (double**) malloc(Nell_block * sizeof(double*));
    Fk1_Mag = (double**) malloc(Nell_block * sizeof(double*));
    Fk2_Mag = (double**) malloc(Nell_block * sizeof(double*));
    for (int i = 0; i < Nell_block; i++) 
    {
      k1[i] = (double*) malloc(Nchi * sizeof(double));
      k2[i] = (double*) malloc(Nchi * sizeof(double));
      Fk1[i] = (double*) malloc(Nchi * sizeof(double));
      Fk2[i] = (double*) malloc(Nchi * sizeof(double));
      Fk1_Mag[i] = (double*) malloc(Nchi * sizeof(double));
      Fk2_Mag[i] = (double*) malloc(Nchi * sizeof(double));
    }
    chi_ar = (double*) malloc(Nchi * sizeof(double));
  }
  for (int i=0; i < Nchi; i++)
  { // chi_min and chi_max are cosmology dependent
    chi_ar[i] = chi_min * exp(dlnchi * i); 
  }
  #pragma omp parallel for collapse(2)
  for (int i=0; i<Nell_block; i++)
  {
    for (int j=0; j<Nchi; j++)
    {
      k1[i][j] = 0.0;
      k2[i][j] = 0.0;
      Fk1[i][j] = 0.0;
      Fk2[i][j] = 0.0;
      Fk1_Mag[i][j] = 0.0;
      Fk2_Mag[i][j] = 0.0;
    }
  }
  
  {
    const double zmin = tomo.clustering_zmin[ni];
    const double zmax = tomo.clustering_zmax[ni];
    f_chi_for_Psi_cl(chi_ar, Nchi, f1_chi, ni, zmin, zmax);
    f_chi_for_Psi_cl_RSD(chi_ar, Nchi, f1_chi_RSD, ni, zmin, zmax);
    f_chi_for_Psi_cl_Mag(chi_ar, Nchi, f1_chi_Mag, ni, zmax);
  }
  if (ni != nj)
  { 
    const double zmin = tomo.clustering_zmin[nj];
    const double zmax = tomo.clustering_zmax[nj];
    f_chi_for_Psi_cl(chi_ar, Nchi, f2_chi, nj, zmin, zmax);
    f_chi_for_Psi_cl_RSD(chi_ar, Nchi, f2_chi_RSD, nj, zmin, zmax);
    f_chi_for_Psi_cl_Mag(chi_ar, Nchi, f2_chi_Mag, nj, zmax);
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
    for (int i=0; i<Nell_block; i++)
    {
      ell_ar[i] = i + i_block * Nell_block; 
    }
    i_block++;
    if (L >= limits.LMAX_NOLIMBER - Nell_block)
    { //Xiao: break before memory leak in next iteration
      break;
    }
    L = i_block * Nell_block - 1;

    cfftlog_ells(chi_ar, f1_chi, Nchi, &cfg, ell_ar, Nell_block, k1, Fk1);
    cfftlog_ells_increment(chi_ar, f1_chi_RSD, Nchi, &cfg_RSD, ell_ar, Nell_block, k1, Fk1);   
    cfftlog_ells(chi_ar, f1_chi_Mag, Nchi, &cfg_Mag, ell_ar, Nell_block, k1, Fk1_Mag);    
    if (ni != nj)
    {
      cfftlog_ells(chi_ar, f2_chi, Nchi, &cfg, ell_ar, Nell_block, k2, Fk2);
      cfftlog_ells_increment(chi_ar, f2_chi_RSD, Nchi, &cfg_RSD, ell_ar, Nell_block, k2, Fk2);
      cfftlog_ells(chi_ar, f2_chi_Mag, Nchi, &cfg_Mag, ell_ar, Nell_block, k2, Fk2_Mag);
    }

    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
    {
      double init = p_lin(k1[0][0]*real_coverH0, 1.0);
      init = C_gg_tomo_limber_nointerp((double) ell_ar[0], ni, nj, 0, 1);
    }
    #pragma GCC diagnostic pop
    #pragma GCC diagnostic pop
    #pragma omp parallel for
    for (int i=0; i<Nell_block; i++)
    {
      double cl_temp = 0.;
      for (int j=0; j<Nchi; j++)
      {
        const double ell_prefactor = ell_ar[i] * (ell_ar[i] + 1.);
        Fk1[i][j] += gbias.b_mag[ni]*ell_prefactor*Fk1_Mag[i][j]/(k1[i][j]*k1[i][j]);
        if (ni != nj)
        {
          Fk2[i][j] += gbias.b_mag[nj]*ell_prefactor*Fk2_Mag[i][j]/(k2[i][j]*k2[i][j]);
        }
        // ------------------------------------------------------------------------------------
        const double k1cH0 = k1[i][j] * real_coverH0;
        const double PK = p_lin(k1cH0, 1.0);
        const double k1cH03 = k1cH0*k1cH0*k1cH0;
        cl_temp += (ni == nj) ? Fk1[i][j]*Fk1[i][j]*k1cH03*PK : Fk1[i][j]*Fk2[i][j]*k1cH03*PK;
      }
      Cl[ell_ar[i]] = cl_temp * dlnk * 2. / M_PI + 
         C_gg_tomo_limber_nointerp((double) ell_ar[i], ni, nj, 0, 0)
        -C_gg_tomo_limber_nointerp((double) ell_ar[i], ni, nj, 1, 0);
    }
    dev = Cl[L]/C_gg_tomo_limber_nointerp((double) L, ni, nj, 0, 0) - 1;
  }
  L++;

  Cl[limits.LMAX_NOLIMBER] = C_gg_tomo_limber((double) limits.LMAX_NOLIMBER, ni, nj);
  #pragma omp parallel for
  for (int l=L; l<limits.LMAX_NOLIMBER; l++)
  {
    Cl[l] = (l > limits.LMIN_tab) ? C_gg_tomo_limber((double) l, ni, nj) :
      C_gg_tomo_limber_nointerp((double) l, ni, nj, use_linear_ps_limber, 0);
  }
}

// --------------------------------------------------------------------------------------
// Galaxy-Galaxylensing
// --------------------------------------------------------------------------------------


void f_chi_for_Psi_sh(double *const chi, int Nchi, double *const fchi, const int nj,
const double zmax)
{
  if (nj < -1 || nj > tomo.shear_Nbin - 1)
  {
    log_fatal("error in selecting bin number nj = %d", nj);
    exit(1);
  }
  const double real_coverH0 = cosmology.coverH0/cosmology.h0;
  for (int i=0; i<Nchi; i++)
  {
    const double a = a_chi(chi[i] / real_coverH0) ;
    const double z = 1.0/a - 1.0;
    const double fK = f_K(chi[i]/real_coverH0);
    const double wkappa = W_kappa(a, fK, nj);
    fchi[i] = (wkappa/fK/(real_coverH0*real_coverH0))*growfac(a); // unit [Mpc^-2]
    if (z > zmax)
    {
      fchi[i] = 0.;
    }
  }
}

void f_chi_for_Psi_sh_IA(double *const chi, int Nchi, double *const fchi, const int nj,
const double zmin, const double zmax)
{ // TODO: ADD ALL IA POSSIBILITIES
  if (nj < -1 || nj > tomo.shear_Nbin - 1)
  {
    log_fatal("error in selecting bin number nj = %d", nj);
    exit(1);
  }
  const double real_coverH0 = cosmology.coverH0 / cosmology.h0;
  for (int i=0; i<Nchi; i++)
  {
    // first convert unit of chi from Mpc to c/H0
    const double a = a_chi(chi[i] / real_coverH0) ;
    const double z = 1./a - 1.;
    struct chis chidchi = chi_all(a);
    const double hoverh0 = hoverh0v2(a, chidchi.dchida);
    const double fK = f_K(chi[i]/real_coverH0);
    const double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac(a)*
      nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);
    const double tmp1 = W_source(a, (double) nj, hoverh0);
    const double wsource = (tmp1 > 0.) ? tmp1 : 0.;
    fchi[i] = -wsource*norm/fK/(real_coverH0*real_coverH0)*growfac(a); // unit [Mpc^-2]
    if ((z<zmin) || (z>zmax))
    {
      fchi[i] = 0.;
    }
  }
}

void C_gl_tomo(int L, int nl, int ns, double *const Cl, double dev, double tolerance)
{
  if (nl < -1 || nl > tomo.clustering_Nbin - 1 || ns < -1 || ns > tomo.clustering_Nbin - 1)
  {
    log_fatal("error in selecting bin number (ni, nj) = [%d,%d]", nl, ns);
    exit(1);
  }

  static double** k1;
  static double** k2;
  static double** Fk1;
  static double** Fk2;
  static double** Fk1_Mag;
  static double** Fk2_Mag;
  static double* chi_ar;

  const int Nell_block = Ntable.NL_Nell_block;
  const int Nchi = Ntable.NL_Nchi;  
  int ell_ar[Nell_block];
  double f1_chi[Nchi];
  double f1_chi_RSD[Nchi];
  double f1_chi_Mag[Nchi];
  double f2_chi[Nchi];
  double f2_chi_IA_ar[Nchi];
  const double real_coverH0 = cosmology.coverH0 / cosmology.h0;
  const double chi_min = chi(1./(1.0 + 0.002))*real_coverH0; // DIMENSIONELESS
  const double chi_max = chi(1./(1.0 + 4.0))*real_coverH0; // DIMENSIONELESS
  const double dlnchi = log(chi_max/chi_min) / ((double) Nchi - 1.0);
  const double dlnk = dlnchi;

  if (k1 == 0)
  {
    k1 = (double**) malloc(Nell_block * sizeof(double*));
    k2 = (double**) malloc(Nell_block * sizeof(double*));
    Fk1 = (double**) malloc(Nell_block * sizeof(double*));
    Fk2 = (double**) malloc(Nell_block * sizeof(double*));
    Fk1_Mag = (double**) malloc(Nell_block * sizeof(double*));
    Fk2_Mag = (double**) malloc(Nell_block * sizeof(double*));
    for (int i=0; i<Nell_block; i++)
    {
      k1[i] = (double*) malloc(Nchi * sizeof(double));
      k2[i] = (double*) malloc(Nchi * sizeof(double));
      Fk1[i] = (double*) malloc(Nchi * sizeof(double));
      Fk2[i] = (double*) malloc(Nchi * sizeof(double));
      Fk1_Mag[i] = (double*) malloc(Nchi * sizeof(double));
      Fk2_Mag[i] = (double*) malloc(Nchi * sizeof(double));
    }
    chi_ar = (double*) malloc(Nchi * sizeof(double));
  }
  for (int i=0; i<Nchi; i++)
  {
    // chi_min and chi_max may be cosmology dependent
    chi_ar[i] = chi_min * exp(dlnchi * i);
  }
  #pragma omp parallel for
  for (int i=0; i<Nell_block; i++)
  {
    for (int j=0; j<Nchi; j++)
    {
      k1[i][j] = 0.0;
      k2[i][j] = 0.0;
      Fk1[i][j] = 0.0;
      Fk2[i][j] = 0.0;
      Fk1_Mag[i][j] = 0.0;
      Fk2_Mag[i][j] = 0.0;
    }
  }

  {
    const double zmin = tomo.clustering_zmin[nl];
    const double zmax = tomo.clustering_zmax[nl];
    f_chi_for_Psi_cl(chi_ar, Nchi, f1_chi, nl, zmin, zmax);
    f_chi_for_Psi_cl_RSD(chi_ar, Nchi, f1_chi_RSD, nl, zmin, zmax);
    f_chi_for_Psi_cl_Mag(chi_ar, Nchi, f1_chi_Mag, nl, zmax);
  }
  {
    const double zmin = tomo.shear_zmin[ns];
    const double zmax = tomo.shear_zmax[ns];
    f_chi_for_Psi_sh(chi_ar, Nchi, f2_chi, ns, zmax);
    f_chi_for_Psi_sh_IA(chi_ar, Nchi, f2_chi_IA_ar, ns, zmin, zmax);
  }

  for (int j=0; j<Nchi; j++)
  {
    f2_chi[j] += f2_chi_IA_ar[j];
  }

  int i_block = 0;

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
  cfg_RSD.N_pad = 200;
  cfg_RSD.N_extrap_low = 0;
  cfg_RSD.N_extrap_high = 0;

  config cfg_Mag;
  cfg_Mag.nu = 1.;
  cfg_Mag.c_window_width = 0.25;
  cfg_Mag.derivative = 0;
  cfg_Mag.N_pad = 1000;
  cfg_Mag.N_extrap_low = 0;
  cfg_Mag.N_extrap_high = 0;

  config cfg_shear;
  cfg_shear.nu = 1.;
  cfg_shear.c_window_width = 0.25;
  cfg_shear.derivative = 0;
  cfg_shear.N_pad = 1000.;
  cfg_shear.N_extrap_low = 0;
  cfg_shear.N_extrap_high = 0;

  while ((fabs(dev) > tolerance) && (L < limits.LMAX_NOLIMBER))
  {
    for (int i=0; i<Nell_block; i++)
    {
      ell_ar[i] = i + i_block*Nell_block;
    }
    i_block++;
    if (L >= limits.LMAX_NOLIMBER - Nell_block)
    { //Xiao: break before memory leak in next iteration
      break;
    }
    L = i_block*Nell_block -1 ;

    cfftlog_ells(chi_ar, f1_chi, Nchi, &cfg, ell_ar, Nell_block, k1, Fk1);
    cfftlog_ells_increment(chi_ar, f1_chi_RSD, Nchi, &cfg_RSD, ell_ar, Nell_block, k1, Fk1);
    cfftlog_ells(chi_ar, f1_chi_Mag, Nchi, &cfg_Mag, ell_ar, Nell_block, k1, Fk1_Mag);
    cfftlog_ells(chi_ar, f2_chi, Nchi, &cfg_shear, ell_ar, Nell_block, k2, Fk2);

    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
    {
      double init = p_lin(k1[0][0]*real_coverH0, 1.0);
      init = C_gs_tomo_limber_nointerp((double) ell_ar[0], nl, ns, 1, 1);
    }
    #pragma GCC diagnostic pop
    #pragma GCC diagnostic pop
    #pragma omp parallel for
    for (int i=0; i<Nell_block; i++)
    {
      double cl_temp = 0.;
      for (int j=0; j<Nchi; j++)
      {
        const double ell_prefactor = ell_ar[i]*(ell_ar[i]+1);
        Fk1[i][j] += gbias.b_mag[nl]*(ell_prefactor/(k1[i][j]*k1[i][j])*Fk1_Mag[i][j]);
        double ell_prefactor2 = (ell_ar[i]-1.)*ell_ar[i]*(ell_ar[i]+1.)*(ell_ar[i]+2.);
        if (ell_prefactor2 <= 0.)
        {
          ell_prefactor2 = 0.;
        }
        else
        {
          ell_prefactor2 = sqrt(ell_prefactor2);
        }
        Fk2[i][j] *= (ell_prefactor2/(k1[i][j]*k1[i][j]));
        // ------------------------------------------------------------------------------------
        const double k1_cH0 = k1[i][j] * real_coverH0;
        cl_temp += Fk1[i][j]*Fk2[i][j]*k1_cH0*k1_cH0*k1_cH0*p_lin(k1_cH0, 1.0);
      }
      Cl[ell_ar[i]] = cl_temp * dlnk * 2./M_PI +
        C_gs_tomo_limber_nointerp((double) ell_ar[i], nl, ns, 0, 0)
       -C_gs_tomo_limber_nointerp((double) ell_ar[i], nl, ns, 1, 0);
    }
    dev = Cl[L]/C_gs_tomo_limber_nointerp((double) L, nl, ns, 0, 0) - 1;
  }
  L++;

  Cl[limits.LMAX_NOLIMBER] = C_gs_tomo_limber((double) limits.LMAX_NOLIMBER, nl, ns);
  #pragma omp parallel for
  for (int l=L; l<limits.LMAX_NOLIMBER; l++)
  {
    Cl[l] = (l > limits.LMIN_tab) ? C_gs_tomo_limber((double) l, nl, ns) :
      C_gs_tomo_limber_nointerp((double) l, nl, ns, use_linear_ps_limber, 0);
  }      
}
