#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../cfftlog/cfftlog.h"

#include "bias.h"
#include "basics.h"
#include "cosmo3D.h"
#include "cluster_util.h"
#include "cosmo2D_cluster.h"
#include "radial_weights.h"
#include "recompute.h"
#include "redshift_spline.h"
#include "structs.h"

#include "log.c/src/log.h"

static int GSL_WORKSPACE_SIZE = 250;
static int use_linear_ps_limber = 0; /* 0 or 1 */
static double w_l_min = 0.0001;
static double w_l_max = 5.0e6; 

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Correlation Functions (real space) - flat sky
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double w_gammat_cluster_tomo_flatsky(double theta, int nl, int ni, int nj, int limber) 
{ // nl = lambda_obs bin, ni = cluster redshift bin, nj = source redshift bin
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double** table;
 
  const int ntheta = Ntable.N_thetaH;
  const int nlsize = Cluster.N200_Nbin;
  const int ngammat_size = tomo.cgl_Npowerspectra;
  const int NSIZE = nlsize*ngammat_size;

  const double l_min = w_l_min;
  const double l_max = w_l_max;
  const double loglmax = log(l_max);
  const double loglmin = log(l_min);
  const double dlnl = (loglmax-loglmin)/(1.0*ntheta - 1.);
  const double lnrc = 0.5*(loglmax + loglmin);
  const double nc = ntheta/2.0 + 1;

  const double logthetamin = (nc-ntheta+1)*dlnl-lnrc;
  const double logthetamax = nc*dlnl-lnrc;
  const double dlntheta = (lnthetamax - lnthetamin)/((double) ntheta);
  const double lntheta = log(theta);

  if (table == 0)
  {
    table  = create_double_matrix(0, NSIZE, 0, ntheta - 1);
  }

  if (recompute_cs(C, G, N))
  {
    if(limber == 1)
    {
      typedef fftw_complex fftwZ;

      // go to log-Fourier-space
      fftwZ** flP = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
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

        // Power spectrum on logarithmic bins
        for(int i=0; i<nlsize; i++)
        { // loop in lambda_obs - (cannot thread loops in lambda_obs  - see cluster_utils)
          { 
            const int j = 0;
            const int ZC = ZC(j);
            const int ZS = ZSC(j);
            const int q = nlsize*i + j;
            {
              const double p = 0;
              const double l = exp(lnrc + (p - nc)*dlnl);
              lP[q][p] = l*C_gs_tomo_limber(l, i, ZC, ZS);
            }
            #pragma omp parallel for
            for(int p=1; p<ntheta; p++)
            {
              const double l = exp(lnrc + (p - nc)*dlnl);
              lP[q][p] = l*C_gs_tomo_limber(l, i, ZC, ZS);
            }
          }
          #pragma omp parallel for
          for(int j=1; j<ngammat_size; j++)
          { 
            const int ZC = ZC(j);
            const int ZS = ZSC(j);
            const int q = nlsize*i + j;
            for(int p=0; p<ntheta; p++)
            {
              const double l = exp(lnrc + (p - nc)*dlnl);
              lP[j][i] = l*C_gs_tomo_limber(l, i, ZC, ZS);
            }
          }    
        }

        #pragma omp parallel for
        for (int j=0; j<NSIZE; j++)
        { // Execute FFTW in parallel (thread-safe)
          fftw_execute(plan[j]); 
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
        lP[j] = (double*) malloc(ntheta*sizeof(double));
        kernel[j] = (fftwZ*) fftw_malloc(sizeof(fftwZ));
        conv[j] = (fftwZ*) fftw_malloc((ntheta/2+1)*sizeof(fftwZ));
        plan[j] = fftw_plan_dft_c2r_1d(ntheta, conv[j], lP[j], FFTW_ESTIMATE);
        tab[j] = create_double_matrix(0, 1, 0, ntheta - 1);
      }

      #pragma omp parallel for
      for (int j=0; j<NSIZE; j++)
      { 
        double arg[2];
        arg[0] = 0; // bias
        arg[1] = 2; // order of Bessel function

        for(int i=0; i<(ntheta/2+1); i++)
        {
          const double kk = 2*constants.pi*i/(dlnl*ntheta);
          hankel_kernel_FT(kk, kernel[j], arg, 2);
          conv[j][i][0] = flP[j][i][0]*kernel[j][0][0] - flP[j][i][1]*kernel[j][0][1];
          conv[j][i][1] = flP[j][i][1]*kernel[j][0][0] + flP[j][i][0]*kernel[j][0][1];
        }

        // force Nyquist- and 0-frequency-components to be double
        conv[j][0][1] = 0;
        conv[j][ntheta/2][1] = 0;

        fftw_execute(plan[j]); // Execute FFTW in parallel (thread-safe)

        for(int k=0; k<ntheta; k++)
        {
          const double t = exp((nc-k)*dlnl-lnrc); 
          tab[j][0][ntheta-k-1] = lP[j][k]/(t*2*constants.pi*ntheta);
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
        free_double_matrix(tab[j], 0, 1, 0, ntheta);
      }
      free(flP);
      free(lP);
      free(conv);
      free(kernel);
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
  
  if (!test_zoverlap_c(ni, nj)) 
  {
    return 0.0;
  }

  if (lntheta < lnthetamin || lntheta > lnthetamax)
  {
    const double theta = exp(lntheta);
    const double theta_min = exp(lnthetamin);
    const double theta_max = exp(lnthetamax);
    log_fatal("theta = %e outside look-up table range [%e, %e]", theta, theta_min, theta_max);
    exit(1);
  }

  const int q = nlsize*nl + N_cgl(ni, nj);
  if(q > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }
  return interpol(table[q], ntheta, lnthetamin, lnthetamax, dlntheta, lntheta, 0, 0);
}

double w_cc_tomo_flatsky(double theta, int nl1, int nl2, int ni, int nj, int limber)
{ // nl{1,2} = lambda_obs bins, n{i,j} = cluster redshift bins
  static cosmopara C;
  static nuisancepara N;
  static double** table;

  const int nlsize = Cluster.N200_Nbin;
  const int nccl_size = tomo.cluster_Nbin; // cross redshift bin not supported so not using
                                           // tomo.cc_clustering_Npowerspectra
  const int NSIZE = nlsize*nlsize*nccl_size;
  const int ntheta = Ntable.N_thetaH;

  const double l_min = w_l_min;
  const double l_max = w_l_max;
  const double loglmax = log(l_max);
  const double loglmin = log(l_min);
  const double dlnl = (loglmax-loglmin)/(1.0*ntheta-1.);
  const double lnrc = 0.5*(loglmax + loglmin);
  const double nc = ntheta/2+1;

  const double logthetamin = (nc-ntheta+1)*dlnl-lnrc;
  const double logthetamax = nc*dlnl-lnrc;
  const double dlntheta = (lnthetamax - lnthetamin)/((double) ntheta);
  const double lntheta = log(theta);

  if (table == 0)
  {
    table = create_double_matrix(0, NSIZE - 1, 0, ntheta - 1); 
  } 
  if (recompute_cc(C, N))
  {
    if (limber == 1)
    {
      // --------------------------------------------------------------------
      // Cocoa: code extracted (& adapted) from xipm_via_hankel (begins)
      // --------------------------------------------------------------------
      typedef fftw_complex fftwZ;

      // go to log-Fourier-space
      fftwZ** flP = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
      for (int j=0; j<NSIZE; j++)
      {
        flP[j] = fftw_malloc((ntheta/2+1)*sizeof(fftwZ));
      }
      { 
        double** lP = (double**) malloc(sizeof(double*)*NSIZE);
        fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
        for (int j=0; j<NSIZE; j++)
        {
          lP[j] = (double*) malloc(ntheta*sizeof(double));
          plan[j] = fftw_plan_dft_r2c_1d(ntheta, lP[j], flP[j], FFTW_ESTIMATE);
        }

        for (int i=0; i<nlsize; i++)  
        { // loop in lambda_obs - (cannot thread loops in lambda_obs  - see cluster_utils)
          for (int j=0; j<nlsize; j++)  
          { // loop in lambda_obs - (cannot thread loops in lambda_obs  - see cluster_utils)
            {
              const int k = 0;
              const int ZCCL1 = k; // cross redshift bin not supported so not using ZCCL1(k)
              const int ZCCL2 = k; // cross redshift bin not supported so not using ZCCL2(k)
              const int q = i*nlsize + j*nlsize*nlsize + k;
              {
                const int p = 0;
                const double l = exp(lnrc + (p - nc)*dlnl);
                lP[q][p] = l*C_cc_tomo_limber(l, i, j, ZCCL1, ZCCL2);
              }
              #pragma omp parallel for
              for(int p=1; p<ntheta; p++)
              {
                const double l = exp(lnrc + (p - nc)*dlnl);
                lP[q][p] = l*C_cc_tomo_limber(l, i, j, k, k);
              }
            }
            #pragma omp parallel for
            for (int k=1; k<nccl_size; k++)  
            {
              const int q = i*nlsize + j*nlsize*nlsize + k;
              const int ZCCL1 = k; // cross redshift bin not supported so not using ZCCL1(k)
              const int ZCCL2 = k; // cross redshift bin not supported so not using ZCCL2(k)
              for(int p=0; p<ntheta; p++)
              {
                const double l = exp(lnrc + (p - nc)*dlnl);
                lP[q][p] = l*C_cc_tomo_limber(l, i, j, ZCCL1, ZCCL2);
              }
            }
          }
        }

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
        lP[j] = (double*) malloc(ntheta*sizeof(double));
        kernel[j] = (fftwZ*) fftw_malloc(sizeof(fftwZ));
        conv[j] = (fftwZ*) fftw_malloc((ntheta/2+1)*sizeof(fftwZ));
        plan[j] = fftw_plan_dft_c2r_1d(ntheta, conv[j], lP[j], FFTW_ESTIMATE);
        tab[j] = create_double_matrix(0, 1, 0, ntheta - 1);
      }
      #pragma omp parallel for
      for (int j=0; j<NSIZE; j++)
      {
        double arg[2];
        arg[0] = 0; // bias
        arg[1] = 0; // order of Bessel function

        for(int i=0; i<(ntheta/2+1); i++)
        {
          const double kk = 2*constants.pi*i/(dlnl*ntheta);
          hankel_kernel_FT(kk, kernel[j], arg, 2);
          conv[j][i][0] = flP[j][i][0]*kernel[j][0][0] - flP[j][i][1]*kernel[j][0][1];
          conv[j][i][1] = flP[j][i][1]*kernel[j][0][0] + flP[j][i][0]*kernel[j][0][1];
        }

        // force Nyquist- and 0-frequency-components to be double
        conv[j][0][1] = 0;
        conv[j][ntheta/2][1] = 0;

        fftw_execute(plan[j]); // Execute FFTW in parallel (thread-safe)

        for(int k=0; k<ntheta; k++)
        {
          const double t = exp((nc-k)*dlnl-lnrc); 
          tab[j][0][ntheta-k-1] = lP[j][k]/(t*2*constants.pi*ntheta);
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
        free_double_matrix(tab[j], 0, 1, 0, ntheta);
      }
      free(flP);
      free(lP);
      free(conv);
      free(kernel);
      free(plan);
      free(tab);
      // ------------------------------------------------------------------------
      // Cocoa: code extracted (& adapted) from twopoint_via_hankel (ends)
      // ------------------------------------------------------------------------
    }
    else
    {
      log_fatal("NonLimber not implemented");
      exit(1);
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }

  if (ni != nj) 
  {
    log_fatal("cross-tomography not supported");
    exit(1);
  }

  if (lntheta < lnthetamin || lntheta > lnthetamax)
  {
    const double theta = exp(lntheta);
    const double theta_min = exp(lnthetamin);
    const double theta_max = exp(lnthetamax);
    log_fatal("theta = %e outside look-up table range [%e, %e]", theta, theta_min, theta_max);
    exit(1);
  }

  const int q = nlsize*nl1 + nlsize*nlsize*nl2 + ni; // cross redshift bin not supported so not
                                                     // using N_CCL(ni, nj) instead of ni
  if(q > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }
  return interpol(table[q], ntheta, lnthetamin, lnthetamax, dlntheta, lntheta, 0, 0);
}

double w_cg_tomo_flatsky(double theta, int nl, int ni, int nj, int limber)
{ // nl = lambda_obs bin, ni = cluster bin, nj = galaxy bin
  static cosmopara C;
  static galpara G;
  static nuisancepara N;
  static double** table;

  const int nlsize = Cluster.N200_Nbin;
  const int nisize = tomo.cluster_Nbin;
  const int njsize = tomo.clustering_Nbin;
  const int NSIZE = nlsize*nisize*njsize;
  const int ntheta = Ntable.N_thetaH;

  const double l_min = w_l_min;
  const double l_max = w_l_max;
  const double loglmax = log(l_max);
  const double loglmin = log(l_min);
  const double dlnl = (loglmax-loglmin)/(1.0*ntheta-1.);
  const double lnrc = 0.5*(loglmax + loglmin);
  const double nc = ntheta/2+1;

  const double logthetamin = (nc-ntheta+1)*dlnl-lnrc;
  const double logthetamax = nc*dlnl-lnrc;
  const double dlntheta = (lnthetamax - lnthetamin)/((double) ntheta);
  const double lntheta = log(theta);

  if (table == 0)
  {
    table  = create_double_matrix(0, NSIZE - 1, 0, ntheta - 1); 
  } 
  if (recompute_cg(C, G, N))
  {
    if (limber == 1)
    {
      typedef fftw_complex fftwZ;

      // go to log-Fourier-space
      fftwZ** flP = (fftwZ**) malloc(sizeof(fftwZ*)*NSIZE);
      for (int j=0; j<NSIZE; j++)
      {
        flP[j] = fftw_malloc((ntheta/2+1)*sizeof(fftwZ));
      }
      { 
        double** lP = (double**) malloc(sizeof(double*)*NSIZE);
        fftw_plan* plan = (fftw_plan*) malloc(sizeof(fftw_plan)*NSIZE);
        for (int j=0; j<NSIZE; j++)
        {
          lP[j] = (double*) malloc(ntheta*sizeof(double));
          plan[j] = fftw_plan_dft_r2c_1d(ntheta, lP[j], flP[j], FFTW_ESTIMATE);
        }

        // Power spectrum on logarithmic bins (begins)
        for (int i=0; i<nlsize; i++)  
        { // cannot thread loop on lambda_obs (see cluster_util)
          {
            const int j = 0;
            {
              const int k = 0;
              const int q = i*nlsize + j*nlsize*nisize + k;
              {
                const int p = 0;
                const double l = exp(lnrc + (p - nc)*dlnl);
                lP[q][p] = l*C_cg_tomo_limber(l, i, j, k);
              }
              #pragma omp parallel for
              for(int p=1; p<ntheta; p++)
              {
                const double l = exp(lnrc + (p - nc)*dlnl);
                lP[q][p] = l*C_cg_tomo_limber(l, i, j, k);
              }
            }
            #pragma omp parallel for
            for (int k=1; k<njsize; k++)  
            {
              const int q = i*nlsize + j*nlsize*nisize + k;
              for(int p=0; p<ntheta; p++)
              {
                const double l = exp(lnrc + (p - nc)*dlnl);
                lP[q][p] = l*C_cg_tomo_limber(l, i, j, k);
              }
            }
          }
          #pragma omp parallel for
          for (int j=1; j<nisize; j++)  
          {
            for (int k=0; k<njsize; k++)  
            {
              const int q = i*nlsize + nlsize*nisize*j + k;
              for(int p=0; p<ntheta; p++)
              {
                const double l = exp(lnrc + (p - nc)*dlnl);
                lP[q][p] = l*C_cg_tomo_limber(l, i, j, k);
              }
            }
          }
        }

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
        lP[j] = (double*) malloc(ntheta*sizeof(double));
        kernel[j] = (fftwZ*) fftw_malloc(sizeof(fftwZ));
        conv[j] = (fftwZ*) fftw_malloc((ntheta/2+1)*sizeof(fftwZ));
        plan[j] = fftw_plan_dft_c2r_1d(ntheta, conv[j], lP[j], FFTW_ESTIMATE);
        tab[j] = create_double_matrix(0, 1, 0, ntheta-1);
      }

      #pragma omp parallel for
      for (int j=0; j<NSIZE; j++)
      {
        double arg[2];
        arg[0] = 0; // bias
        arg[1] = 0; // order of Bessel function

        for(int i=0; i<(ntheta/2+1); i++)
        {
          const double kk = 2*constants.pi*i/(dlnl*ntheta);
          hankel_kernel_FT(kk, kernel[j], arg, 2);
          conv[j][i][0] = flP[j][i][0]*kernel[j][0][0] - flP[j][i][1]*kernel[j][0][1];
          conv[j][i][1] = flP[j][i][1]*kernel[j][0][0] + flP[j][i][0]*kernel[j][0][1];
        }

        // force Nyquist- and 0-frequency-components to be double
        conv[j][0][1] = 0;
        conv[j][ntheta/2][1] = 0;

        fftw_execute(plan[j]); // Execute FFTW in parallel (thread-safe)

        for(int k=0; k<ntheta; k++)
        {
          const double t = exp((nc-k)*dlnl-lnrc); 
          tab[j][0][ntheta-k-1] = lP[j][k]/(t*2*constants.pi*ntheta);
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
        free_double_matrix(tab[j], 0, 1, 0, ntheta);
      }
      free(flP);
      free(lP);
      free(conv);
      free(kernel);
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

  if (lntheta < lnthetamin || lntheta > lnthetamax)
  {
    const double theta = exp(lntheta);
    const double theta_min = exp(lnthetamin);
    const double theta_max = exp(lnthetamax);
    log_fatal("theta = %e outside look-up table range [%e, %e]", theta, theta_min, theta_max);
    exit(1);
  }

  const int q = nlsize*nl + nlsize*nisize*ni + nj;
  if(i > NSIZE-1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }
  return interpol(table[q], ntheta, lnthetamin, lnthetamax, dlntheta, lntheta, 0, 0);
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Limber Approximation (Angular Power Spectrum)
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// cluster lensing 
// ----------------------------------------------------------------------------
// nl = lambda_obs bin, ni = cluster redshift bin, nj = source redshift bin

double int_for_C_cs_tomo_limber_nointerp(double a, void* params)
{ 
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double *) params;

  const int nl = (int) ar[0];
  const int ni = (int) ar[1];
  const int nj = (int) ar[2];
  const double ell = ar[3] + 0.5;
  const int use_linear_ps = (int) ar[4];

  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  
  const double fK = f_K(chidchi.chi);
  const double k = ell / fK;
   
  double res = W_cluster(ni, a, chidchi.chi, hoverh0)*W_kappa(a, fK, nj)
  if (res != 0) 
  {
    res = res*binned_p_cm(k, a, nl, use_linear_ps);
  }
  return res*chidchi.dchida/(fK*fK);
}

double C_cs_tomo_limber_nointerp(double l, int nl, int ni, int nj, int use_linear_ps)
{
  double ar[5] = {(double) nl, (double) ni, (double) nj, l, (double) use_linear_ps};
  
  const double zmin = tomo.cluster_zmin[ni];
  const double zmax = tomo.cluster_zmax[ni];
  const double amin = 1./(1. + zmax);
  const double amax = 1./(1. + zmin);
  
  if(zmin > zmax) 
  {
    return 0.;
  }
  else 
  {
    return int_gsl_integrate_low_precision(int_for_C_cs_tomo_limber_nointerp, ar, amin, amax, 
      NULL, GSL_WORKSPACE_SIZE);
  }
}

double C_cs_tomo_limber(double l, int nl, int ni, int nj)
{
  static cosmopara C;
  static galpara G;
  static nuisancepara N;
  static double** table_Cl = 0;

  const int nell = Ntable.N_ell;
  const int nlsize = Cluster.N200_Nbin;
  const int ngammat_size = tomo.cgl_Npowerspectra;
  const int NSIZE = nlsize*ngammat_size;

  const double lnlmin = log(limits.P_2_s_min);
  const double lnlmax = log(limits.P_2_s_max);
  const double dlnl = (lnlmax - lnlmin)/(nell-1);

  if (table_Cl == 0)
  {
    table_Cl = create_double_matrix(0, NSIZE - 1, 0, nell - 1);     
  }

  if (recompute_cs(C, G, N))
  {
    for(int i=0; i<nlsize; i++)
    { // loop in lambda_obs - (cannot thread loops in lambda_obs  - see cluster_utils)
      {
        const int j = 0;
        const int ZC = ZC(j);
        const int ZS = ZSC(j);
        const int q = nlsize*i + j;
        {
          const int p=0;
          const double lnl = lnlmin + p*dlnl;
          const double l = exp(lnl);
          table_Cl[q][p] = log(C_cs_tomo_limber_nointerp(l, i, ZC, ZS, use_linear_ps_limber));
        }
        #pragma omp parallel for
        for (int p=1; p<nell; p++)
        {
          const double lnl = lnlmin + p*dlnl;
          const double l = exp(lnl);
          table_Cl[q][p] = log(C_cs_tomo_limber_nointerp(l, i, ZC, ZS, use_linear_ps_limber));
        }
      }
      #pragma omp parallel for
      for(int j=1; j<ngammat_size; j++)
      {
        const int ZC = ZC(j);
        const int ZS = ZSC(j);
        const int q = nlsize*i + j;
        for (int p=0; p<nell; p++)
        {
          const double lnl = lnlmin + p*dlnl;
          const double l = exp(lnl);
          table_Cl[q][p] = log(C_cs_tomo_limber_nointerp(l, i, ZC, ZS, use_linear_ps_limber));
        }
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  
  if (!test_zoverlap_c(ni, nj)) 
  {
    return 0.0;
  }

  const double lnl = log(l);
  if (lnl < lnlmin || lnl > lnlmax)
  {
    const double l = exp(log(l));
    log_fatal("l = %e outside look-up table range [%e, %e]", l, exp(lnlmin), exp(lnlmax));
    exit(1);
  }

  const int q = nlsize*nl + N_cgl(ni, nj);
  if(q > NSIZE - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }

  const double f1 = exp(interpol(table_Cl[q], nell, lnlmin, lnlmax, dlnl, lnl, 1, 1));
  if (isnan(f1)) 
  {
    return 0.0;
  }
  else
  {
    return f1;
  }
}

// ---------------------------------------------------------------------------------------------
// Cluster clustering 
// ---------------------------------------------------------------------------------------------
// nl{1,2} = lambda_obs bins, n{i,j} = cluster redshift bins

double int_for_C_cc_tomo_limber_nointerp(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;
   
  const int nl1 = (int) ar[0];
  const int nl2 = (int) ar[1];
  const int ni = (int) ar[2];
  const int nj = (int) ar[3];
  const double ell = ar[4] + 0.5;
  const int use_linear_ps = (int) ar[5];
  
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  
  const double fK = f_K(chidchi.chi);
  const double k  = ell/fK;
  
  double res = W_cluster(ni, a, chidchi.chi, hoverh0)*W_cluster(nj, a, chidchi.chi, hoverh0);
  if (res != 0)
  {
    res = res*binned_p_cc(k, a, nl1, nl2, ni, nj, use_linear_ps);
  }

  return res*chidchi.dchida/(fK*fK);
}

double C_cc_tomo_limber_nointerp(double l, int nl1, int nl2, int ni, int nj, int use_linear_ps)
{ 
  double ar[5] = {(double) nl1, (double) nl2, (double) ni, (double) nj, l, (double) use_linear_ps};  
  const double zmin = fmax(tomo.cluster_zmin[ni], tomo.cluster_zmin[nj]);
  const double zmax = fmin(tomo.cluster_zmax[ni], tomo.cluster_zmax[nj]);
  const double amin = 1./(1. + zmax);
  const double amax = 1./(1. + zmin);

  if(zmin > zmax) 
  {
    return 0.;
  }
  else
  {
    return int_gsl_integrate_low_precision(int_for_C_cc_tomo_limber_nointerp, (void*) ar, 
      amin, amax, NULL, GSL_WORKSPACE_SIZE);
  }
}

double C_cc_tomo_limber(double l, int nl1, int nl2, int ni, int nj) 
{
  static cosmopara C;
  static nuisancepara N;
  static double** table_Cl = 0;

  if (ni != nj)
  {
    log_fatal("ni != nj tomography not supported");
    exit(1);
  }  

  const int nell = Ntable.N_ell;
  const int nlsize = Cluster.N200_Nbin;
  const int nccl_size = tomo.cluster_Nbin; // cross redshift bin not supported so not using
                                           // tomo.cc_clustering_Npowerspectra
  const int NSIZE = nccl_size*nlsize*nlsize;
  
  const double loglmin = log(limits.P_2_s_min);
  const double loglmax = log(limits.P_2_s_max);
  const double dl = (loglmax - loglmin)/(nell - 1); 

  if (table_Cl == 0)
  { 
    table_Cl = create_double_matrix(0, NSIZE - 1, 0, nell);     
  }

  if (recompute_cc(C, N))
  {
    for (int i=0; j<nlsize; i++) 
    { // loop in lambda_obs - (cannot thread loops in lambda_obs  - see cluster_utils)
      for (int j=0; j<i+1; j++) 
      { // loop in lambda_obs - (cannot thread loops in lambda_obs  - see cluster_utils)
        {
          const int k = 0;     
          const int ZCCL1 = k; // cross redshift bin not supported so not using ZCCL1(k)
          const int ZCCL2 = k; // cross redshift bin not supported so not using ZCCL2(k)
          const int q = nlsize*i + nlsize*nlsize*j + k;
          const int qstar = nlsize*j + nlsize*nlsize*i + k;
          {
            const int p = 0;
            const double lnl = loglminn + p*dl;
            const double l = exp(lnl); 
            table_Cl[q][p] = log(C_cc_tomo_limber_nointerp(l, i, j, ZCCL1, ZCCL2, use_linear_ps_limber));
            table_Cl[qstar][p] = table_Cl[qstar][p];
          }
          #pragma omp parallel for
          for (int p=1; p<nell; ++p)
          {
            const double lnl = loglminn + p*dl;
            const double l = exp(lnl); 
            table_Cl[q][p] = log(C_cc_tomo_limber_nointerp(l, i, j, ZCCL1, ZCCL2, use_linear_ps_limber));
            table_Cl[qstar][p] = table_Cl[q][p];
          }
        }
        #pragma omp parallel for
        for (int k=1; k<nccl_size; j++) // loop in ni (= nj)
        {
          const int ZCCL1 = k; // cross redshift bin not supported so not using ZCCL1(k)
          const int ZCCL2 = k; // cross redshift bin not supported so not using ZCCL2(k)
          const int q = nlsize*i + nlsize*nlsize*j + k;
          const int qstar = nlsize*j + nlsize*nlsize*i + k;
          for (int p=0; p<nell; ++p)
          {
            const double lnl = loglminn + p*dl;
            const double l = exp(lnl); 
            table_Cl[q][p] = log(C_cc_tomo_limber_nointerp(l, i, j, ZCCL1, ZCCL2, use_linear_ps_limber));
            table_Cl[qstar][p] = table_Cl[q][p];
          }
        }
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

  const int q = nlsize*nl1 + nlsize*nlsize*nl2 + ni; // cross redshift bin not supported so not
                                                     // using N_CCL(ni, nj) instead of ni
  if(q > NSIZE-1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }

  const double f1 = exp(interpol(table_Cl[q], nell, loglmin, loglmax, dl, lnl, 1, 1));
  if (isnan(f1)) 
  {
    return 0.0
  }
  return f1;
}


// ---------------------------------------------------------------------------------------------
// cluster x galaxy clustering
// ---------------------------------------------------------------------------------------------
// nl = lambda_obs bin, ni = cluster redshift bin, nj = galaxy redshift bin

double int_for_C_cg_tomo_limber_nointerp(double a, void* params)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;

  const int nl = (int) ar[0];
  const int ni = (int) ar[1];
  const int nj = (int) ar[2];
  const double ell = ar[3] + 0.5;
  const int use_linear_ps = (int) ar[4]

  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  
  const  double fK  = f_K(chidchi.chi);
  const double k = ell/fK;

  double res = W_cluster(ni, a, chidchi.chi, hoverh0) * W_HOD(a, nj, hoverh0)

  if (res != 0)
  {
    res = res*binned_p_cg(k, a, nl1, ni, nj, use_linear_ps);
  }
  return res*chidchi.dchida/(fK*fK);
}

double C_cg_tomo_limber_nointerp(double l, int nl, int ni, int nj, int use_linear_ps) 
{
   double ar[5] = {(double) nl, (double) ni, (double) nj, l, (double) use_linear_ps};
   const double zmin = fmax(tomo.cluster_zmin[ni], tomo.clustering_zmin[nj]);
   const double zmax = fmin(tomo.cluster_zmax[ni], tomo.clustering_zmax[nj]);
   const double amin = 1./(1. + zmax);
   const double amax = 1./(1 + zmin);
   
   if (zmin > zmax) 
   {
      return 0.;
   }
   else
   {
      return int_gsl_integrate_low_precision(int_for_C_cg_tomo_limber_nointerp, (void*) ar, 
        amin, amax, NULL, GSL_WORKSPACE_SIZE);
   }
}

double C_cg_tomo_limber(double l, int nl, int ni, int nj)
{
  static cosmopara C;
  static galpara G;
  static nuisancepara N;
  static double** table_Cl = 0;

  const int nell = Ntable.N_ell;
  const int nlsize = Cluster.N200_Nbin;
  const int nisize = tomo.cluster_Nbin;
  const int njsize = tomo.clustering_Nbin;
  const int NSIZE = nlsize*nisize*njsize;

  const double loglmin = log(limits.P_2_s_min);
  const double loglmax = log(limits.P_2_s_max);
  const double dl = (loglmax - loglmin)/(nell-1);

  if (table_Cl == 0)
  { 
    table_Cl = create_double_matrix(0, NSIZE - 1, 0, nell);      
  }
  if (recompute_cg(C, N))
  {
    for (int i=0; i<nlsize; i++) 
    { // loop in lambda_obs - (cannot thread loops in lambda_obs  - see cluster_utils)
      {
        const int j = 0;
        {
          const int k = 0;
          const int q = i*nlsize + nlsize*nisize*j + k;
          {
            const int p = 0;
            const double lnl = loglmin + p*dl;
            const double l = exp(lnl);
            table_Cl[q][p] = log(C_cg_tomo_limber_nointerp(l, i, j, k, use_linear_ps_limber));
          }
          #pragma omp parallel for
          for (int p=1; p<nell; ++p)
          {
            const double lnl = loglmin + p*dl;
            const double l = exp(lnl);
            table_Cl[q][p] = log(C_cg_tomo_limber_nointerp(l, i, j, k, use_linear_ps_limber));
          }
        }
        #pragma omp parallel for
        for (int k=1; k<njsize; k++)  
        {
          const int q = i*nlsize + nlsize*nisize*j + k;
          for (int p=0; p<nell; ++p)
          {
            const double lnl = loglmin + p*dl;
            const double l = exp(lnl);
            table_Cl[q][p] = log(C_cg_tomo_limber_nointerp(l, i, j, k, use_linear_ps_limber));
          }
        }
      }
      #pragma omp parallel for
      for (int j=1; j<nisize; j++)  
      {
        for (int k=0; k<njsize; k++)  
        {
          const int q = i*nlsize + nlsize*nisize*j + k;
          for (int p=0; p<nell; ++p)
          {
            const double lnl = loglmin + p*dl;
            const double l = exp(lnl);
            table_Cl[q][p] = log(C_cg_tomo_limber_nointerp(l, i, j, k, use_linear_ps_limber));
          }
        }
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }
   
  const double lnl = log(l);
  if (lnl < lnlmin || lnl > lnlmax)
  {
    const double l = exp(l);
    log_fatal("l = %e outside look-up table range [%e, %e]", l, exp(lnlmin), exp(lnlmax));
    exit(1);
  } 

  const int q = nlsize*nl + nlsize*nisize*ni + nj;
  if(q > NSIZE-1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }
  const double f1 = exp(interpol(table_Cl[q], nell, loglmin, loglmax, dl, lnl, 1, 1));
  if (isnan(f1)) 
  {
    return 0.0;
  }
  return f1;
}


// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// cluster number counts
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// nl = lambda_obs bin, ni = cluster redshift bin

double binned_average_number_counts(int nl, double z)
{ // def: eq 3 of https://arxiv.org/pdf/1810.09456.pdf; nl = lambda_obs bin, nz = redshift bin
  double param[2] = {(double) nl, z};
  const double mmin = 12.0/0.4342944819; // log(pow(10.,12.)) --- 0.4342944819 = log10(e)
  const double mmax = 15.9/0.4342944819; // log(pow(10.,15.9)) --- 0.4342944819 = log10(e)
  
  return int_gsl_integrate_low_precision(int_dndlogM_times_binned_P_lambda_obs_given_M, 
    (void*) param, mmin, mmax, NULL, GSL_WORKSPACE_SIZE);
}

double int_projected_average_number_counts(double a, void* params)
{
  double* ar = (double*) params;   

  if(!(a>0)) 
  {
    log_fatal("a > 0 not true");
    exit(1);
  }
  const double z = 1./a-1 ;
  const int nl = (int) ar[0];
  const int interpolate_survey_area = (int) ar[1];
  const double norm = interpolate_survey_area > 0 ? get_area(z) : survey.area;  
  struct chis chidchi = chi_all(a);  
  const double fK = f_K(chidchi.chi);
  
  return fK*fK*chidchi.dchida*binned_average_number_counts(nl, z)*norm;
}

double projected_average_number_counts(int nl, int ni)
{ // nl = lambda_obs bin, ni = cluster redshift bin
  const int interpolate_survey_area = survey.area < 0 ? 1 : 0;
  double params[2] = {(double) nl, interpolate_survey_area};

  const double amin = 1./(1 + tomo.cluster_zmax[ni]);
  const double amax = 1./(1 + tomo.cluster_zmin[ni]);
  const double tmp = (4.0*M_PI/41253.0);

  return tmp*int_gsl_integrate_low_precision(int_projected_average_number_counts, 
    (void*) params, amin, amax, NULL, GSL_WORKSPACE_SIZE);
}

