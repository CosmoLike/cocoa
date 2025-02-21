#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_spline.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../cfftlog/cfftlog.h"

#include "bias.h"
#include "basics.h"
#include "cosmo3D.h"
#include "cluster_util.h"
#include "cosmo2D.h"
#include "cosmo2D_cluster.h"
#include "radial_weights.h"
#include "recompute.h"
#include "redshift_spline.h"
#include "structs.h"

#include "log.c/src/log.h"

static int GSL_WORKSPACE_SIZE = 250;
static int use_linear_ps_limber = 0; /* 0 or 1 */
static int INCLUDE_MAG_IN_C_CC_NONLIMBER = 0; /* 0 or 1 */
static int INCLUDE_MAG_IN_C_CG_NONLIMBER = 0; /* 0 or 1 */

// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
// Correlation Functions (real Space) - Full Sky - bin average
// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------

double w_gammat_cluster_tomo(const int nt, const int nl, const int ni, const int nj, 
const int limber)
{ // nt = theta bin, nl = lambda_obs bin, ni = cluster redshift bin, nj = source redshift bin
  if (like.Ntheta == 0)
  {
    log_fatal("like.Ntheta not initialized");
    exit(1);
  }
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double** Pl = 0;
  static double* w_vec = 0;

  const int nell = limits.LMAX;
  const int ntheta = like.Ntheta;
  const int nlsize = Cluster.N200_Nbin;
  const int ngammat_size = tomo.cgl_Npowerspectra;
  const int NSIZE = nlsize*ngammat_size;

  if (Pl == 0)
  {    
    Pl = (double**) malloc(sizeof(double*)*ntheta);
    for (int i=0; i<ntheta; i++)
    {
      Pl[i] = (double*) malloc(sizeof(double)*nell);
    }
    w_vec = (double*) malloc(sizeof(double)*NSIZE*ntheta); 
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
      double* Pmin = (double*) malloc(sizeof(double)*(nell + 1));
      double* Pmax = (double*) malloc(sizeof(double)*(nell + 1));
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
  if (recompute_cs(C, G, N))
  {
    double** Cl = malloc(NSIZE*sizeof(double*));
    for (int i=0; i<NSIZE; i++)
    {
      Cl[i] = calloc(nell, sizeof(double));
    }
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
    {
      const int ZC = ZCL(0);
      const int ZSC = ZCS(0);
      double init_static_vars_only = C_cs_tomo_limber(limits.LMIN_tab + 1, 0, ZC, ZSC);
    }
    #pragma GCC diagnostic pop
    #pragma GCC diagnostic pop
    if (limber == 1)
    { 
      #pragma omp parallel for collapse(3) 
      for (int i=0; i<nlsize; i++)
      {  
        for (int j=0; j<ngammat_size; j++)
        { 
          for (int l=2; l<nell; l++)
          {
            const int ZC = ZCL(j);
            const int ZSC = ZCS(j);
            const int q = i*ngammat_size + j;
            Cl[q][l] = (l > limits.LMIN_tab) ? C_cs_tomo_limber(l, i, ZC, ZSC) :
              C_cs_tomo_limber_nointerp(l, i, ZC, ZSC, use_linear_ps_limber, 0);
          }
        }    
      } 
    }
    else
    { 
      log_fatal("NonLimber not implemented");
      exit(1);     
    }
    #pragma omp parallel for collapse(3)
    for (int i=0; i<nlsize; i++)
    {
      for (int j=0; j<ngammat_size; j++)
      { 
        for (int p=0; p<ntheta; p++)
        {
          const int nz = i*ngammat_size + j;
          const int q  = nz*ntheta + p;
          w_vec[q] = 0;
          for (int l=1; l<nell; l++)
          {
            w_vec[q] += Pl[p][l]*Cl[nz][l];
          }
        }
      }
    }
    for (int i=0; i<NSIZE; i++)
    {
      free(Cl[i]);
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
  if (nl < 0 || nl > nlsize - 1)
  {
    log_fatal("error in selecting bin number nl = %d (max %d)", nl, nlsize);
    exit(1); 
  } 
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d (max %d)", ni, redshift.clustering_nbin);
    exit(1); 
  } 
  if (nj < 0 || nj > redshift.shear_nbin - 1)
  {
    log_fatal("error in selecting bin number nj = %d (max %d)", nj, redshift.shear_nbin);
    exit(1); 
  }
  const int ntomo = N_cgl(ni, nj);
  if (!(ntomo>0))
  {
    return 0.0;
  }
  {
    const int q = (nl*ngammat_size + ntomo)*ntheta + nt;
    if (q > NSIZE*ntheta - 1)
    {
      log_fatal("internal logic error in selecting bin number");
      exit(1);
    }
    return w_vec[q];
  }
}


double w_cc_tomo(const int nt, const int nl1, const int nl2, const int ni, const int nj, 
const int limber)
{ // nt = theta bin , nl{1,2} = lambda_obs bins, n{i,j} = cluster redshift bins
  if (like.Ntheta == 0)
  {
    log_fatal("like.Ntheta not initialized");
    exit(1);
  }
  static cosmopara C;
  static nuisancepara N;
  static double** Pl = 0;
  static double* w_vec = 0;

  const int nlsize = Cluster.N200_Nbin;
  const int nccl_size = redshift.clustering_nbin; // cross redshift bin not supported so not using
                                           // tomo.cc_clustering_Npowerspectra

  const int nell = limits.LMAX;
  const int ntheta = like.Ntheta;
  const int NSIZE = nlsize*nlsize*nccl_size;

  if (Pl == 0)
  {
    Pl = (double**) malloc(sizeof(double*)*ntheta);
    for (int i=0; i<ntheta; i++)
    {
      Pl[i] = (double*) malloc(sizeof(double)*nell);
    }
    w_vec = (double*) malloc(sizeof(double)*NSIZE*ntheta);
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
      Pmin[i] = (double*) malloc(sizeof(double)*(nell + 1));
      Pmax[i] = (double*) malloc(sizeof(double)*(nell + 1));
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
      const double tmp = (1.0 / (xmin[i] - xmax[i]))*(1.0 / (4.0 * M_PI));
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
  if (recompute_cc(C, N))
  {
    double** Cl = malloc(NSIZE*sizeof(double*));
    for (int i=0; i<NSIZE; i++)
    {
      Cl[i] = calloc(nell, sizeof(double));  
    }
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
    {
      const int i = 0;
      const int j = 0;
      const int ZCCL1 = 0;
      const int ZCCL2 = 0;
      double init_static_vars_only = C_cc_tomo_limber(limits.LMIN_tab + 1, i, j, ZCCL1, ZCCL2);
    }
    #pragma GCC diagnostic pop
    #pragma GCC diagnostic pop
    if (limber == 1)
    { 
      #pragma omp parallel for collapse(3)
      for (int i=0; i<nlsize; i++) 
      {
        for (int k=0; k<nccl_size; k++)
        {
          for (int l=1; l<nell; l++)
          {
            for (int j=i; j<nlsize; j++)
            {
              const int q     = i*nlsize*nccl_size + j*nccl_size + k;
              const int qstar = j*nlsize*nccl_size + i*nccl_size + k;
              const int ZCCL1 = k; // cross redshift bin not supported so not using ZCCL1(k)
              const int ZCCL2 = k; // cross redshift bin not supported so not using ZCCL2(k)
              Cl[q][l] = (l > limits.LMIN_tab) ? C_cc_tomo_limber(l, i, j, ZCCL1, ZCCL2) :
                C_cc_tomo_limber_nointerp(l, i, j, ZCCL1, ZCCL2, use_linear_ps_limber, 0);
              Cl[qstar][l] = Cl[q][l];
            }
          } 
        }
      }
    }
    else
    {
      for (int i=0; i<nlsize; i++)  // NON LIMBER PART
      { 
        for (int j=i; j<nlsize; j++)  
        {
          for (int k=0; k<nccl_size; k++) // can't openmp - C_cc_tomo - FFTW calls
          {
            const int L = 1;
            const double tol = 0.01;      // required fractional accuracy in C(l)
            const double dev = 10.0 * tol;// will be diff exact vs Limber init to large value in
                                          // order to start while loop
            const int ZCCL1 = k;          // cross redshift bin not supported so not using ZCCL1(k)
            const int ZCCL2 = k;          // cross redshift bin not supported so not using ZCCL2(k)
            const int q = i*nlsize*nccl_size + j*nccl_size + k;
            C_cc_tomo(L, i, j, ZCCL1, ZCCL2, Cl[q], dev, tol); // can't openmp
          }
        }
      }
      #pragma omp parallel for collapse(3)
      for (int i=0; i<nlsize; i++) // LIMBER PART
      {
        for (int k=0; k<nccl_size; k++)
        {
          for (int l=limits.LMAX_NOLIMBER+1; l<nell; l++)
          {
            for (int j=i; j<nlsize; j++)
            {
              const int q = i*nlsize*nccl_size + j*nccl_size + k;
              const int ZCCL1 = k; // cross redshift bin not supported so not using ZCCL1(k)
              const int ZCCL2 = k; // cross redshift bin not supported so not using ZCCL2(k)
              Cl[q][l] = C_cc_tomo_limber(l, i, j, ZCCL1, ZCCL2);
            }
          }
        }
      }
      #pragma omp parallel for collapse(3)
      for (int i=0; i<nlsize; i++)
      {
        for (int k=0; k<nccl_size; k++)
        {
          for (int l=1; l<nell; l++)
          {
            for (int j=i; j<nlsize; j++)
            {
              const int q     = i*nlsize*nccl_size + j*nccl_size + k;
              const int qstar = j*nlsize*nccl_size + i*nccl_size + k;
              Cl[qstar][l] = Cl[q][l];
            }
          }
        }
      } 
    }
    #pragma omp parallel for collapse(4)
    for (int i=0; i<nlsize; i++)
    { 
      for (int j=0; j<nlsize; j++)  
      {
        for (int k=0; k<nccl_size; k++)  
        {
          for (int p=0; p<ntheta; p++)
          {
            const int nz = i*nlsize*nccl_size + j*nccl_size + k;
            const int q = nz*ntheta + p;
            w_vec[q] = 0;
            for (int l=1; l<nell; l++)
            {
              w_vec[q] += Pl[p][l]*Cl[nz][l];
            }
          }
        }
      }
    }
    for (int i=0; i<NSIZE; i++)
    {
      free(Cl[i]);
    }
    free(Cl);
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  if (ni != nj)
  {
    log_fatal("ni != nj cross tomography not supported");
    exit(1);
  } 
  if (nt < 0 || nt > like.Ntheta - 1)
  {
    log_fatal("error in selecting bin number nt = %d (max %d)", nt, like.Ntheta);
    exit(1); 
  }
  if (nl1 < 0 || nl1 > nlsize - 1)
  {
    log_fatal("error in selecting bin number nl1 = %d (max %d)", nl1, nlsize);
    exit(1); 
  } 
  if (nl2 < 0 || nl2 > nlsize - 1)
  {
    log_fatal("error in selecting bin number nl2 = %d (max %d)", nl2, nlsize);
    exit(1); 
  } 
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d (max %d)", ni, redshift.clustering_nbin);
    exit(1); 
  } 
  if (nj < 0 || nj > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number nj = %d (max %d)", nj, redshift.clustering_nbin);
    exit(1); 
  }
  const int q = (nl1*nlsize*nccl_size + nl2*nccl_size + ni)*ntheta + nt; // cross redshift bin not 
                                                                         // supported so not using 
                                                                         // N_CCL(ni, nj)
  if (q > NSIZE*ntheta - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  }
  return w_vec[q];
}

double w_cg_tomo(const int nt, const int nl, const int ni, const int nj, const int limber)
{
  if (like.Ntheta == 0)
  {
    log_fatal("like.Ntheta not initialized");
    exit(1);
  }
  static cosmopara C;
  static galpara G;
  static nuisancepara N;
  static double** Pl = 0;
  static double* w_vec = 0;

  const int ntheta = like.Ntheta;
  const int nlsize = Cluster.N200_Nbin;
  const int ncg_size = tomo.cg_clustering_Npowerspectra;
  const int nell = limits.LMAX;
  const int NSIZE = nlsize*ncg_size;

  if (Pl == 0)
  {
    Pl = (double**) malloc(sizeof(double*)*ntheta);
    for (int i=0; i<ntheta; i++)
    {
      Pl[i] = (double*) malloc(sizeof(double)*nell);
    }
    w_vec = (double*) malloc(sizeof(double)*NSIZE*ntheta);
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
      Pmin[i] = (double*) malloc(sizeof(double)*(nell + 1));
      Pmax[i] = (double*) malloc(sizeof(double)*(nell + 1));
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
      const double tmp = (1.0/(xmin[i] - xmax[i]))*(1. / (4. * M_PI));
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
  if (recompute_cg(C, G, N))
  {
    double** Cl = malloc(NSIZE*sizeof(double*));
    for (int i=0; i<NSIZE; i++)
    {
      Cl[i] = calloc(nell, sizeof(double));
    }
    if (limber == 1)
    {
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wunused-variable"
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
      {
        const double Z1 = ZCGCL1(0);
        const double Z2 = ZCGCL2(0);
        double init_static_vars_only = C_cg_tomo_limber(limits.LMIN_tab + 1, 0, Z1, Z2);
      }
      #pragma GCC diagnostic pop
      #pragma GCC diagnostic pop
      #pragma omp parallel for collapse(3)
      for (int i=0; i<nlsize; i++)
      {
        for (int j=0; j<ncg_size; j++)
        {
          for (int l=1; l<nell; l++)
          {
            const int q = i*ncg_size+ j;
            const double Z1 = ZCGCL1(j);
            const double Z2 = ZCGCL2(j);
            Cl[q][l] = (l > limits.LMIN_tab) ? C_cg_tomo_limber(l, i, Z1, Z2) :
              C_cg_tomo_limber_nointerp(l, i, Z1, Z2, use_linear_ps_limber, 0);
          }
        }
      }
    }
    else
    {
      log_fatal("NonLimber not implemented");
      exit(1);
    }
    #pragma omp parallel for collapse(3)
    for (int i=0; i<nlsize; i++)
    {
      for (int j=0; j<ncg_size; j++)
      {
        for (int p=0; p<ntheta; p++)
        {
          const int nz = i*ncg_size + j;
          const int q  = nz*ntheta + p;
          w_vec[q] = 0;
          for (int l=1; l<nell; l++)
          {
            w_vec[q] += Pl[p][l]*Cl[nz][l];
          }
        }
      }
    }
    for (int i=0; i<NSIZE; i++)
    {
      free(Cl[i]);
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
  if (nl < 0 || nl > nlsize - 1)
  {
    log_fatal("error in selecting bin number nl = %d (max %d)", nl, nlsize);
    exit(1); 
  } 
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d (max %d)", ni, redshift.clustering_nbin);
    exit(1); 
  } 
  if (nj < 0 || nj > tomo.clustering_Nbin - 1)
  {
    log_fatal("error in selecting bin number nj = %d (max %d)", nj, tomo.clustering_Nbin);
    exit(1); 
  }
  const int ntomo = N_CGCL(ni, nj);
  if (!(ntomo>0))
  {
    return 0.0;
  }
  else
  {
    const int q = nl*ncg_size + ntomo;
    if (q > NSIZE-1)
    {
      log_fatal("internal logic error in selecting bin number");
      exit(1);
    }
    return w_vec[q];
  }
}

// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
// Correlation Functions (real space) - flat sky
// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------

double w_gammat_cluster_tomo_flatsky(const double theta, const int nl, const int ni, const int nj,
const int limber)
{ // nl = lambda_obs bin, ni = cluster redshift bin, nj = source redshift bin
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double** table;
 
  const int ntheta = Ntable.N_thetaH;
  const int nlsize = Cluster.N200_Nbin;
  const int ngammat_size = tomo.cgl_Npowerspectra;
  const int NSIZE = nlsize*ngammat_size;
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
  if (recompute_cs(C, G, N))
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

      // ------------------------------------------------------------------------
      // Power spectrum on logarithmic bins begins
      // ------------------------------------------------------------------------
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wunused-variable"
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
      {
        const int i = 0;
        const int j = 0;
        const int ZC = ZCL(j);
        const int ZSC = ZCS(j);
        double init_static_vars_only = C_cs_tomo_limber(limits.LMIN_tab + 1, i, ZC, ZSC);
      }
      #pragma GCC diagnostic pop
      #pragma GCC diagnostic pop
      if (limber == 1)
      { 
        #pragma omp parallel for collapse (3)
        for (int i=0; i<nlsize; i++)
        {
          for (int j=0; j<ngammat_size; j++)
          { 
            for (int p=0; p<ntheta; p++)
            {
              const int ZC = ZCL(j);
              const int ZSC = ZCS(j);
              const int q = i*ngammat_size + j;
              const double l = exp(lnrc + (p - nc)*dlnl);
              lP[q][p] = (l > limits.LMIN_tab) ? l*C_cs_tomo_limber(l, i, ZC, ZSC) :
                  l*C_cs_tomo_limber_nointerp(l, i, ZC, ZSC, use_linear_ps_limber, 0);
            }
          } 
        }
      }
      else
      {
        log_fatal("NonLimber not implemented");
        exit(1);
      }
      // ------------------------------------------------------------------------
      // Power spectrum on logarithmic ends
      // ------------------------------------------------------------------------
      
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
      tab[j] = (double**) malloc(sizeof(double*)*1);
      for (int i=0; i<1; i++)
      {
        tab[j][i] = (double*) malloc(sizeof(double)*ntheta);
      }
    }
    #pragma omp parallel for
    for (int j=0; j<NSIZE; j++)
    { 
      double arg[2];
      arg[0] = 0; // bias
      arg[1] = 2; // order of Bessel function - J2

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
        const double t = exp((nc-k)*dlnl-lnrc); 
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
      free(tab[j]);
    }
    free(flP);
    free(lP);
    free(conv);
    free(kernel);
    free(plan);
    free(tab);
    update_galpara(&G);
    update_cosmopara(&C);
    update_nuisance(&N);
  } 
  if (nl < 0 || nl > nlsize - 1)
  {
    log_fatal("error in selecting bin number nl = %d (max %d)", nl, nlsize);
    exit(1); 
  } 
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d (max %d)", ni, redshift.clustering_nbin);
    exit(1); 
  } 
  if (nj < 0 || nj > redshift.shear_nbin - 1)
  {
    log_fatal("error in selecting bin number nj = %d (max %d)", nj, redshift.shear_nbin);
    exit(1); 
  }
  const double lntheta = log(theta);
  if (lntheta < lnthetamin || lntheta > lnthetamax)
  {
    const double theta = exp(lntheta);
    const double theta_min = exp(lnthetamin);
    const double theta_max = exp(lnthetamax);
    log_fatal("theta = %e outside look-up table range [%e,%e]", theta, theta_min, theta_max);
    exit(1);
  }
  const int ntomo = N_cgl(ni, nj);
  if (!(ntomo>0))
  {
    return 0.0;
  }
  else
  {
    const int q = nl*ngammat_size + ntomo;
    if (q > NSIZE - 1)
    {
      log_fatal("internal logic error in selecting bin number");
      exit(1);
    }
    return interpol(table[q], ntheta, lnthetamin, lnthetamax, dlntheta, lntheta, 0, 0);
  }
}

double w_cc_tomo_flatsky(const double theta, const int nl1, const int nl2, const int ni, 
const int nj, const int limber)
{ // nl{1,2} = lambda_obs bins, n{i,j} = cluster redshift bins
  static cosmopara C;
  static nuisancepara N;
  static double** table;

  const int nlsize = Cluster.N200_Nbin;
  const int nccl_size = redshift.clustering_nbin; // cross redshift bin not supported so not using
                                           // tomo.cc_clustering_Npowerspectra
  const int NSIZE = nlsize*nlsize*nccl_size;
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
  if (recompute_cc(C, N))
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

      // ------------------------------------------------------------------------
      // Power spectrum on logarithmic bins begins
      // ------------------------------------------------------------------------
      double* ll = (limber == 1) ? NULL : calloc(limits.LMAX_NOLIMBER, sizeof(double));
      if (limber != 1)
      { 
        for (int i=0; i<limits.LMAX_NOLIMBER; i++)
        {
          ll[i] = i;
        }
      }
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wunused-variable"
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
      {
        const int i = 0;
        const int j = 0;
        double init_static_vars_only = C_cc_tomo_limber(limits.LMIN_tab + 1, i, j, 0, 0);
      }
      #pragma GCC diagnostic pop
      #pragma GCC diagnostic pop
      for (int i=0; i<nlsize; i++) 
      { 
        gsl_spline*** fCL_NL = 
          (limber == 1) ? NULL : (gsl_spline***) malloc(sizeof(gsl_spline**)*nlsize);
        double*** Cl_NL = (double***) malloc(sizeof(double**)*nlsize);
        if (limber != 1)
        { 
          for (int j=i; j<nlsize; j++)  
          { 
            fCL_NL[j] = (gsl_spline**) malloc(sizeof(gsl_spline*)*nccl_size);
            Cl_NL[j] = (double**) malloc(sizeof(double*)*nccl_size);
            for (int k=0; k<nccl_size; k++)
            {
              Cl_NL[j][k] = calloc(limits.LMAX_NOLIMBER, sizeof(double));
            }
            const int L = 1;
            const double tol = 0.0075;        // required fractional accuracy in C(l)
            const double dev = 10. * tol;     // will be diff exact vs Limber init to large
                                              // value in order to start while loop
            for (int k=0; k<nccl_size; k++)
            {
              const int ZCCL1 = k; // cross redshift bin not supported so not using ZCCL1(k)
              const int ZCCL2 = k; // cross redshift bin not supported so not using ZCCL1(k)
              C_cc_tomo(L, i, j, ZCCL1, ZCCL2, Cl_NL[j][k], dev, tol);
              
              const gsl_interp_type* T = gsl_interp_linear;
              fCL_NL[j][k] = gsl_spline_alloc(T, limits.LMAX_NOLIMBER);
              if (fCL_NL[j][k] == NULL)
              {
                log_fatal("fail allocation");
                exit(1);
              }
            }
          }
          #pragma omp parallel for collapse(2)
          for (int j=i; j<nlsize; j++)  
          {
            for (int k=0; k<nccl_size; k++)
            {
              int status = gsl_spline_init(fCL_NL[j][k], ll, Cl_NL[j][k], limits.LMAX_NOLIMBER);
              if (status) 
              {
                log_fatal(gsl_strerror(status));
                exit(1);
              }
            }
          }
        }    
        #pragma omp parallel for collapse(3)
        for (int j=i; j<nlsize; j++)  
        {
          for (int k=0; k<nccl_size; k++)  
          {
            for (int p=0; p<ntheta; p++)
            {
              const int q     = i*nlsize*nccl_size + j*nccl_size + k;
              const int qstar = j*nlsize*nccl_size + i*nccl_size + k;
              const int ZCCL1 = k; // cross redshift bin not supported so not using ZCCL1(k)
              const int ZCCL2 = k; // cross redshift bin not supported so not using ZCCL2(k)
              const double l = exp(lnrc + (p - nc)*dlnl);
              if (limber == 1 || (limber != 1 && l > limits.LMAX_NOLIMBER - 1))
              {
                lP[q][p] = (l > limits.LMIN_tab) ? l*C_cc_tomo_limber(l, i, j, ZCCL1, ZCCL2) :
                  l*C_cc_tomo_limber_nointerp(l, i, j, ZCCL1, ZCCL2, use_linear_ps_limber, 0); 
              }
              else
              {
                double CLNL;
                int status = gsl_spline_eval_e(fCL_NL[j][k], l, NULL, &CLNL);
                if (status) 
                {
                  log_fatal(gsl_strerror(status));
                  exit(1);
                }
                lP[q][p] = l*CLNL;
              }
              lP[qstar][p] = lP[q][p];
            }
          }
        }
        if (limber != 1)
        {
          for (int j=i; j<nlsize; j++) 
          {
            for (int k=0; k<nccl_size; k++)
            {
              free(Cl_NL[j][k]);
              gsl_spline_free(fCL_NL[j][k]);
            }
          }
          for (int j=i; j<nlsize; j++)
          {
            free(Cl_NL[j]);
            free(fCL_NL[j]);
          }
          free(Cl_NL);
          free(fCL_NL);
        }
      } 
      if (limber != 1)
      { 
        free(ll);
      }
      // ------------------------------------------------------------------------
      // Power spectrum on logarithmic bins ends
      // ------------------------------------------------------------------------
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
      tab[j] = (double**) malloc(sizeof(double*)*1);
      for (int i=0; i<1; i++)
      {
        tab[j][i] = (double*) malloc(sizeof(double)*ntheta);
      }
    }
    #pragma omp parallel for
    for (int j=0; j<NSIZE; j++)
    {
      double arg[2];
      arg[0] = 0; // bias
      arg[1] = 0; // order of Bessel function - J0

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
        const double t = exp((nc - k)*dlnl - lnrc); 
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
      free(tab[j]);
    }
    free(flP);
    free(lP);
    free(conv);
    free(kernel);
    free(plan);
    free(tab);
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  if (ni != nj) 
  {
    log_fatal("cross-tomography not supported");
    exit(1);
  }
  if (nl1 < 0 || nl1 > nlsize - 1)
  {
    log_fatal("error in selecting bin number nl1 = %d (max %d)", nl1, nlsize);
    exit(1); 
  } 
  if (nl2 < 0 || nl2 > nlsize - 1)
  {
    log_fatal("error in selecting bin number nl2 = %d (max %d)", nl2, nlsize);
    exit(1); 
  } 
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d (max %d)", ni, redshift.clustering_nbin);
    exit(1); 
  } 
  if (nj < 0 || nj > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number nj = %d (max %d)", nj, redshift.clustering_nbin);
    exit(1); 
  }
  const double lntheta = log(theta);
  if (lntheta < lnthetamin || lntheta > lnthetamax)
  {
    const double theta = exp(lntheta);
    const double theta_min = exp(lnthetamin);
    const double theta_max = exp(lnthetamax);
    log_fatal("theta = %e outside look-up table range [%e,%e]", theta, theta_min, theta_max);
    exit(1);
  }
  const int q = nl1*nlsize*nccl_size + nl2*nccl_size + ni; // cross redshift bin not supported so
                                                           // not using N_CCL(ni, nj) instead of ni
  if (q > NSIZE - 1)
  {
    log_fatal("internal logic error in selecting bin number");
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
  const int ncg_size = tomo.cg_clustering_Npowerspectra;
  const int NSIZE = nlsize*ncg_size;
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
        // ----------------------------------------------------------------------------------
        // Power spectrum on logarithmic bins (begins)
        // ----------------------------------------------------------------------------------
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-variable"
        {
          const double Z1 = ZCGCL1(0);
          const double Z2 = ZCGCL2(0);
          double init_static_vars_only = C_cg_tomo_limber(limits.LMIN_tab + 1, 0, Z1, Z2);
        }
        #pragma GCC diagnostic pop
        #pragma omp parallel for collapse(3)
        for (int i=0; i<nlsize; i++)  
        {
          for (int j=0; j<ncg_size; j++)
          {
            for (int p=0; p<ntheta; p++)
            {
              const int q = i*ncg_size+ j;
              const double Z1 = ZCGCL1(j);
              const double Z2 = ZCGCL2(j);
              const double l = exp(lnrc + (p - nc)*dlnl);
              lP[q][p] = (l > limits.LMIN_tab) ? l*C_cg_tomo_limber(l, i, Z1, Z2) :
                l*C_cg_tomo_limber_nointerp(l, i, Z1, Z2, use_linear_ps_limber, 0);
            }
          }
        }
        // ----------------------------------------------------------------------------------
        // Power spectrum on logarithmic bins (ends)
        // ----------------------------------------------------------------------------------
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
        tab[j] = (double**) malloc(sizeof(double*)*1);
        for (int i=0; i<1; i++)
        {
          tab[j][i] = (double*) malloc(sizeof(double)*ntheta);
        }
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
          const double t = exp((nc-k)*dlnl-lnrc); 
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
        free(tab[j]);
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
    update_galpara(&G);
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  if (nl < 0 || nl > nlsize - 1)
  {
    log_fatal("error in selecting bin number nl = %d (max %d)", nl, nlsize);
    exit(1); 
  } 
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d (max %d)", ni, redshift.clustering_nbin);
    exit(1); 
  } 
  if (nj < 0 || nj > tomo.clustering_Nbin - 1)
  {
    log_fatal("error in selecting bin number nj = %d (max %d)", nj, tomo.clustering_Nbin);
    exit(1); 
  }
  const double lntheta = log(theta);
  if (lntheta < lnthetamin || lntheta > lnthetamax)
  {
    const double theta = exp(lntheta);
    const double theta_min = exp(lnthetamin);
    const double theta_max = exp(lnthetamax);
    log_fatal("theta = %e outside look-up table range [%e,%e]", theta, theta_min, theta_max);
    exit(1);
  }
  const int ntomo = N_CGCL(ni, nj);
  if (!(ntomo>0))
  {
    return 0.0;
  }
  else
  {
    const int q = nl*ncg_size + ntomo;
    if (q > NSIZE-1)
    {
      log_fatal("internal logic error in selecting bin number");
      exit(1);
    }
    return interpol(table[q], ntheta, lnthetamin, lnthetamax, dlntheta, lntheta, 0, 0);
  }
}

// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
// Limber Approximation (Angular Power Spectrum)
// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// cluster lensing 
// ----------------------------------------------------------------------------
// nl = lambda_obs bin, ni = cluster redshift bin, nj = source redshift bin

double int_for_C_cs_tomo_limber(double a, void* params)
{ 
  if (!(a>0) || !(a<1))
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
  const int include_1h_term = 1;

  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;  
  const double PCM = binned_p_cm(k, a, nl, ni, include_1h_term, use_linear_ps);

  const double WCL = W_cluster(ni, a, hoverh0);
  const double WK = W_kappa(a, fK, nj);
  
  return WCL*WK*PCM*chidchi.dchida/(fK*fK);
}

double C_cs_tomo_limber_nointerp(const double l, const int nl, const int ni, const int nj, 
const int use_linear_ps, const int init_static_vars_only)
{
  if (nl < 0 || nl > Cluster.N200_Nbin - 1)
  {
    log_fatal("error in selecting bin number nl = %d (max %d)", nl, Cluster.N200_Nbin);
    exit(1); 
  } 
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d (max %d)", ni, redshift.clustering_nbin);
    exit(1); 
  } 
  if (nj < 0 || nj > redshift.shear_nbin - 1)
  {
    log_fatal("error in selecting bin number nj = %d (max %d)", nj, redshift.shear_nbin);
    exit(1); 
  }
  double ar[5] = {(double) nl, (double) ni, (double) nj, l, (double) use_linear_ps};
  const double zmin = tomo.cluster_zmin[ni];
  const double zmax = tomo.cluster_zmax[ni];
  if ((zmin > zmax) || !(zmin>0) || !(zmax>0))
  {
    log_fatal("error in selecting redshift range (zmin,zmax) = [%e,%e]", zmin, zmax);
    exit(1); 
  } 
  const double amin = 1./(1. + zmax);
  const double amax = 1./(1. + zmin);
  return (init_static_vars_only == 1) ? int_for_C_cs_tomo_limber(amin, (void*) ar) :
    int_gsl_integrate_low_precision(int_for_C_cs_tomo_limber, (void*) ar, amin, amax, NULL, 
      GSL_WORKSPACE_SIZE);
}

double C_cs_tomo_limber(const double l, const int nl, const int ni, const int nj)
{
  static cosmopara C;
  static galpara G;
  static nuisancepara N;
  static double** table = 0;

  const int nlsize = Cluster.N200_Nbin;
  const int ngammat_size = tomo.cgl_Npowerspectra;
  const int NSIZE = nlsize*ngammat_size;
  const int nell = Ntable.N_ell;
  const double lnlmin = log(fmax(limits.LMIN_tab, 1.0));
  const double lnlmax = log(limits.LMAX);
  const double dlnl = (lnlmax - lnlmin)/((double) nell - 1.0);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for (int i=0; i<NSIZE; i++)
    {
      table[i] = (double*) malloc(sizeof(double)*nell);
    }
  }
  if (recompute_cs(C, G, N))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      const int ZC = ZCL(0);
      const int ZS = ZCS(0);
      const double ll = exp(lnlmin);
      double init_static_vars_only = 
        C_cs_tomo_limber_nointerp(ll, 0, ZC, ZS, use_linear_ps_limber, 1);
    }
    #pragma GCC diagnostic pop
    #pragma omp parallel for collapse(3)
    for (int i=0; i<nlsize; i++)
    { 
      for (int j=0; j<ngammat_size; j++)
      {
        for (int p=0; p<nell; p++)
        {
          const int ZC = ZCL(j);
          const int ZS = ZCS(j);
          const int q = i*ngammat_size + j;
          const double lnl = lnlmin + p*dlnl;
          const double ll = exp(lnl);
          table[q][p] = log(C_cs_tomo_limber_nointerp(ll, i, ZC, ZS, use_linear_ps_limber, 0));
        }
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  if (nl < 0 || nl > nlsize - 1)
  {
    log_fatal("error in selecting bin number nl = %d (max %d)", nl, nlsize);
    exit(1); 
  } 
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d (max %d)", ni, redshift.clustering_nbin);
    exit(1); 
  } 
  if (nj < 0 || nj > redshift.shear_nbin - 1)
  {
    log_fatal("error in selecting bin number nj = %d (max %d)", nj, redshift.shear_nbin);
    exit(1); 
  }
  const double lnl = log(l);
  if (lnl < lnlmin)
  {
    log_warn("l = %e < l_min = %e. Extrapolation adopted", l, exp(lnlmin));
  }
  const int ntomo = N_cgl(ni, nj);
  if (!(ntomo>0)) 
  {
    return 0.0;
  }
  else
  {
    const int q = nl*ngammat_size + ntomo;
    if (q > NSIZE - 1)
    {
      log_fatal("internal logic error in selecting bin number");
      exit(1);
    }
    const double f1 = exp(interpol(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 1, 1));
    return isnan(f1) ? 0.0 : f1;
  }
}

// ---------------------------------------------------------------------------------------------
// Cluster clustering 
// ---------------------------------------------------------------------------------------------

double int_for_C_cc_tomo_limber(double a, void* params)
{ // nl{1,2} = lambda_obs bins, n{i,j} = cluster redshift bins
  if (!(a>0) || !(a<1))
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

  const double res = W_cluster(ni, a, hoverh0)*W_cluster(nj, a, hoverh0);
  const double PCC = binned_p_cc(k, a, nl1, nl2, use_linear_ps);
  return (res == 0.0) ? 0.0 : res*PCC*chidchi.dchida/(fK*fK);
}

double C_cc_tomo_limber_nointerp(const double l, const int nl1, const int nl2, const int ni, 
const int nj, const int use_linear_ps, const int init_static_vars_only)
{ // nl{1,2} = lambda_obs bins, n{i,j} = cluster redshift bins
  if (nl1 < 0 || nl1 > Cluster.N200_Nbin - 1)
  {
    log_fatal("error in selecting bin number nl1 = %d (max %d)", nl1, Cluster.N200_Nbin);
    exit(1); 
  } 
  if (nl2 < 0 || nl2 > Cluster.N200_Nbin - 1)
  {
    log_fatal("error in selecting bin number nl2 = %d (max %d)", nl2, Cluster.N200_Nbin);
    exit(1); 
  } 
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d (max %d)", ni, redshift.clustering_nbin);
    exit(1); 
  } 
  if (nj < 0 || nj > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number nj = %d (max %d)", nj, redshift.clustering_nbin);
    exit(1); 
  } 
  double ar[6] = {(double) nl1, (double) nl2, (double) ni, (double) nj, l, (double) use_linear_ps};  
  const double zmin = fmax(tomo.cluster_zmin[ni], tomo.cluster_zmin[nj]);
  const double zmax = fmin(tomo.cluster_zmax[ni], tomo.cluster_zmax[nj]);
  if ((zmin > zmax) || !(zmin>0) || !(zmax>0))
  {
    log_fatal("error in selecting redshift range (zmin,zmax) = [%e,%e]", zmin, zmax);
    exit(1); 
  } 
  const double amin = 1.0/(1.0 + zmax);
  const double amax = 1.0/(1.0 + zmin);
  return (init_static_vars_only == 1) ? int_for_C_cc_tomo_limber(amin, (void*) ar) :
    int_gsl_integrate_low_precision(int_for_C_cc_tomo_limber, (void*) ar, amin, amax, NULL,
      GSL_WORKSPACE_SIZE);
}

double C_cc_tomo_limber(const double l, const int nl1, const int nl2, const int ni, 
const int nj) 
{ // nl{1,2} = lambda_obs bins, n{i,j} = cluster redshift bins
  static cosmopara C;
  static nuisancepara N;
  static double** table = 0;

  const int nlsize = Cluster.N200_Nbin;
  const int nccl_size = tomo.cc_clustering_Npowerspectra; 
  const int NSIZE = nlsize*nlsize*nccl_size;
  const int nell = Ntable.N_ell;
  const double lnlmin = log(fmax(limits.LMIN_tab, 1.0));
  const double lnlmax = log(limits.LMAX);
  const double dlnl = (lnlmax - lnlmin)/((double) nell - 1.0); 

  if (table == 0)
  { 
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for (int i=0; i<NSIZE; i++)
    {
      table[i] = (double*) malloc(sizeof(double)*nell);
    }   
  }
  if (recompute_cc(C, N))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
    {
      double init_static_vars_only =
        C_cc_tomo_limber_nointerp(exp(lnlmin), 0, 0, 0, 0, use_linear_ps_limber, 1);
    }
    #pragma GCC diagnostic pop
    #pragma GCC diagnostic pop
    #pragma omp parallel for collapse(3)
    for (int i=0; i<nlsize; i++)
    {
      for (int k=0; k<nccl_size; k++)
      {
        for (int p=0; p<nell; ++p)
        {
          for (int j=i; j<nlsize; j++)
          {
            const int ZCCL1  = k; // cross redshift bin not supported so not using ZCCL1(k)
            const int ZCCL2  = k; // cross redshift bin not supported so not using ZCCL2(k)
            const int q      = i*nlsize*nccl_size + j*nccl_size + k;
            const int qstar  = j*nlsize*nccl_size + i*nccl_size + k;
            const double lnl = lnlmin + p*dlnl;
            table[q][p] = 
              log(C_cc_tomo_limber_nointerp(exp(lnl), i, j, ZCCL1, ZCCL2, use_linear_ps_limber, 0));
            table[qstar][p] = table[q][p];
          }
        }
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  if (ni != nj)
  {
    log_fatal("ni != nj tomography not supported");
    exit(1);
  } 
  if (nl1 < 0 || nl1 > nlsize - 1)
  {
    log_fatal("error in selecting bin number nl1 = %d (max %d)", nl1, nlsize);
    exit(1); 
  } 
  if (nl2 < 0 || nl2 > nlsize - 1)
  {
    log_fatal("error in selecting bin number nl2 = %d (max %d)", nl2, nlsize);
    exit(1); 
  } 
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d (max %d)", ni, redshift.clustering_nbin);
    exit(1); 
  } 
  if (nj < 0 || nj > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number nj = %d (max %d)", nj, redshift.clustering_nbin);
    exit(1); 
  } 
  const double lnl = log(l);
  if (lnl < lnlmin)
  {
    log_warn("l = %e < l_min = %e. Extrapolation adopted", l, exp(lnlmin));
  }
  const int q = nl1*nlsize*nccl_size + nl2*nccl_size + ni; // cross redshift bin not supported so 
                                                           // not using N_CCL(ni, nj)
  if (q > NSIZE - 1)
  {
    log_fatal("internal logic error in selecting bin number");
    exit(1);
  }
  const double f1 = exp(interpol(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 1, 1));
  return isnan(f1) ? 0.0 : f1;
}


// ---------------------------------------------------------------------------------------------
// cluster x galaxy clustering
// ---------------------------------------------------------------------------------------------

double int_for_C_cg_tomo_limber(double a, void* params)
{ // nl = lambda_obs bin, ni = cluster redshift bin, nj = galaxy redshift bin
  if (!(a>0) || !(a<1))
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  double* ar = (double*) params;

  const int nl = (int) ar[0];
  const int ni = (int) ar[1];
  const int nj = (int) ar[2];
  const double ell = ar[3] + 0.5;
  const int use_linear_ps = (int) ar[4];

  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const  double fK  = f_K(chidchi.chi);
  const double k = ell/fK;
  const double PCG = binned_p_cg(k, a, nl, nj, use_linear_ps);
  const double tmp = W_cluster(ni, a, hoverh0)*W_HOD(a, nj, hoverh0);
  return tmp*PCG*chidchi.dchida/(fK*fK);
}

double C_cg_tomo_limber_nointerp(const double l, const int nl, const int ni, const int nj, 
const int use_linear_ps, const int init_static_vars_only)
{ // nl = lambda_obs bin, ni = cluster redshift bin, nj = galaxy redshift bin
  if (nl > Cluster.N200_Nbin - 1)
  {
    log_fatal("error in selecting bin number nl = %d (max %d)", nl, Cluster.N200_Nbin);
    exit(1); 
  } 
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d (max %d)", ni, redshift.clustering_nbin);
    exit(1); 
  } 
  if (nj < 0 || nj > tomo.clustering_Nbin - 1)
  {
    log_fatal("error in selecting bin number nj = %d (max %d)", nj, tomo.clustering_Nbin);
    exit(1); 
  }
  double ar[5] = {(double) nl, (double) ni, (double) nj, l, (double) use_linear_ps};
  const double zmin = fmax(tomo.cluster_zmin[ni], tomo.clustering_zmin[nj]);
  const double zmax = fmin(tomo.cluster_zmax[ni], tomo.clustering_zmax[nj]);
  if ((zmin > zmax) || !(zmin>0) || !(zmax>0))
  {
    log_fatal("error in selecting redshift range (zmin,zmax) = [%e,%e]", zmin, zmax);
    exit(1); 
  } 
  const double amin = 1.0/(1.0 + zmax);
  const double amax = 1.0/(1.0 + zmin);
  if (init_static_vars_only == 1)
  {
    return int_for_C_cg_tomo_limber(amin, (void*) ar);
  }
  else
  {
    return int_gsl_integrate_low_precision(int_for_C_cg_tomo_limber, (void*) ar, amin, amax,
      NULL, GSL_WORKSPACE_SIZE);
  }
}

double C_cg_tomo_limber(const double l, const int nl, const int ni, const int nj)
{ // nl = lambda_obs bin, ni = cluster redshift bin, nj = galaxy redshift bin
  static cosmopara C;
  static galpara G;
  static nuisancepara N;
  static double** table = 0;

  const int nlsize = Cluster.N200_Nbin;
  const int ncg_size = tomo.cg_clustering_Npowerspectra;
  const int NSIZE = nlsize*ncg_size;
  const int nell = Ntable.N_ell;
  const double lnlmin = log(fmax(limits.LMIN_tab, 1.0));
  const double lnlmax = log(limits.LMAX);
  const double dlnl = (lnlmax - lnlmin)/((double) nell - 1.0);

  if (table == 0)
  { 
    table = (double**) malloc(sizeof(double*)*NSIZE);
    for (int i=0; i<NSIZE; i++)
    {
      table[i] = (double*) malloc(sizeof(double)*nell);
    }      
  }
  if (recompute_cg(C, G, N))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
    {
      const double Z1 = ZCGCL1(0);
      const double Z2 = ZCGCL2(0);
      const double ll = exp(lnlmin);
      double init_static_vars_only = 
        C_cg_tomo_limber_nointerp(ll, 0, Z1, Z2, use_linear_ps_limber, 1);
    }
    #pragma GCC diagnostic pop
    #pragma omp parallel for collapse(3)
    for (int i=0; i<nlsize; i++) 
    {
      for (int j=0; j<ncg_size; j++)
      {
        for (int p=0; p<nell; ++p)
        {
          const double Z1 = ZCGCL1(j);
          const double Z2 = ZCGCL2(j);
          const int q = i*ncg_size + j;
          const double lnl = lnlmin + p*dlnl;
          const double ll = exp(lnl);
          table[q][p] = log(C_cg_tomo_limber_nointerp(ll, i, Z1, Z2, use_linear_ps_limber, 0));
        }
      }
    }
    update_galpara(&G);
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  if (nl < 0 || nl > nlsize - 1)
  {
    log_fatal("error in selecting bin number nl = %d (max %d)", nl, nlsize);
    exit(1); 
  } 
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d (max %d)", ni, redshift.clustering_nbin);
    exit(1); 
  } 
  if (nj < 0 || nj > tomo.clustering_Nbin - 1)
  {
    log_fatal("error in selecting bin number nj = %d (max %d)", nj, tomo.clustering_Nbin);
    exit(1); 
  }
  const double lnl = log(l);
  if (lnl < lnlmin)
  {
    log_warn("l = %e < l_min = %e. Extrapolation adopted", l, exp(lnlmin));
  }
  const int ntomo = N_CGCL(ni, nj);
  if (!(ntomo>0))
  {
    return 0.0;
  }
  else
  {
    const int q = nl*ncg_size + ntomo;
    if (q > NSIZE-1)
    {
      log_fatal("internal logic error in selecting bin number");
      exit(1);
    }
    const double f1 = exp(interpol(table[q], nell, lnlmin, lnlmax, dlnl, lnl, 1, 1));
    return isnan(f1) ? 0.0 : f1;
  }
}

// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------
// Non Limber (Angular Power Spectrum)
// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------

void f_chi_for_Psi_cluster_cl(double *const chi, const int Nchi, double *const fchi, const int ni,
const int nl, const double zmin, const double zmax)
{
  const double real_coverH0 = cosmology.coverH0 / cosmology.h0; // unit Mpc
  {
    const int i = 0;
    const double a = a_chi(chi[i]/real_coverH0 /* convert unit to c/H0 */);
    const double z = 1./a - 1.;
    const double hoh0 = hoverh0(a);
    const double tmp1 = zdistr_cluster(z, ni);
    const double pf = (tmp1 < 0.) ? 0 : tmp1; ; // get rid of unphysical negatives
    fchi[i] = chi[i]*pf*growfac(a)*weighted_B1(nl, z)*hoh0/real_coverH0;
    if ((z < zmin) || (z > zmax))
    {
      fchi[i] = 0.0;
    }
  }
  #pragma omp parallel for
  for (int i=1; i<Nchi; i++) 
  {
    const double a = a_chi(chi[i]/real_coverH0 /* convert unit to c/H0 */);
    const double z = 1./a - 1.;
    const double hoh0 = hoverh0(a);
    const double tmp1 = zdistr_cluster(z, ni);
    const double pf = (tmp1 < 0.) ? 0 : tmp1; // get rid of unphysical negatives
    fchi[i] = chi[i]*pf*growfac(a)*weighted_B1(nl, z)*hoh0/real_coverH0;
    if ((z < zmin) || (z > zmax))
    {
      fchi[i] = 0.0;
    }
  }
}

void f_chi_for_Psi_cluster_cl_RSD(double *const chi, const int Nchi, double *const fchi,
const int ni, const int nl __attribute__((unused)), const double zmin, const double zmax)
{
  const double real_coverH0 = cosmology.coverH0 / cosmology.h0;
  {
    const int i = 0;
    const double a = a_chi(chi[i]/real_coverH0 /* convert unit to c/H0 */);
    const double z = 1./a - 1.;
    const double hoh0 = hoverh0(a);
    const double tmp1 = zdistr_cluster(z, ni);
    const double pf = (tmp1 < 0.) ? 0 : tmp1; // get rid of unphysical negatives
    struct growths tmp2 = growfac_all(a);
    fchi[i] = -chi[i]*pf*tmp2.D*tmp2.f*hoh0/real_coverH0;
    if ((z < zmin) || (z > zmax))
    {
      fchi[i] = 0.;
    }
  }
  #pragma omp parallel for
  for (int i=1; i<Nchi; i++)
  {
    const double a = a_chi(chi[i]/real_coverH0 /* convert unit to c/H0 */);
    const double z = 1./a - 1.;
    const double hoh0 = hoverh0(a);
    const double tmp1 = zdistr_cluster(z, ni);
    const double pf = (tmp1 < 0.) ? 0 : tmp1; // get rid of unphysical negatives
    struct growths tmp2 = growfac_all(a);
    fchi[i] = -chi[i]*pf*tmp2.D*tmp2.f*hoh0/real_coverH0;
    if ((z < zmin) || (z > zmax))
    {
      fchi[i] = 0.;
    }
  }
}

void f_chi_for_Psi_cluster_cl_Mag(double *const chi, const int Nchi, double *const fchi,
const int ni, const int nl, const double zmax)
{
  const double real_coverH0 = cosmology.coverH0 / cosmology.h0;
  {
    const int i = 0;
    const double a = a_chi(chi[i]/real_coverH0 /* convert unit to c/H0 */);
    const double z = 1.0/a - 1.0;
    const double fK = f_K(chi[i]/real_coverH0);
    const double wmag = W_mag_cluster(a, fK, ni, nl);
    const double window_M = wmag/fK/(real_coverH0*real_coverH0);
    fchi[i] = window_M*growfac(a); // unit [Mpc^-2]    
    if (z > zmax)
    {
      fchi[i] = 0.;
    }
  }
  #pragma omp parallel for
  for (int i=1; i<Nchi; i++)
  {
    const double a = a_chi(chi[i]/real_coverH0 /* convert unit to c/H0 */);
    const double z = 1.0/a - 1.0;
    const double fK = f_K(chi[i]/real_coverH0);
    const double wmag = W_mag_cluster(a, fK, ni, nl);
    const double window_M = wmag/fK/(real_coverH0*real_coverH0);
    fchi[i] = window_M*growfac(a); // unit [Mpc^-2]    
    if (z > zmax)
    {
      fchi[i] = 0.;
    }
  }
}

void C_cc_tomo(int L, const int nl1, const int nl2, const int ni, const int nj, double *const Cl,
double dev, const double tol)
{ // nl{1,2} = lambda_obs bins, n{i,j} = cluster redshift bins
  if (ni != nj)
  {
    log_fatal("cross-spectrum not supported");
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
  const double real_coverH0 = cosmology.coverH0 / cosmology.h0;
  const double chi_min = chi(1.0/(1.0 + 0.002))*real_coverH0; // DIMENSIONELESS
  const double chi_max = chi(1.0/(1.0 + 4.0))*real_coverH0;    // DIMENSIONELESS
  const double dlnchi = log(chi_max / chi_min) / ((double) Nchi - 1.0);
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
  { // chi_min and chi_max may be cosmology dependent
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

  const double zmin = tomo.cluster_zmin[ni];
  const double zmax = tomo.cluster_zmax[ni];
  f_chi_for_Psi_cluster_cl(chi_ar, Nchi, f1_chi, ni, nl1, zmin, zmax);
  f_chi_for_Psi_cluster_cl_RSD(chi_ar, Nchi, f1_chi_RSD, ni, nl1, zmin, zmax);
  if (INCLUDE_MAG_IN_C_CC_NONLIMBER == 1)
  {
    f_chi_for_Psi_cluster_cl_Mag(chi_ar, Nchi, f1_chi_Mag, ni, nl1, zmax);
  }
  if (nl1 != nl2) 
  {
    f_chi_for_Psi_cluster_cl(chi_ar, Nchi, f2_chi, ni, nl2, zmin, zmax);
    f_chi_for_Psi_cluster_cl_RSD(chi_ar, Nchi, f1_chi_RSD, ni, nl2, zmin, zmax);
    if (INCLUDE_MAG_IN_C_CC_NONLIMBER == 1)
    {
      f_chi_for_Psi_cluster_cl_Mag(chi_ar, Nchi, f1_chi_Mag, ni, nl2, zmax);
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
    
  while ((fabs(dev) > tol) & (L < limits.LMAX_NOLIMBER))
  { 
    for (int i=0; i<Nell_block; i++)
    {
      ell_ar[i]=i+i_block*Nell_block;
    }
    i_block++;
    if (L >= limits.LMAX_NOLIMBER - Nell_block)
    { //Xiao: break before memory leak in next iteration
      break;
    }
    L = i_block*Nell_block - 1;

    cfftlog_ells(chi_ar, f1_chi, Nchi, &cfg, ell_ar, Nell_block, k1, Fk1);
    cfftlog_ells_increment(chi_ar, f1_chi_RSD, Nchi, &cfg_RSD, ell_ar, Nell_block, k1, Fk1);
    if (nl1 != nl2)
    {
      cfftlog_ells(chi_ar, f2_chi, Nchi, &cfg, ell_ar, Nell_block, k2, Fk2);
      cfftlog_ells_increment(chi_ar, f2_chi_RSD, Nchi, &cfg_RSD, ell_ar, Nell_block, k2, Fk2);
    }
    if (INCLUDE_MAG_IN_C_CC_NONLIMBER == 1)
    {
      cfftlog_ells(chi_ar, f1_chi_Mag, Nchi, &cfg_Mag, ell_ar, Nell_block, k1, Fk1_Mag);
      if (nl1 != nl2) 
      {
        cfftlog_ells(chi_ar, f2_chi_Mag, Nchi, &cfg_Mag, ell_ar, Nell_block, k2, Fk2_Mag);
      }
    }
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
    {
      double init_static_vars_only = p_lin(k1[0][0]*real_coverH0, 1.0);
      init_static_vars_only = C_cc_tomo_limber_nointerp((double) ell_ar[0], nl1, nl2, ni, ni, 0, 1);
    }
    #pragma GCC diagnostic pop
    #pragma GCC diagnostic pop
    #pragma omp parallel for
    for (int i=0; i<Nell_block; i++)
    {
      double cl_temp = 0.0;
      for (int j=0; j<Nchi; j++)
      {
        if (INCLUDE_MAG_IN_C_CC_NONLIMBER == 1)
        {
          const double ell_prefactor = ell_ar[i]*(ell_ar[i]+1.);
          Fk1[i][j] += ell_prefactor*Fk1_Mag[i][j]/(k1[i][j]*k1[i][j]); 
          if (nl1 != nl2) 
          {
            Fk2[i][j] += ell_prefactor*Fk2_Mag[i][j]/(k2[i][j]*k2[i][j]) ;
          }
        }      
        // ------------------------------------------------------------------------------------
        const double k1_cH0 = k1[i][j] * real_coverH0;
        const double PK = p_lin(k1_cH0, 1.0);
        const double k1_cH03 = k1_cH0*k1_cH0*k1_cH0;
        cl_temp += (nl1 == nl2) ? Fk1[i][j]*Fk1[i][j]*k1_cH03*PK : Fk1[i][j]*Fk2[i][j]*k1_cH03*PK;
      }
      Cl[ell_ar[i]] = cl_temp * dlnk * 2./M_PI + 
        C_cc_tomo_limber_nointerp((double) ell_ar[i], nl1, nl2, ni, ni, 0, 0) 
       -C_cc_tomo_limber_nointerp((double) ell_ar[i], nl1, nl2, ni, ni, 1, 0);
    }
    dev = Cl[L]/C_cc_tomo_limber_nointerp((double) L, nl1, nl2, ni, ni, 0, 0) - 1;
  }   
  L++;
  
  Cl[limits.LMAX_NOLIMBER] = C_cc_tomo_limber((double) limits.LMAX_NOLIMBER, nl1, nl2, ni, ni);
  #pragma omp parallel for
  for (int l=L; l<limits.LMAX_NOLIMBER; l++)
  {
    Cl[l] = (l > limits.LMIN_tab) ? C_cc_tomo_limber((double) l, nl1, nl2, ni, ni) :
      C_cc_tomo_limber_nointerp((double) l, nl1, nl2, ni, ni, use_linear_ps_limber, 0);
  }
}

// ------------------------------------------------------------------------------------

void C_cg_tomo(int L, const int nl, const int ni, const int nj, double *const Cl, double dev,
const double tol)
{ // nl = lambda_obs bin, ni = cluster redshift bin, nj = galaxy redshift bin
  if (ni != nj)
  {
    log_fatal("cross-spectrum not supported");
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
  const double real_coverH0 = cosmology.coverH0 / cosmology.h0;
  const double chi_min = chi(1./(1.0 + 0.002))*real_coverH0; // DIMENSIONELESS
  const double chi_max = chi(1./(1.0 + 4.0))*real_coverH0;    // DIMENSIONELESS
  const double dlnchi = log(chi_max/chi_min) / ((double) Nchi - 1.0);
  const double dlnk = dlnchi;

  if (k1 == 0)
  { // COCOA: no need to create/destroy arrays with same size at every call
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

  const double zmin = fmax(tomo.cluster_zmin[ni], tomo.clustering_zmin[nj]);
  const double zmax = fmin(tomo.cluster_zmax[ni], tomo.clustering_zmax[nj]);

  f_chi_for_Psi_cluster_cl(chi_ar, Nchi, f1_chi, ni, nl, zmin, zmax);
  f_chi_for_Psi_cluster_cl_RSD(chi_ar, Nchi, f1_chi_RSD, ni, nl, zmin, zmax);
  if (INCLUDE_MAG_IN_C_CG_NONLIMBER)
  {
    f_chi_for_Psi_cluster_cl_Mag(chi_ar, Nchi, f1_chi_Mag, ni, nl, zmax);
  }
  
  f_chi_for_Psi_cl(chi_ar, Nchi, f2_chi, nj, zmin, zmax);
  f_chi_for_Psi_cl_RSD(chi_ar, Nchi, f2_chi_RSD, nj, zmin, zmax);
  if (INCLUDE_MAG_IN_C_CG_NONLIMBER)
  {
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
    
  while ((fabs(dev) > tol) & (L < limits.LMAX_NOLIMBER))
  { 
    for (int i=0; i<Nell_block; i++)
    {
      ell_ar[i]=i+i_block*Nell_block;
    } 
    i_block++;  
    if (L >= limits.LMAX_NOLIMBER - Nell_block)
    { //Xiao: break before memory leak in next iteration
      break;
    }
    L = i_block*Nell_block - 1;

    cfftlog_ells(chi_ar, f1_chi, Nchi, &cfg, ell_ar, Nell_block, k1, Fk1);
    cfftlog_ells_increment(chi_ar, f1_chi_RSD, Nchi, &cfg_RSD, ell_ar, Nell_block, k1, Fk1);
    cfftlog_ells(chi_ar, f2_chi, Nchi, &cfg, ell_ar, Nell_block, k2, Fk2);
    cfftlog_ells_increment(chi_ar, f2_chi_RSD, Nchi, &cfg_RSD, ell_ar, Nell_block, k2, Fk2);
    if (INCLUDE_MAG_IN_C_CG_NONLIMBER)
    {
      cfftlog_ells(chi_ar, f1_chi_Mag, Nchi, &cfg_Mag, ell_ar, Nell_block, k1, Fk1_Mag);
      cfftlog_ells(chi_ar, f2_chi_Mag, Nchi, &cfg_Mag, ell_ar, Nell_block, k2, Fk2_Mag);
    }

    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
    {
      double init_static_vars_only = p_lin(k1[0][0]*real_coverH0, 1.0);
      init_static_vars_only = C_cg_tomo_limber_nointerp((double) ell_ar[0], nl, ni, nj, 0, 1);
    }
    #pragma GCC diagnostic pop
    #pragma GCC diagnostic pop
    #pragma omp parallel for
    for (int i=0; i<Nell_block; i++)
    {
      double cl_temp = 0.;
      for (int j=0; j<Nchi; j++)
      {
        if (INCLUDE_MAG_IN_C_CG_NONLIMBER)
        {
          const double ell_prefactor = ell_ar[i]*(ell_ar[i]+1.);
          Fk1[i][j] += (ell_prefactor/(k1[i][j]*k1[i][j])*Fk1_Mag[i][j]);
          Fk2[i][j] += gbias.b_mag[nj]*(ell_prefactor/(k2[i][j]*k2[i][j])*Fk2_Mag[i][j]);
        }
        // ------------------------------------------------------------------------------------
        const double k1_cH0 = k1[i][j] * real_coverH0;
        cl_temp += Fk1[i][j]*Fk2[i][j]*k1_cH0*k1_cH0*k1_cH0 *p_lin(k1_cH0, 1.0);
      }
      Cl[ell_ar[i]] = cl_temp * dlnk * 2./M_PI + 
        C_cg_tomo_limber_nointerp((double) ell_ar[i], nl, ni, nj, 0, 0) 
       -C_cg_tomo_limber_nointerp((double) ell_ar[i], nl, ni, nj, 1, 0);
    }
    dev = Cl[L]/C_cg_tomo_limber_nointerp((double) L, nl, ni, nj, 0, 0) - 1;
  }   
  L++;

  Cl[limits.LMAX_NOLIMBER] = C_cg_tomo_limber((double) limits.LMAX_NOLIMBER, nl, ni, nj);
  #pragma omp parallel for
  for (int l=L; l<limits.LMAX_NOLIMBER; l++)
  {
    Cl[l] = (l > limits.LMIN_tab) ? C_cg_tomo_limber((double) l, nl, ni, nj) :
      C_cg_tomo_limber_nointerp((double) l, nl, ni, nj, use_linear_ps_limber, 0);
  }
}

// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
// cluster number counts
// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
// nl = lambda_obs bin, ni = cluster redshift bin

double int_for_binned_N(double a, void* params)
{
  if (!(a>0))
  {
    log_fatal("a > 0 not true");
    exit(1);
  }
  double* ar = (double*) params;   
  const int nl = (int) ar[0];
  const int nz = (int) ar[1];
  const int interpolate_survey_area = (int) ar[2];
  const double z = 1.0/a - 1.0 ;
  const double dzda = 1.0/(a*a); 
  const double norm = get_area(z, interpolate_survey_area);  

  double tmp_param[1] = {(double) nz};
  return dV_cluster(z, (void*) tmp_param)*dzda*binned_Ndensity(nl, z)*norm;
}

double binned_N_nointerp(const int nl, const int nz, const int interpolate_survey_area, 
const int init_static_vars_only)
{
  double params[3] = {(double) nl, (double) nz, interpolate_survey_area};
  const double tmp = 4.0*M_PI/41253.0;
  const double amin = 1.0/(1.0 + tomo.cluster_zmax[nz]);
  const double amax = 1.0/(1.0 + tomo.cluster_zmin[nz]);
  return (init_static_vars_only == 1) ? int_for_binned_N(amin, (void*) params) :
    tmp*int_gsl_integrate_low_precision(int_for_binned_N, (void*) params, amin, amax, NULL, 
      GSL_WORKSPACE_SIZE);
}

double binned_N(const int nl, const int nz)
{
  static cosmopara C;
  static nuisancepara N;
  static double** table;

  const int N_l = Cluster.N200_Nbin;
  const int N_z = redshift.clustering_nbin;

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*N_l);
    for (int i=0; i<N_l; i++)
    {
      table[i] = (double*) malloc(sizeof(double)*N_z);
    }
  }
  if (recompute_clusters(C, N))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
    {
      double init_static_vars_only = binned_N_nointerp(0, 0, Cluster.interpolate_survey_area, 1);
    }
    #pragma GCC diagnostic pop
    #pragma GCC diagnostic pop
    #pragma omp parallel for collapse(2)
    for (int i=0; i<N_l; i++)
    {
      for (int j=0; j<N_z; j++)
      {
        table[i][j] = binned_N_nointerp(i, j, Cluster.interpolate_survey_area, 0);
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  if (nl < 0 || nl > N_l - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }
  if (nz < 0 || nz > N_z - 1)
  {
    log_fatal("error in selecting bin number");
    exit(1);
  }
  return table[nl][nz];
}