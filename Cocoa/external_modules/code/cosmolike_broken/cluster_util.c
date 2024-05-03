#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>

#include "basics.h"
#include "bias.h"
#include "cosmo3D.h"
#include "cluster_util.h"
#include "halo.h"
#include "pt_cfastpt.h"
#include "recompute.h"
#include "radial_weights.h"
#include "redshift_spline.h"
#include "structs.h"
#include "tinker_emulator.h"

#include "log.c/src/log.h"

static int GSL_WORKSPACE_SIZE = 250;
static double M_PIVOT = 5E14; // Msun/h

static int has_b2_galaxies()
{
  int res = 0;
  for(int i=0; i<tomo.clustering_Nbin; i++) 
  {
    if (gbias.b2[i])
    {
      res = 1;
    }
  }
  return res;
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// BUZZARD binned_P_lambda_obs_given_M
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

double buzzard_P_lambda_obs_given_M(const double obs_lambda, const double M, const double z)
{        
  const double lnlm0 = nuisance.cluster_MOR[0]; 
  const double Alm = nuisance.cluster_MOR[1];
  const double sigma_lnlm_intrinsic = nuisance.cluster_MOR[2];
  const double Blm = nuisance.cluster_MOR[3];
  
  const double Mpiv = M_PIVOT; //Msun/h
  
  const double lnlm = lnlm0 + (Alm)*log(M/Mpiv) + Blm*log((1 + z)/1.45);
  
  const double sigma_total = (lnlm>0) ? 
    sqrt(sigma_lnlm_intrinsic*sigma_lnlm_intrinsic+(exp(lnlm)-1.)/exp(2*lnlm)) : 
    sqrt(sigma_lnlm_intrinsic*sigma_lnlm_intrinsic);
  
  const double x = 1.0/2.0*(log(obs_lambda)-lnlm)*(log(obs_lambda)-lnlm)/pow(sigma_total,2.0);
  
  return exp(-x)/M_SQRTPI/M_SQRT2/sigma_total/obs_lambda;
}

static double buzzard_P_lambda_obs_given_M_wrapper(double obs_lambda, void* params)
{      
  double* ar = (double*) params;
  const double M = ar[0];
  const double z = ar[1];
  return buzzard_P_lambda_obs_given_M(obs_lambda, M, z);
}

// \int_(bin_lambda_obs_min)^(bin_lambda_obs_max) \dlambda_obs P(\lambda_obs|M)
// (see for example https://arxiv.org/pdf/1810.09456.pdf - eq3 qnd 6) 
double buzzard_binned_P_lambda_obs_given_M(const int nl, const double M, const double z, 
  const int init_static_vars_only)
{
  double params[2] = {M, z};
  const double bin_lambda_obs_min = Cluster.N_min[nl];
  const double bin_lambda_obs_max = Cluster.N_max[nl];
  return (init_static_vars_only == 1) ? 
    buzzard_P_lambda_obs_given_M_wrapper(bin_lambda_obs_min, (void*) params) :
    int_gsl_integrate_medium_precision(buzzard_P_lambda_obs_given_M_wrapper,
      (void*) params, bin_lambda_obs_min, bin_lambda_obs_max, NULL, GSL_WORKSPACE_SIZE);
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// SDSS binned_P_lambda_obs_given_M
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

// Cocoa: we try to avoid reading of files in the cosmolike_core code 
// Cocoa: (reading is done in the C++/python interface)
void setup_SDSS_P_true_lambda_given_mass(int* io_nintrinsic_sigma, double** io_intrinsic_sigma, 
int* io_natsrgm, double** io_atsrgm, double** io_alpha, double** io_sigma, int io)
{
  static int nintrinsic_sigma;
  static int natsrgm;  
  static double* intrinsic_sigma = NULL;
  static double* atsrgm = NULL;
  static double* alpha = NULL;
  static double* sigma = NULL;

  if (io == 1) // IO == 1 IMPLES THAT IO_XXX WILL COPIED TO LOCAL XXX
  {
    if (intrinsic_sigma != NULL)
    {
      free(intrinsic_sigma);
    }
    else if (atsrgm != NULL)
    {
      free(atsrgm);
    }
    else if (alpha != NULL)
    {
      free(alpha);
    }
    else if (sigma != NULL)
    {
      free(sigma);
    }

    nintrinsic_sigma = (*io_nintrinsic_sigma);
    natsrgm = (*io_natsrgm);

    if (!(nintrinsic_sigma > 5) || !(natsrgm > 5))
    {
      log_fatal("array to small for 2D interpolation");
      exit(1);
    }

    intrinsic_sigma = (double*) malloc(nintrinsic_sigma*sizeof(double));
    atsrgm = (double*) malloc(natsrgm*sizeof(double));
    alpha = (double*) malloc(natsrgm*nintrinsic_sigma*sizeof(double));
    sigma = (double*) malloc(natsrgm*nintrinsic_sigma*sizeof(double));
    if (intrinsic_sigma == NULL || atsrgm == NULL || alpha == NULL || sigma == NULL)
    {
      log_fatal("fail allocation");
      exit(1);
    }
    
    for (int i=0; i<nintrinsic_sigma; i++)
    {
      intrinsic_sigma[i] = (*io_intrinsic_sigma)[i];
      for (int j=0; j<natsrgm; j++)
      {
        if (i == 0) 
        {
          atsrgm[j] = (*io_atsrgm)[j];
        }
        alpha[i*nintrinsic_sigma + j] = (*io_alpha)[i*nintrinsic_sigma + j];
        sigma[i*nintrinsic_sigma + j] = (*io_sigma)[i*nintrinsic_sigma + j];
      }
    }
  }
  else  // IO != 1 IMPLES THAT LOCAL XXX WILL BE COPIED TO IO_XXX
  {
    if (intrinsic_sigma == NULL || atsrgm == NULL || alpha == NULL || sigma == NULL ||
        io_intrinsic_sigma == NULL || io_atsrgm == NULL || io_alpha == NULL || io_sigma == NULL)
    {
      log_fatal("array/pointers not allocated\n");
      exit(1);
    }
    if ((*io_intrinsic_sigma) != NULL)
    {
      free((*io_intrinsic_sigma));
      (*io_intrinsic_sigma) = NULL;
    }
    else if ((*io_atsrgm) != NULL)
    {
      free((*io_atsrgm));
      (*io_atsrgm) = NULL;
    }
    else if ((*io_alpha) != NULL)
    {
      free((*io_alpha));
      (*io_alpha) = NULL;
    }
    else if ((*io_sigma) != NULL)
    {
      free((*io_sigma));
      (*io_sigma) = NULL;
    }

    (*io_intrinsic_sigma) = (double*) malloc(nintrinsic_sigma*sizeof(double));
    (*io_atsrgm) = (double*) malloc(natsrgm*sizeof(double));
    (*io_alpha) = (double*) malloc(nintrinsic_sigma*natsrgm*sizeof(double));
    (*io_sigma) = (double*) malloc(nintrinsic_sigma*natsrgm*sizeof(double));
    if (*io_intrinsic_sigma == NULL || *io_atsrgm == NULL || *io_alpha == NULL || *io_sigma == NULL)
    {
      log_fatal("fail allocation");
      exit(1);
    }

    for (int i=0; i<nintrinsic_sigma; i++)
    {
      (*io_intrinsic_sigma)[i] = intrinsic_sigma[i];
      for (int j=0; j<natsrgm; j++)
      {
        if (i == 0) 
        {
          (*io_atsrgm)[j] = atsrgm[j];
        }
        (*io_alpha)[i*nintrinsic_sigma + j] = alpha[i*nintrinsic_sigma + j];
        (*io_sigma)[i*nintrinsic_sigma + j] = sigma[i*nintrinsic_sigma + j];
      }
    }
  }
}

double SDSS_P_true_lambda_given_mass(const double true_lambda, const double mass, const double z)
{ // SKEW-NORMAL APPROXIMATION eq B1 of https://arxiv.org/pdf/1810.09456.pdf
  static int first = 0;
  static gsl_spline2d* falpha = NULL; // skewness of the skew-normal distribution
  static gsl_spline2d* fsigma = NULL; // variance of the skew-normal distribution
  if (first == 0)
  {
    first = 1;
    int nintrinsic_sigma;
    int natsrgm;
    double** intrinsic_sigma;
    double** atsrgm;
    double** tmp_alpha;
    double** tmp_sigma;
    double* alpha;
    double* sigma;

    intrinsic_sigma = (double**) malloc(1*sizeof(double*));
    atsrgm = (double**) malloc(1*sizeof(double*));
    tmp_alpha = (double**) malloc(1*sizeof(double*));
    tmp_sigma = (double**) malloc(1*sizeof(double*));
    if (intrinsic_sigma == NULL || atsrgm == NULL || tmp_alpha == NULL || tmp_sigma == NULL)
    {
      log_fatal("fail allocation");
      exit(1);
    }
    else 
    {
      (*intrinsic_sigma) = NULL;
      (*atsrgm) = NULL;
      (*tmp_alpha) = NULL;
      (*tmp_sigma) = NULL; 
    }
   
    setup_SDSS_P_true_lambda_given_mass(&nintrinsic_sigma, intrinsic_sigma, &natsrgm, atsrgm, 
      tmp_alpha, tmp_sigma, 0);

    const gsl_interp2d_type* T = gsl_interp2d_bilinear;
    falpha = gsl_spline2d_alloc(T, nintrinsic_sigma, natsrgm);
    fsigma = gsl_spline2d_alloc(T, nintrinsic_sigma, natsrgm);
    if (falpha == NULL || fsigma == NULL)
    {
      log_fatal("fail allocation");
      exit(1);
    }
    { // we don't want to guess the appropriate GSL z array ordering in z = f(x, y) BEGINS
      alpha = (double*) malloc(nintrinsic_sigma*natsrgm*sizeof(double));
      sigma = (double*) malloc(nintrinsic_sigma*natsrgm*sizeof(double));
      if (alpha == NULL || sigma == NULL)
      {
        log_fatal("fail allocation");
        exit(1);
      }
      for (int i=0; i<nintrinsic_sigma; i++) 
      {
        for (int j=0; j<natsrgm; j++) 
        {
          int status = 0;
          status = gsl_spline2d_set(falpha, alpha, i, j, (*tmp_alpha)[i*nintrinsic_sigma+j]);
          if (status) 
          {
            log_fatal(gsl_strerror(status));
            exit(1);
          }
          status = gsl_spline2d_set(fsigma, sigma, i, j, (*tmp_sigma)[i*nintrinsic_sigma+j]);
          if (status) 
          {
            log_fatal(gsl_strerror(status));
            exit(1);
          }
        }
      }
    } // we don't want to guess the appropriate GSL z array ordering in z = f(x, y) ENDS

    int status = 0;
    status = gsl_spline2d_init(falpha, intrinsic_sigma[0], atsrgm[0], alpha, 
      nintrinsic_sigma, natsrgm);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
    
    status = gsl_spline2d_init(fsigma, intrinsic_sigma[0], atsrgm[0], sigma, 
      nintrinsic_sigma, natsrgm);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }

    free(intrinsic_sigma[0]);
    free(atsrgm[0]);
    free(tmp_alpha[0]);
    free(tmp_sigma[0]);
    free(intrinsic_sigma);
    free(atsrgm);
    free(tmp_alpha);
    free(tmp_sigma);
    free(alpha); // GSL SPLINE 2D copies the array
    free(sigma); // GSL SPLINE 2D copies the array
  }

  const double mass_min = pow(10.0, nuisance.cluster_MOR[0]);
  const double mass_M1 = pow(10.0, nuisance.cluster_MOR[1]);
  const double intrinsic_alpha = nuisance.cluster_MOR[2];
  const double intrinsic_sigma = nuisance.cluster_MOR[3]; //intrisic scatter, mass-richness relation
  
  // atsrgm = average true satellite richness given mass
  const double tmp = (mass - mass_min)/(mass_M1 - mass_min);
  double atsrgm = pow(tmp, intrinsic_alpha)*pow(((1 + z)/1.45), nuisance.cluster_MOR[4]); 
  if (atsrgm > 160) 
  {
    atsrgm = 160;
  }
  else if (atsrgm < 1) 
  {
    atsrgm = 1;
  }
  
  int status = 0;

  double alpha = 0.0;
  status = gsl_spline2d_eval_e(falpha, intrinsic_sigma, atsrgm, NULL, NULL, &alpha);
  if (status) 
  {
    log_fatal(gsl_strerror(status));
    exit(1);
  }
  double sigma = 0.0;
  status = gsl_spline2d_eval_e(fsigma, intrinsic_sigma, atsrgm, NULL, NULL, &sigma);
  if (status) 
  {
    log_fatal(gsl_strerror(status));
    exit(1);
  }

  const double y = 1.0/(M_SQRT2*abs(sigma));
  const double x = (true_lambda - atsrgm)*y;
  const double result1 = exp(-x*x)*y/M_SQRTPI;
  gsl_sf_result result2;
  status = gsl_sf_erfc_e(-alpha*x, &result2);
  if (status) 
  {
    log_fatal(gsl_strerror(status));
    exit(1);
  }  
  return result1*result2.val;
}

// Cocoa: we try to avoid reading of files in the cosmolike_core code 
// Cocoa: (reading is done in the C++/python interface)
void setup_SDSS_P_lambda_obs_given_true_lambda(int* io_nz, double** io_z, int* io_nlambda, 
double** io_lambda, double** io_tau, double** io_mu, double** io_sigma, double** io_fmask, 
double** io_fprj, int io)
{
  static int nz;
  static int nlambda; 
  static double* z = NULL; 
  static double* lambda = NULL;
  static double* tau = NULL;
  static double* mu = NULL;
  static double* sigma = NULL;
  static double* fmask = NULL;
  static double* fprj = NULL;

  if (io == 1) // IO == 1 IMPLES THAT IO_XXX WILL COPIED TO LOCAL XXX
  {
    if (z != NULL)
    {
      free(z);
    }
    else if (lambda != NULL)
    {
      free(lambda);
    }
    else if (tau != NULL)
    {
      free(tau);
    }
    else if (mu != NULL)
    {
      free(mu);
    }
    else if (sigma != NULL)
    {
      free(sigma);
    }
    else if (fmask != NULL)
    {
      free(fmask);
    }
    else if (fprj != NULL)
    {
      free(fprj);
    }

    nz = (*io_nz);
    nlambda = (*io_nlambda);
    if (!(nz > 5) || !(nlambda > 5))
    {
      log_fatal("array to small for 2D interpolation");
      exit(1);
    }

    z =       (double*) malloc(nz*sizeof(double));
    lambda =  (double*) malloc(nlambda*sizeof(double));
    tau =     (double*) malloc(nz*nlambda*sizeof(double));
    mu =      (double*) malloc(nz*nlambda*sizeof(double));
    sigma =   (double*) malloc(nz*nlambda*sizeof(double));
    fmask =   (double*) malloc(nz*nlambda*sizeof(double));
    fprj =    (double*) malloc(nz*nlambda*sizeof(double));
    if (z == NULL || lambda == NULL || tau == NULL || mu == NULL || sigma == NULL || 
      fmask == NULL || fprj == NULL) 
    {
      log_fatal("fail allocation");
      exit(1);
    }
    for (int i=0; i<nz; i++)
    {
      z[i] = (*io_z)[i];
      for (int j=0; j<nlambda; j++)
      {
        if (i == 0) 
        {
          lambda[j] = (*io_lambda)[j];
        }
        tau[i*nz + j] = (*io_tau)[i*nz + j];
        mu[i*nz + j] = (*io_mu)[i*nz + j];
        sigma[i*nz + j] = (*io_sigma)[i*nz + j];
        fmask[i*nz + j] = (*io_fmask)[i*nz + j];
        fprj[i*nz + j] = (*io_fprj)[i*nz + j];
      }
    }
  }
  else
  {
    // IO != 1 IMPLES THAT LOCAL H(Z) WILL BE COPIED TO IO_chi(Z)
    if (z == NULL || lambda == NULL || tau == NULL || mu == NULL || sigma == NULL || fmask == NULL 
        || fprj == NULL || io_lambda == NULL || io_tau == NULL || io_mu == NULL || io_sigma == NULL
        || io_fmask == NULL || io_fprj == NULL)
    {
      log_fatal("array/pointer not allocated");
      exit(1);
    }
    if ((*io_z) != NULL)
    {
      free((*io_z));
      (*io_z) = NULL;
    }
    else if ((*io_lambda) != NULL)
    {
      free((*io_lambda));
      (*io_lambda) = NULL;
    }
    else if ((*io_tau) != NULL)
    {
      free((*io_tau));
      (*io_tau) = NULL;
    }
    else if ((*io_mu) != NULL)
    {
      free((*io_mu));
      (*io_mu) = NULL;
    }
    else if ((*io_sigma) != NULL)
    {
      free((*io_sigma));
      (*io_sigma) = NULL;
    }
    else if ((*io_fmask) != NULL)
    {
      free((*io_fmask));
      (*io_fmask) = NULL;
    }   
    else if ((*io_fprj) != NULL)
    {
      free((*io_fprj));
      (*io_fprj) = NULL;
    }
   
    (*io_z)       = (double*) malloc(nz*sizeof(double));
    (*io_lambda)  = (double*) malloc(nlambda*sizeof(double));
    (*io_tau)     = (double*) malloc(nz*nlambda*sizeof(double));
    (*io_mu)      = (double*) malloc(nz*nlambda*sizeof(double));
    (*io_sigma)   = (double*) malloc(nz*nlambda*sizeof(double));
    (*io_fmask)   = (double*) malloc(nz*nlambda*sizeof(double));
    (*io_fprj)    = (double*) malloc(nz*nlambda*sizeof(double));

    if ((*io_z) == NULL || (*io_lambda) == NULL || (*io_tau) == NULL || (*io_mu) == NULL ||
       (*io_sigma) == NULL || (*io_fmask) == NULL || (*io_fprj) == NULL) 
    {
      log_fatal("fail allocation");
      exit(1);
    }
    for (int i=0; i<nz; i++)
    {
      (*io_z)[i] = z[i];
      for (int j=0; j<nlambda; j++)
      {
        if (i == 0) 
        {
          (*io_lambda)[j] = lambda[j];
        }
        (*io_tau)[i*nz + j] = tau[i*nz + j];
        (*io_mu)[i*nz + j]  = mu[i*nz + j];
        (*io_sigma)[i*nz + j] = sigma[i*nz + j];
        (*io_fmask)[i*nz + j] = fmask[i*nz + j];
        (*io_fprj)[i*nz + j] = fprj[i*nz + j];
      }
    }
  }
}

double SDSS_P_lambda_obs_given_true_lambda(const double observed_lambda, const double true_lambda, 
const double zz) 
{
  static int first = 0;
  static gsl_spline2d* ftau;
  static gsl_spline2d* fmu;
  static gsl_spline2d* fsigma;
  static gsl_spline2d* ffmask;
  static gsl_spline2d* ffprj;

  if (first == 0)
  {
    first = 1;
    int nz;
    int nlambda;
    double** tmp_tau;
    double** tmp_mu;
    double** tmp_sigma;
    double** tmp_fmask;
    double** tmp_fprj;
    double** z;
    double** lambda;
    double* tau;
    double* mu;
    double* sigma;
    double* fmask;
    double* fprj;

    z           = (double**) malloc(1*sizeof(double*));
    lambda      = (double**) malloc(1*sizeof(double*));
    tmp_tau     = (double**) malloc(1*sizeof(double*));
    tmp_mu      = (double**) malloc(1*sizeof(double*));
    tmp_sigma   = (double**) malloc(1*sizeof(double*));
    tmp_fmask   = (double**) malloc(1*sizeof(double*));
    tmp_fprj    = (double**) malloc(1*sizeof(double*));
    if (z == NULL || lambda == NULL || tmp_tau == NULL || tmp_mu == NULL || tmp_sigma == NULL
        || tmp_fmask == NULL || tmp_fprj == NULL)
    {
      log_fatal("fail allocation");
      exit(1);
    }
    else 
    {
      (*z) = NULL;
      (*lambda) = NULL;
      (*tmp_tau) = NULL;
      (*tmp_mu) = NULL; 
      (*tmp_sigma) = NULL; 
      (*tmp_fmask) = NULL; 
      (*tmp_fprj) = NULL; 
    }
   
    setup_SDSS_P_lambda_obs_given_true_lambda(&nz, z, &nlambda, lambda, tmp_tau, tmp_mu, tmp_sigma, 
      tmp_fmask, tmp_fprj, 0);

    const gsl_interp2d_type *T = gsl_interp2d_bilinear;
    ftau = gsl_spline2d_alloc(T, nz, nlambda);
    fmu = gsl_spline2d_alloc(T, nz, nlambda);
    fsigma = gsl_spline2d_alloc(T, nz, nlambda);
    ffmask = gsl_spline2d_alloc(T, nz, nlambda);    
    ffprj = gsl_spline2d_alloc(T, nz, nlambda);
    if (ftau == NULL || fmu == NULL || fsigma == NULL || ffmask == NULL || ffprj == NULL)
    {
      log_fatal("fail allocation");
      exit(1);
    }

    // we don't want to guess the appropriate GSL z array ordering in z = f(x, y)
    tau     = (double*) malloc(nz*nlambda*sizeof(double));
    mu      = (double*) malloc(nz*nlambda*sizeof(double));
    sigma   = (double*) malloc(nz*nlambda*sizeof(double));
    fmask   = (double*) malloc(nz*nlambda*sizeof(double));
    fprj    = (double*) malloc(nz*nlambda*sizeof(double));
    if (tau == NULL || mu == NULL || sigma == NULL || fmask == NULL || fprj == NULL)
    {
      log_fatal("fail allocation");
      exit(1);
    }

    for (int i=0; i<nz; i++) 
    {
      for (int j=0; j<nlambda; j++) 
      {
        int status = 0;
        status = gsl_spline2d_set(ftau, tau, i, j, (*tmp_tau)[i*nz + j]);
        if (status) 
        {
          log_fatal(gsl_strerror(status));
          exit(1);
        }
        status = gsl_spline2d_set(fmu, mu, i, j, (*tmp_mu)[i*nz + j]);
        if (status) 
        {
          log_fatal(gsl_strerror(status));
          exit(1);
        }
        status = gsl_spline2d_set(fsigma, sigma, i, j, (*tmp_sigma)[i*nz + j]);
        if (status) 
        {
          log_fatal(gsl_strerror(status));
          exit(1);
        }
        status = gsl_spline2d_set(ffmask, fmask, i, j, (*tmp_fmask)[i*nz + j]);
        if (status) 
        {
          log_fatal(gsl_strerror(status));
          exit(1);
        }
        status = gsl_spline2d_set(ffprj, fprj, i, j, (*tmp_fprj)[i*nz + j]);
        if (status) 
        {
          log_fatal(gsl_strerror(status));
          exit(1);
        }
      }
    }
    int status = 0;
    status = gsl_spline2d_init(ftau, (*z), (*lambda), tau, nz, nlambda);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
    status = gsl_spline2d_init(fmu, (*z), (*lambda), mu, nz, nlambda);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
    status = gsl_spline2d_init(fsigma, (*z), (*lambda), sigma, nz, nlambda);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
    status = gsl_spline2d_init(ffmask, (*z), (*lambda), fmask, nz, nlambda);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
    status = gsl_spline2d_init(ffprj, (*z), (*lambda), fprj, nz, nlambda);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
    free(z[0]);
    free(lambda[0]);
    free(tmp_tau[0]);
    free(tmp_mu[0]);
    free(tmp_sigma[0]);
    free(tmp_fmask[0]);
    free(tmp_fprj[0]);
    free(z);
    free(lambda);
    free(tmp_tau);
    free(tmp_mu);
    free(tmp_sigma);
    free(tmp_fmask);
    free(tmp_fprj);
    free(tau);   // GSL SPLINE 2D copies the array
    free(mu);    // GSL SPLINE 2D copies the array
    free(sigma); // GSL SPLINE 2D copies the array
    free(fmask); // GSL SPLINE 2D copies the array
    free(fprj);  // GSL SPLINE 2D copies the array
  }
  
  double tau, mu, sigma, fmask, fprj;
  {
    int status = 0;
    status = gsl_spline2d_eval_e(ftau, zz, true_lambda, NULL, NULL, &tau);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }  
    status = gsl_spline2d_eval_e(fmu, zz, true_lambda, NULL, NULL, &mu);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }  
    status = gsl_spline2d_eval_e(fsigma, zz, true_lambda, NULL, NULL, &sigma);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
    status = gsl_spline2d_eval_e(ffmask, zz, true_lambda, NULL, NULL, &fmask);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
    status = gsl_spline2d_eval_e(ffprj, zz, true_lambda, NULL, NULL, &fprj);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
  }

  const double x = 1.0/(M_SQRT2*abs(sigma));
  const double y = 1.0/true_lambda;
  const double j = exp(0.5*tau*(2*mu+tau*sigma*sigma-2*observed_lambda));

  double r0 = exp(-1.*(observed_lambda-mu)*(observed_lambda-mu)*x*x);
  r0 *= (1-fmask)*(1-fprj)*x/M_SQRTPI;

  double r1 = 0.5*((1-fmask)*fprj*tau + fmask*fprj*y)*j;
  r1 *= gsl_sf_erfc((mu+tau*sigma*sigma-observed_lambda)*x);

  double r2 = 0.5*fmask*y;
  r2 *= gsl_sf_erfc((mu-observed_lambda-true_lambda)*x) - gsl_sf_erfc((mu-observed_lambda)*x);

  double r3 = 0.5*fmask*fprj*y*exp(-tau*true_lambda)*j;
  r3 *= gsl_sf_erfc((mu+tau*sigma*sigma-observed_lambda-true_lambda)*x);
  
  return r0 + r1 + r2 - r3;
}

static double SDSS_P_lambda_obs_lambda_true_given_M(double true_lambda, 
const double observed_lambda, const double M, const double z)
{
  const double r1 = SDSS_P_lambda_obs_given_true_lambda(observed_lambda, true_lambda, z);
  const double r2 = SDSS_P_true_lambda_given_mass(true_lambda, M, z);
  return r1*r2;
}

double SDSS_P_lambda_obs_lambda_true_given_M_wrapper(double true_lambda, void* params)
{
  double* ar = (double*) params;
  const double M = ar[0];
  const double z = ar[1];
  const double observed_lambda = ar[2];
  return SDSS_P_lambda_obs_lambda_true_given_M(true_lambda, observed_lambda, M, z);
}


double SDSS_P_lambda_obs_given_M(const double observed_lambda, const double M, const double z, 
const int init_static_vars_only)
{  
  double params_in[3] = {M, z, observed_lambda}; 
  
  const double true_lambda_min = limits.SDSS_P_lambda_obs_given_M_true_lambda_min;
  const double true_lambda_max = limits.SDSS_P_lambda_obs_given_M_true_lambda_max;
  
  return (init_static_vars_only == 1) ? 
    SDSS_P_lambda_obs_lambda_true_given_M_wrapper(true_lambda_min, (void*) params_in) :
    int_gsl_integrate_medium_precision(SDSS_P_lambda_obs_lambda_true_given_M_wrapper, 
      (void*) params_in, true_lambda_min, true_lambda_max, NULL, GSL_WORKSPACE_SIZE);
}

static double SDSS_P_lambda_obs_given_M_wrapper(double observed_lambda, void* params)
{
  double* ar = (double*) params; 
  const double M = ar[0];
  const double z = ar[1];
  const int init_static_vars_only = (int) ar[2];
  return SDSS_P_lambda_obs_given_M(observed_lambda, M, z, init_static_vars_only);
}

double SDSS_binned_P_lambda_obs_given_M(const int nl, const double M, const double z, 
const int init_static_vars_only)
{
  double params[3] = {M, z, (double) init_static_vars_only};
  
  const int nl_min = Cluster.N_min[nl];
  const int nl_max = Cluster.N_max[nl];
  
  return (init_static_vars_only == 1) ?
    SDSS_P_lambda_obs_given_M_wrapper((double) nl_min, (void*) params) :
    int_gsl_integrate_medium_precision(SDSS_P_lambda_obs_given_M_wrapper, (void*) params, 
      (double) nl_min, (double) nl_max, NULL, GSL_WORKSPACE_SIZE);
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// INTERFACE - binned_P_lambda_obs_given_M
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

// \int_(bin_lambda_obs_min)^(bin_lambda_obs_max) \dlambda_obs P(\lambda_obs|M)
// (see for example https://arxiv.org/pdf/1810.09456.pdf - eq 3 and eq 6) 
double binned_P_lambda_obs_given_M_nointerp(const int nl, const double M, const double z, 
const int init_static_vars_only)
{
  if (strcmp(Cluster.model, "SDSS") == 0) 
  {
    return SDSS_binned_P_lambda_obs_given_M(nl, M, z, init_static_vars_only);
  }
  else if (strcmp(Cluster.model, "BUZZARD") == 0)
  {
    return buzzard_binned_P_lambda_obs_given_M(nl, M, z, init_static_vars_only);
  }
  else
  {
    log_fatal("Cluster.model not implemented");
    exit(1);
  }
}

double binned_P_lambda_obs_given_M(const int nl, const double M, const double z) 
{
  static cosmopara C;
  static nuisancepara N;
  static double*** table = 0;

  const int N_l = Cluster.N200_Nbin;
  
  const int N_M  = Ntable.binned_P_lambda_obs_given_M_size_M_table;
  const double log_M_min = limits.cluster_util_log_M_min;
  const double log_M_max = limits.cluster_util_log_M_max;
  const double dlogM = (log_M_max - log_M_min)/((double) N_M - 1.0);
  
  const int N_z = Ntable.binned_P_lambda_obs_given_M_size_z_table;
  const double zmin = limits.binned_P_lambda_obs_given_M_zmin_table; 
  const double zmax = limits.binned_P_lambda_obs_given_M_zmax_table;
  const double dz = (zmax - zmin)/((double) N_z - 1.0);
  
  if (table == 0)
  {
    table = (double***) malloc(sizeof(double**)*N_l);
    for(int i=0; i<N_l; i++)
    {
      table[i] = (double**) malloc(sizeof(double*)*N_z);
      for(int j=0; j<N_z; i++)
      {
        table[i][j] = (double*) malloc(sizeof(double)*N_M);
      }
    }
  }
  if (recompute_clusters(C, N))
  { 
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      const double MM = pow(10.0, limits.cluster_util_log_M_min);
      double init_static_vars_only = binned_P_lambda_obs_given_M_nointerp(0, MM, zmin, 1);
    }
    #pragma GCC diagnostic pop
    #pragma omp parallel for collapse(3)
    for (int i=0; i<N_l; i++) 
    {
      for (int j=0; j<N_z; j++) 
      {
        for (int k=0; k<N_M; k++) 
        {
          const double zz = zmin + j*dz;
          const double MM = pow(10.0, limits.cluster_util_log_M_min + k*dlogM);
          table[i][j][k] = binned_P_lambda_obs_given_M_nointerp(i, MM, zz, 0);
        }
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
  if (z < zmin || z > zmax)
  {
    log_fatal("z = %e outside look-up table range [%e,%e]", z, zmin, zmax);
    exit(1);
  } 
  const double logM = log10(M);
  if (logM < limits.cluster_util_log_M_min || logM > limits.cluster_util_log_M_max)
  {
    log_fatal("logM = %e outside look-up table range [%e,%e]", logM, 
      limits.cluster_util_log_M_min, limits.cluster_util_log_M_max);
    exit(1);
  } 
  return interpol2d(table[nl], N_z, zmin, zmax, dz, z, N_M, log_M_min, log_M_max, dlogM,
    logM, 1.0, 1.0);
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// BINNED CLUSTER MASS FUNCTION 
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

double dndlnM_times_binned_P_lambda_obs_given_M(double lnM, void* params)
{
  double* ar = (double*) params; 
  const int nl = (int) ar[0];
  if (nl < 0 || nl > Cluster.N200_Nbin - 1)
  {
    log_fatal("invalid bin input nl = %d", nl);
    exit(1);
  }
  const double z = ar[1];
  if (!(z>0))
  {
    log_fatal("invalid redshift input z = %d", z);
    exit(1);
  } 
  const double a = 1.0/(1.0 + z);
  const double M = exp(lnM);

  double mfunc; 
  if (Cluster.hmf_model == 0)
  {
    mfunc = massfunc(M, a);
  }
  else if (Cluster.hmf_model == 1)
  {
    mfunc = tinker_emulator_massfunc(M, a);
  }
  else
  {
    log_fatal("massfunc model %i not implemented", Cluster.hmf_model);
    exit(1); 
  }
  return mfunc*M*binned_P_lambda_obs_given_M(nl, M, z);
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// CLUSTER BIAS
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

double BSF(const double M)
{
  double BS;
  if (Cluster.N_SF == 1) 
  {
    BS = nuisance.cluster_selection[0]; 
  }
  else if (Cluster.N_SF == 2)
  {
    BS = 
    nuisance.cluster_selection[0]*pow((M/M_PIVOT), nuisance.cluster_selection[1]);
  }
  else
  {
    log_fatal("Cluster selection bias model %i not implemented", Cluster.bias_model);
    exit(0); 
  }
  return BS;
}

double B1_x_BSF(const double M, const double a)
{ // cluster bias including selection bias
  double tmp_B1;
  if (Cluster.bias_model == 0) 
  {
    tmp_B1 = B1(M, a);
  } 
  else if (Cluster.bias_model == 1) 
  {
    tmp_B1 = tinker_emulator_B1(M, a);
  }
  else 
  {
    log_fatal("Cluster bias model %i not implemented", Cluster.bias_model);
    exit(0); 
  }
  return tmp_B1 * BSF(M);
}

double B2_x_BSF(const double M, const double a)
{ // cluster bias (b2) including selection bias
  double tmp_B1;
  if (Cluster.bias_model == 0) 
  {
    tmp_B1 = B1(M, a);
  } 
  else if (Cluster.bias_model == 1) 
  {
    tmp_B1 = tinker_emulator_B1(M, a);
  }
  else 
  {
    log_fatal("Cluster bias model %i not implemented", Cluster.bias_model);
    exit(0); 
  }
  return b2_from_b1(tmp_B1) * BSF(M);
}

double B1M1_x_BSF(const double M, const double a)
{
  double B1M1;
  if (Cluster.bias_model == 0) 
  {
    B1M1 = B1(M, a) - 1;
  } 
  else if (Cluster.bias_model == 1) 
  {
    B1M1 = tinker_emulator_B1(M, a) - 1;
  }
  else 
  {
    log_fatal("Cluster bias model %i not implemented", Cluster.bias_model);
    exit(0); 
  }
  return B1M1 * BSF(M);
}

double int_for_weighted_B1(double lnM, void* params)
{
  double* ar = (double*) params; // {nl, z}
  const double z = ar[1];
  if (!(z>0))
  {
    log_fatal("invalid redshift input z = %d", z);
    exit(1);
  } 
  const double M = exp(lnM);  
  const double a = 1.0/(1.0 + z);
  return B1_x_BSF(M, a)*dndlnM_times_binned_P_lambda_obs_given_M(lnM, params); 
}

double weighted_B1_nointerp(const int nl, const double z, const int init_static_vars_only)
{  
  const double ln_M_min = limits.cluster_util_log_M_min/M_LOG10E;
  const double ln_M_max = limits.cluster_util_log_M_max/M_LOG10E;
 
  double param[2] = {(double) nl, z};

  if (init_static_vars_only == 1)
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      const double r1 = int_for_weighted_B1(ln_M_min, (void*) param);
      const double r2 = dndlnM_times_binned_P_lambda_obs_given_M(ln_M_min, (void*) param);
      return 0.0;
    }
    #pragma GCC diagnostic pop
  }
  else
  {
    const double r1 = int_gsl_integrate_low_precision(int_for_weighted_B1, 
      (void*) param, ln_M_min, ln_M_max, NULL, GSL_WORKSPACE_SIZE);
  
    const double r2 = int_gsl_integrate_low_precision(dndlnM_times_binned_P_lambda_obs_given_M, 
      (void*) param, ln_M_min, ln_M_max, NULL, GSL_WORKSPACE_SIZE);

    return (r2 <= 0) ? 0.0 : r1/r2;
  }
}

double weighted_B1(const int nl, const double z)
{
  static cosmopara C;
  static nuisancepara N;
  static double** table;

  const int N_l = Cluster.N200_Nbin;
  const int N_a = Ntable.N_a;
  const double zmin = fmax(tomo.cluster_zmin[0] - 0.05, 0.01); 
  const double zmax = tomo.cluster_zmax[tomo.cluster_Nbin - 1] + 0.05;
  const double amin = 1.0/(1.0 + zmax); 
  const double amax = 1.0/(1.0 + zmin);
  const double da = (amax - amin)/((double) N_a - 1.);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*N_l);
    for(int i=0; i<N_l; i++)
    {
      table[i] = (double*) malloc(sizeof(double)*N_a);
    }
  }
  if (recompute_clusters(C, N))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"    
    {
      double init_static_vars_only = weighted_B1_nointerp(0, 1.0/amin - 1.0, 1); 
    }
    #pragma GCC diagnostic pop
    #pragma omp parallel for collapse(2)
    for (int n=0; n<N_l; n++)
    {
      for (int i=0; i<N_a; i++)
      {
        const double a = amin + i*da;
        table[n][i] = weighted_B1_nointerp(n, 1.0/a - 1.0, 0);
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  const double a = 1.0/(z + 1.0);
  if (a < amin || a > amax)
  {
    log_fatal("a = %e outside look-up table range [%e,%e]", a, amin, amax);
    exit(1);
  } 
  if (nl < 0 || nl > N_l - 1)
  {
    log_fatal("invalid bin input nl = %d", nl);
    exit(1);
  }
  return interpol(table[nl], N_a, amin, amax, da, a, 0., 0.);
}

double int_for_weighted_B2(double lnM, void* params)
{
  double* ar = (double*) params; //nl, z
  const double z = ar[1];
  if (!(z>0))
  {
    log_fatal("invalid redshift input z = %d", z);
    exit(1);
  }
  const double M = exp(lnM);  
  const double a = 1./(1. + z);
  return B2_x_BSF(M, a)*dndlnM_times_binned_P_lambda_obs_given_M(lnM, params); 
}

double weighted_B2_nointerp(const int nl, const double z, const int init_static_vars_only)
{  
  const double ln_M_min = limits.cluster_util_log_M_min/M_LOG10E;
  const double ln_M_max = limits.cluster_util_log_M_max/M_LOG10E;
 
  double param[2] = {(double) nl, z};

  if (init_static_vars_only == 1)
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      const double r1 = int_for_weighted_B2(ln_M_min, (void*) param);
      const double r2 = dndlnM_times_binned_P_lambda_obs_given_M(ln_M_min, (void*) param);
      return 0.0;
    }
    #pragma GCC diagnostic pop
  }
  else
  {
    const double r1 = int_gsl_integrate_low_precision(int_for_weighted_B2, 
      (void*) param, ln_M_min, ln_M_max, NULL, GSL_WORKSPACE_SIZE);

    const double r2 = int_gsl_integrate_low_precision(dndlnM_times_binned_P_lambda_obs_given_M, 
      (void*) param, ln_M_min, ln_M_max, NULL, GSL_WORKSPACE_SIZE);

    return (r2 <= 0) ? 0.0 : r1/r2; 
  }
}

double weighted_B2(const int nl, const double z)
{
  static cosmopara C;
  static nuisancepara N;
  static double** table;

  const int N_l = Cluster.N200_Nbin;
  const int N_a = Ntable.N_a;
  const double zmin = fmax(tomo.cluster_zmin[0] - 0.05, 0.01); 
  const double zmax = tomo.cluster_zmax[tomo.cluster_Nbin - 1] + 0.05;
  const double amin = 1.0/(1.0 + zmax); 
  const double amax = 1.0/(1.0 + zmin);
  const double da = (amax-amin)/((double) N_a - 1.0);

  if (table == 0) 
  {
    table = (double**) malloc(sizeof(double*)*N_l);
    for(int n=0; n<N_l; n++)
    {
      table[n] = (double*) malloc(sizeof(double)*N_a);
    }
  }
  if (recompute_clusters(C, N))
  {  
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"    
    {
      double init_static_vars_only =  weighted_B2_nointerp(0, 1.0/amin - 1.0, 1);
    }
    #pragma GCC diagnostic pop
    #pragma omp parallel for collapse(2)
    for (int n=0; n<N_l; n++)
    {
      for (int i=0; i<N_a; i++)
      {
        const double a = amin + i*da;
        table[n][i] = weighted_B2_nointerp(n, 1.0/a - 1.0, 0);
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }

  const double a = 1.0/(z + 1.0);
  if (a < amin || a > amax)
  {
    log_fatal("a = %e outside look-up table range [%e,%e]", a, amin, amax);
    exit(1);
  } 
  if (nl < 0 || nl > N_l - 1)
  {
    log_fatal("invalid bin input nl1 = %d", nl);
    exit(1);
  }
  return interpol(table[nl], N_a, amin, amax, da, a, 0., 0.);
}

double int_for_weighted_B1M1(double lnM, void* params)
{
  double* ar = (double*) params;
  const double z = ar[1];
  if (!(z>0))
  {
    log_fatal("invalid redshift input z = %d", z);
    exit(1);
  }
  const double M = exp(lnM);  
  const double a = 1.0/(1.0 + z);
  return B1M1_x_BSF(M, a)*dndlnM_times_binned_P_lambda_obs_given_M(lnM, params);
}

double weighted_B1M1_nointerp(const int nl, const double z, const int init_static_vars_only)
{  
  const double ln_M_min = limits.cluster_util_log_M_min/M_LOG10E;
  const double ln_M_max = limits.cluster_util_log_M_max/M_LOG10E;
  
  double params[2] = {(double) nl, z};

  if (init_static_vars_only == 1)
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      const double r1 = int_for_weighted_B1M1(ln_M_min, (void*) params);
      const double r2 = dndlnM_times_binned_P_lambda_obs_given_M(ln_M_min, (void*) params);
      return 0.0;
    }
    #pragma GCC diagnostic pop
  }
  else
  {
    const double r1 = int_gsl_integrate_low_precision(int_for_weighted_B1M1, 
      (void*) params, ln_M_min, ln_M_max, NULL, GSL_WORKSPACE_SIZE);

    const double r2 = int_gsl_integrate_low_precision(dndlnM_times_binned_P_lambda_obs_given_M, 
      (void*) params, ln_M_min, ln_M_max, NULL, GSL_WORKSPACE_SIZE);
    
    return (r2 <= 0) ? 0.0 : r1/r2; 
  }
}

double weighted_B1M1(const int nl, const double z)
{
  static cosmopara C;
  static nuisancepara N;
  static double** table;

  const int N_l = Cluster.N200_Nbin; 
  const int N_a = Ntable.N_a;
  const double zmin = fmax(tomo.cluster_zmin[0] - 0.05, 0.01); 
  const double zmax = tomo.cluster_zmax[tomo.cluster_Nbin - 1] + 0.05;
  const double amin = 1./(1.0 + zmax); 
  const double amax = 1./(1.0 + zmin);
  const double da = (amax - amin)/((double) N_a - 1.0);

  if (table == 0) 
  {
    table = (double**) malloc(sizeof(double*)*N_l);
    for(int i=0; i<N_l; i++)
    {
      table[i] = (double*) malloc(sizeof(double)*N_a);
    }
  }
  if (recompute_clusters(C, N))
  { 
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"    
    {
      double init_static_vars_only =  weighted_B1M1_nointerp(0, 1.0/amin - 1.0, 1);
    }
    #pragma GCC diagnostic pop   
    #pragma omp parallel for collapse(2)
    for (int n=0; n<N_l; n++)
    {
      for (int i=0; i<N_a; i++)
      {
        const double a = amin + i*da;
        table[n][i] = weighted_B1M1_nointerp(n, 1.0/a - 1.0, 0);
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  const double a = 1.0/(z + 1.0);
  if (a < amin || a > amax)
  {
    log_fatal("a = %e outside look-up table range [%e,%e]", a, amin, amax);
    exit(1);
  } 
  if (nl < 0 || nl > N_l - 1)
  {
    log_fatal("invalid bin input nl = %d", nl);
    exit(1);
  }
  return interpol(table[nl], N_a, amin, amax, da, a, 0., 0.);
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// CLUSTER POWER SPECTRUM
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// nl{1,2} = lambda_obs bins, n{i,j} = cluster redshift bins

double binned_p_cc(const double k, const double a, const int nl1, const int nl2, 
  const int use_linear_ps)
{ // binned in lambda_obs (nl1, nl2 = lambda_obs bin (cluster 1 and 2))
  if (!(a>0) || !(a<1))
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  if (nl1 < 0 || nl1 > Cluster.N200_Nbin - 1 || nl2 < 0 || nl2 > Cluster.N200_Nbin - 1)
  {
    log_fatal("invalid bin input (nl1,nl2) = (%d,%d)", nl1, nl2);
    exit(1);
  }

  const double z = 1.0/a - 1.0;
  const double cluster_B1_1 = weighted_B1(nl1, z);
  const double cluster_B1_2 = (nl1 == nl2) ? cluster_B1_1 : weighted_B1(nl2, z);  
  const double PK = use_linear_ps == 1 ? p_lin(k, a) : Pdelta(k, a);

  double P_1loop = 0;
  if (has_b2_galaxies() || (Cluster.nonlinear_bias > 0 && use_linear_ps != 1))
  {
    const double g = growfac(a);
    const double g4 = g*g*g*g;
    
    const double b1g = cluster_B1_1;
    const double b2g = weighted_B2(nl1, z);
    const double bs2g = bs2_from_b1(b1g);
    
    const double b1c = cluster_B1_2;
    const double b2c = (nl1 == nl2) ? b2g : weighted_B2(nl2, z);
    const double bs2c = (nl1 == nl2) ? bs2g : bs2_from_b1(b1c);
    
    P_1loop = 0.5*(b1c*b2g+b2c*b1g)*PT_d1d2(k) + 0.25*b2g*b2c*PT_d2d2(k) 
      + 0.5*(b1c*bs2g+b1g*bs2c)*PT_d1s2(k) +0.25*(b2c*bs2g+b2g*bs2c)*PT_d2s2(k) 
      + 0.25*(bs2g*bs2c)*PT_s2s2(k)+ 0.5*(b1c*b3nl_from_b1(b1g) + b1g*b3nl_from_b1(b1c))*PT_d1d3(k);
    P_1loop *= g4;
  }
  return cluster_B1_1*cluster_B1_2*PK + P_1loop;
}

double int_for_int_for_binned_p_cc_incl_halo_exclusion(double lnM1, void* params)
{
  static cosmopara C;
  static nuisancepara N;
  static double***** table = 0; // table -> PK_HALO_EXCL interpolation

  const int N_l = Cluster.N200_Nbin;

  const int N_k = Ntable.N_k_halo_exclusion;
  const double ln_k_min = log(limits.halo_exclusion_k_min); 
  const double ln_k_max = log(limits.halo_exclusion_k_max);
  const double dlnk = (ln_k_max - ln_k_min)/((double) N_k - 1.0);
  
  const int N_a = Ntable.N_a_halo_exclusion;
  const double amin = 1.0/(1.0 + tomo.cluster_zmax[tomo.cluster_Nbin - 1]);
  const double amax = 1.0/(1.0 + tomo.cluster_zmin[0]) + 0.01;
  const double da = (amax - amin)/((double) N_a - 1.0);

  const int N_k_hankel = Ntable.N_k_hankel_halo_exclusion;
  const double ln_k_min_hankel = log(limits.halo_exclusion_k_min_hankel);
  const double ln_k_max_hankel = log(limits.halo_exclusion_k_max_hankel);
  const double dlnk_hankel = (ln_k_max_hankel - ln_k_min_hankel)/((double) N_k_hankel - 1.0);

  const int N_R = Ntable.N_R_halo_exclusion;
  const double R_min = limits.halo_exclusion_R_min;
  const double R_max = limits.halo_exclusion_R_max;
  const double dR = (R_max - R_min)/((double) N_R - 1.0);

  typedef fftw_complex fftwZ;
  double arg[2];
  arg[0] = 0;     // bias
  arg[1] = 0.5;   // order of Bessel function 

  if (table == 0)
  {
    table = (double*****) malloc(sizeof(double****)*N_l);
    for(int i=0; i<N_l; i++)
    {          
      table[i] = (double****) malloc(sizeof(double***)*N_l);
      for(int j=0; j<N_l; j++)
      {
        table[i][j] = (double***) malloc(sizeof(double**)*N_a);
        for (int k=0; k<N_a; k++) 
        { 
          table[i][j][k] = (double**) malloc(sizeof(double*)*N_k);
          for (int p=0; p<N_k; p++) 
          {
            table[i][j][k][p] = (double*) malloc(sizeof(double)*N_R);
          } 
        }
      }
    }
  }
  if (recompute_clusters(C, N))
  { 
    // ---------------------------------------------------------------------------------------------
    // Step 1 - memory allocation (except for FFTW plans)
    // ---------------------------------------------------------------------------------------------    
    fftw_plan**** plan  = (fftw_plan****) malloc(sizeof(fftw_plan***)*N_l);
    fftw_plan**** plan1 = (fftw_plan****) malloc(sizeof(fftw_plan***)*N_l);
    fftwZ***** flP      = (fftwZ*****)    malloc(sizeof(fftwZ****)*N_l);
    fftwZ***** conv     = (fftwZ*****)    malloc(sizeof(fftwZ****)*N_l);
    double***** lP      = (double*****)   malloc(sizeof(double****)*N_l);
    for(int i=0; i<N_l; i++)
    {
      plan[i]  = (fftw_plan***) malloc(sizeof(fftw_plan**)*N_l);
      plan1[i] = (fftw_plan***) malloc(sizeof(fftw_plan**)*N_l);
      flP[i]   = (fftwZ****)    malloc(sizeof(fftwZ***)*N_l);
      conv[i]  = (fftwZ****)    malloc(sizeof(fftwZ***)*N_l);
      lP[i]    = (double****)   malloc(sizeof(double***)*N_l);
      for(int j=0; j<N_l; j++)
      {
        plan[i][j]  = (fftw_plan**) malloc(sizeof(fftw_plan*)*N_a);
        plan1[i][j] = (fftw_plan**) malloc(sizeof(fftw_plan*)*N_a);
        flP[i][j]   = (fftwZ***)    malloc(sizeof(fftwZ**)*N_a);
        conv[i][j]  = (fftwZ***)    malloc(sizeof(fftwZ**)*N_a);
        lP[i][j]    = (double***)   malloc(sizeof(double**)*N_a);
        for (int k=0; k<N_a; k++) 
        { 
          plan[i][j][k]  = (fftw_plan*) malloc(sizeof(fftw_plan)*N_R);
          plan1[i][j][k] = (fftw_plan*) malloc(sizeof(fftw_plan)*N_R);    
          flP[i][j][k]   = (fftwZ**)    malloc(sizeof(fftwZ*)*N_R);
          conv[i][j][k]  = (fftwZ**)    malloc(sizeof(fftwZ*)*N_R);
          lP[i][j][k]    = (double**)   malloc(sizeof(double*)*N_R);
          for (int p=0; p<N_R; p++) 
          {
            flP[i][j][k][p]  = (fftwZ*)  malloc(sizeof(fftwZ)*N_k_hankel);
            conv[i][j][k][p] = (fftwZ*)  malloc(sizeof(fftwZ)*N_k_hankel);
            lP[i][j][k][p]   = (double*) malloc(sizeof(double)*N_k_hankel);
          }
        }
      }
    }
    // ---------------------------------------------------------------------------------------------
    // STEP 2: - pk_to_xi
    // ---------------------------------------------------------------------------------------------
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init_static_vars_only  = weighted_B1(0, 1.0/amin - 1.0);     
      double init_static_vars_only2 = Pdelta(exp(ln_k_min_hankel), amin); 
    }
    #pragma GCC diagnostic pop
    #pragma omp parallel for collapse(3)
    for(int i=0; i<N_l; i++)
    {
      for (int q=0; q<N_k_hankel; q++) 
      {
        for (int k=0; k<N_a; k++) 
        { 
          const double aa = amin + k*da;
          const double zz = 1.0/aa - 1.0;
          const double kk = exp(ln_k_min_hankel + q*dlnk_hankel);
          const double B1_1 = weighted_B1(i, zz);
          const double PK = Pdelta(kk, aa);
          for(int j=i; j<N_l; j++)
          {  
            const double B1_2 = (i == j) ? B1_1 : weighted_B1(j, zz);
            lP[i][j][k][0][q] = kk*sqrt(kk)*PK*B1_1*B1_2;
          }
        }
      }
    }
    for(int i=0; i<N_l; i++)
    {
      for(int j=i; j<N_l; j++)
      {
        for (int k=0; k<N_a; k++) 
        { 
          plan[i][j][k][0] = fftw_plan_dft_r2c_1d(N_k_hankel, lP[i][j][k][0], flP[i][j][k][0], FFTW_ESTIMATE);
        }
      }
    }
    #pragma omp parallel for collapse(2)
    for(int i=0; i<N_l; i++)
    {
      for (int k=0; k<N_a; k++) 
      {
        for(int j=i; j<N_l; j++)
        {
          fftw_execute(plan[i][j][k][0]);
        }
      }
    }  
    for(int i=0; i<N_l; i++)
    {
      for(int j=i; j<N_l; j++)
      {
        for (int k=0; k<N_a; k++) 
        { 
          fftw_destroy_plan(plan[i][j][k][0]);
          plan1[i][j][k][0] = fftw_plan_dft_c2r_1d(N_k_hankel, conv[i][j][k][0], lP[i][j][k][0], FFTW_ESTIMATE);
        }
      }
    }  
    #pragma omp parallel for collapse(2)
    for(int i=0; i<N_l; i++)
    {
      for (int k=0; k<N_a; k++) 
      {
        for(int j=i; j<N_l; j++)
        {    
          for(int q=0; q<N_k_hankel/2+1; q++) 
          {            
            fftw_complex krn;
            hankel_kernel_FT_3D(2.0*M_PI*q/((double) N_k_hankel*dlnk_hankel), &krn, arg, 2);
            conv[i][j][k][0][q][0] = flP[i][j][k][0][q][0]*krn[0] - flP[i][j][k][0][q][1]*krn[1];
            conv[i][j][k][0][q][1] = flP[i][j][k][0][q][1]*krn[0] + flP[i][j][k][0][q][0]*krn[1];
          }

          conv[i][j][k][0][0][1] = 0;            
          conv[i][j][k][0][N_k_hankel/2][1] = 0;

          fftw_execute(plan1[i][j][k][0]);

          for (int q=0; q<N_k_hankel; q++) 
          {
            const double kk = exp(ln_k_min_hankel + q*dlnk_hankel);    
            const double rr = 1.0/kk;
            lP[i][j][k][0][q] = pow(2.0*M_PI*rr, -1.5)*lP[i][j][k][0][q]/((double) N_k_hankel); 
          }
        }
      }
    }
    for(int i=0; i<N_l; i++)
    {
      for (int k=0; k<N_a; k++) 
      {
        for(int j=i; j<N_l; j++)
        {
          fftw_destroy_plan(plan1[i][j][k][0]);
        }
      }
    }
    #pragma omp parallel for collapse(4)
    for(int i=0; i<N_l; i++)
    {
      for (int p=1; p<N_R; p++) 
      { 
        for (int k=0; k<N_a; k++) 
        {
          for (int q=0; q<N_k_hankel; q++)
          {
            for(int j=i; j<N_l; j++)
            {
              lP[i][j][k][p][q] = lP[i][j][k][0][q]; // pk_to_xi does not depend on R!
            }
          }
        }
      }
    }
    // ---------------------------------------------------------------------------------------------
    // Step 3: exclusion_filter
    // ---------------------------------------------------------------------------------------------
    #pragma omp parallel for collapse(4)
    for(int i=0; i<N_l; i++)
    {
      for (int q=0; q<N_k_hankel; q++) 
      {
        for (int k=0; k<N_a; k++) 
        {
          for (int p=0; p<N_R; p++) 
          {
            const double kk = exp(ln_k_min_hankel + q*dlnk_hankel);     
            const double rr = 1.0/kk;
            const double R = R_min + p*dR;

            if (Cluster.halo_exclusion_model == 0)
            {
              if (rr < R)
              {
                for(int j=i; j<N_l; j++)
                {
                  lP[i][j][k][p][q] = -1; 
                }
              } 
            }
            else if (Cluster.halo_exclusion_model == 1)
            { // Baldauf 2013 eq. 43
              if (R > 0) 
              {                     
                gsl_sf_result tmp;
                int status = gsl_sf_erf_e(log(rr/R)/(M_SQRT2*0.0387), &tmp);
                if (status)
                {
                  log_fatal(gsl_strerror(status));
                  exit(1);
                }
                for(int j=i; j<N_l; j++)
                {
                  lP[i][j][k][p][q] = 0.5*(1.0 + tmp.val)*(lP[i][j][k][p][q] + 1.0) - 1.0;
                }
              }
            }
            else
            {  
              log_fatal("halo exclusion model: %d  not implemented", Cluster.halo_exclusion_model);
              exit(1);
            }
          }
        }
      }
    }
    // ---------------------------------------------------------------------------------------------
    // Step4 4: - xi_to_pk
    // ---------------------------------------------------------------------------------------------
    #pragma omp parallel for collapse(4)
    for(int i=0; i<N_l; i++)
    {
      for (int q=0; q<N_k_hankel; q++) 
      {
        for (int k=0; k<N_a; k++) 
        { 
          for (int p=0; p<N_R; p++) 
          {
            const double kk = exp(ln_k_min_hankel + q*dlnk_hankel);    
            const double rr = 1.0/kk;
            for(int j=i; j<N_l; j++)
            {
              lP[i][j][k][p][q] = rr*sqrt(rr)*lP[i][j][k][p][q];
            }
          }
        }
      }
    }
    for(int i=0; i<N_l; i++)
    {
      for(int j=i; j<N_l; j++)
      {
        for (int k=0; k<N_a; k++) 
        { 
          for (int p=0; p<N_R; p++) 
          {
            plan[i][j][k][p] = fftw_plan_dft_r2c_1d(N_k_hankel, lP[i][j][k][p], flP[i][j][k][p], FFTW_ESTIMATE);
          }
        }
      }
    }
    #pragma omp parallel for collapse(3)
    for(int i=0; i<N_l; i++)
    {
      for (int p=0; p<N_R; p++) 
      {
        for (int k=0; k<N_a; k++) 
        {
          for(int j=i; j<N_l; j++)
          {    
            fftw_execute(plan[i][j][k][p]);
          }
        }
      }
    }
    for(int i=0; i<N_l; i++)
    {
      for(int j=i; j<N_l; j++)
      {
        for (int k=0; k<N_a; k++) 
        { 
          for (int p=0; p<N_R; p++) 
          {
            fftw_destroy_plan(plan[i][j][k][p]);
            plan1[i][j][k][p] = 
              fftw_plan_dft_c2r_1d(N_k_hankel, conv[i][j][k][p], lP[i][j][k][p], FFTW_ESTIMATE);
          }
        }
      }
    }
    #pragma omp parallel for collapse(3)
    for(int i=0; i<N_l; i++)
    {
      for (int p=0; p<N_R; p++) 
      {
        for (int k=0; k<N_a; k++) 
        {
          for(int j=i; j<N_l; j++)
          { 
            for(int q=0; q<N_k_hankel/2+1; q++) 
            { // TODO: ok for dR < 0 in the denom? (original code is neg)
              const double rr = 2.0*M_PI*q/((double) N_k_hankel*(-1.0*dlnk_hankel)); 

              fftw_complex  krn;
              hankel_kernel_FT_3D(rr, &krn, arg, 2);
              
              conv[i][j][k][p][q][0] = flP[i][j][k][p][q][0]*krn[0] - flP[i][j][k][p][q][1]*krn[1];
              conv[i][j][k][p][q][1] = flP[i][j][k][p][q][1]*krn[0] + flP[i][j][k][p][q][0]*krn[1];
            }

            conv[i][j][k][p][0][1] = 0;            
            conv[i][j][k][p][N_k_hankel/2][1] = 0;

            fftw_execute(plan[i][j][k][p]);

            for (int q=0; q<N_k_hankel; q++) 
            {
              const double kk = exp(ln_k_min_hankel + q*dlnk_hankel); 
              
              lP[i][j][k][p][q] = pow(2.0*M_PI*kk, -1.5)*lP[i][j][k][p][q]/((double) N_k_hankel);
              
              // normalization 
              lP[i][j][k][p][q] *= 2.0*2.0*2.0*M_PI*M_PI*M_PI;

              // for interpolation
              lP[i][j][k][p][q] = log(kk*kk*kk*lP[i][j][k][p][q] + 1.0E5);
            }
          }
        }
      }
    }
    for(int i=0; i<N_l; i++)
    {
      for(int j=i; j<N_l; j++)
      {
        for (int k=0; k<N_a; k++) 
        { 
          for (int p=0; p<N_R; p++) 
          {
            fftw_destroy_plan(plan1[i][j][k][p]);
          }
        }
      }
    }
    // ---------------------------------------------------------------------------------------------
    // Step 5 - Make Table
    // ---------------------------------------------------------------------------------------------
    #pragma omp parallel for collapse(4)
    for(int i=0; i<N_l; i++)
    {
      for (int q=0; q<N_k; q++) 
      {
        for (int k=0; k<N_a; k++) 
        {
          for (int p=0; p<N_R; p++) 
          {
            for(int j=i; j<N_l; j++)
            {
              table[i][j][k][q][p] = interpol(lP[i][j][k][p], N_k_hankel, ln_k_min_hankel, 
                ln_k_max_hankel, dlnk_hankel, ln_k_min + q*dlnk, 1.0, 1.0);          
     
              table[j][i][k][q][p] = table[i][j][k][q][p];
            }
          }
        }
      }
    }
    // ---------------------------------------------------------------------------------------------
    // Step 6 - clean up
    // ---------------------------------------------------------------------------------------------
    for(int i=0; i<N_l; i++)
    {
      for(int j=0; j<N_l; j++)
      {
        for (int k=0; k<N_a; k++) 
        { 
          for (int p=0; p<N_R; p++) 
          {
            fftw_free(flP[i][j][k][p]);
            fftw_free(lP[i][j][k][p]); 
            fftw_free(conv[i][j][k][p]);  
          }
          free(plan[i][j][k]);
          free(plan1[i][j][k]);
          free(flP[i][j][k]); 
          free(conv[i][j][k]);  
          free(lP[i][j][k]);
        }
        free(plan[i][j]);
        free(plan1[i][j]);
        free(flP[i][j]);
        free(conv[i][j]);
        free(lP[i][j]);
      }
      free(plan[i]);
      free(plan1[i]);
      free(flP[i]);
      free(conv[i]);
      free(lP[i]);
    }
    free(plan);
    free(plan1);
    free(lP);
    free(flP);
    free(conv);
    // ---------------------------------------------------------------------------------------------
    // Step 7 Update cosmo and nuisance param classes
    // ---------------------------------------------------------------------------------------------   
    update_cosmopara(&C);
    update_nuisance(&N);
  }

  double* ar = (double*) params;
  const int nl1 = ar[0]; 
  const int nl2 = ar[1];
  const int nk = ar[2]; 
  const int na = ar[3];
  const double lnM2 = ar[4];
  if (nl1 < 0 || nl1 > N_l -1 || nl2 < 0 || nl2 > N_l -1)
  {
    log_fatal("invalid bin input (nl1, nl2) = (%d, %d)", nl1, nl2);
    exit(1);
  }
  if (nk < 0 || nk > N_k - 1)
  {
    log_fatal("invalid bin input nk = %d", nk);
    exit(1);
  }
  if (na < 0 || na > N_a - 1)
  {
    log_fatal("invalid bin input na = %d", na);
    exit(1);
  }
  const double aa = amin + na*da;
  const double kk = exp(ln_k_min + nk*dlnk);  
  const double zz = 1.0/aa - 1.0;
  const double M1 = exp(lnM1);  
  const double M2 = exp(lnM2);

  const double B1_1 = B1_x_BSF(M1, aa);
  const double B1_2 = B1_x_BSF(M2, aa); 

  const double tmp = 0.75/(M_PI*Cluster.delta_exclusion*cosmology.rho_crit*cosmology.Omega_m);
  const double R1 = pow(M1*tmp, 1./3.);
  const double R2 = pow(M2*tmp, 1./3.);
  const double R = (0.5*(R1 + R2));
  if (R < R_min || R > R_max)
  {
    log_fatal("R = %e outside look-up table range [%e,%e]", R, R_min, R_max);
    exit(1);
  } 
  if (kk*R > 1)
  {
    return 0.0;
  }
  else
  {
    const double V_EXCL = (4.0*M_PI/3.0)*3.0*(sin(kk*R) - kk*R*cos(kk*R))/(kk*kk*kk);

    const double PK_HALO_EXCL = (exp(interpol(table[nl1][nl2][na][nk], N_R, R_min, R_max, dR, R,
      1.0, 1.0)) - 1.0E5)/(kk*kk*kk);

    const double PK = (R == 0) ?  Pdelta(kk, aa)*B1_1*B1_2 :
      (PK_HALO_EXCL + V_EXCL)*B1_1*B1_2 - V_EXCL;

    double res1;
    {
      double params_in[2] = {nl1, zz};
      res1 = dndlnM_times_binned_P_lambda_obs_given_M(lnM1, (void*) params_in)*BSF(M1);
    }
    double res2;
    {
      double params_in[2] = {nl2, zz};
      res2 = dndlnM_times_binned_P_lambda_obs_given_M(lnM2, (void*) params_in)*BSF(M2);
    }
    return PK*res1*res2;
  }
}

double int_for_binned_p_cc_incl_halo_exclusion(double lnM, void* params)
{
  double* ar = (double*) params;
  const int nl1 = ar[0]; 
  const int nl2 = ar[1];
  const int nk = ar[2];
  const int na = ar[3];
  const int init_static_vars_only = (int) ar[4];

  const double ln_M_min = limits.cluster_util_log_M_min/M_LOG10E;
  const double ln_M_max = limits.cluster_util_log_M_max/M_LOG10E;

  double params_in[5] = {nl1, nl2, nk, na, lnM};
  
  return (init_static_vars_only == 1) ? 
    int_for_int_for_binned_p_cc_incl_halo_exclusion(ln_M_min, (void*) params_in) :
    int_gsl_integrate_low_precision(int_for_int_for_binned_p_cc_incl_halo_exclusion, 
      (void*) params_in, ln_M_min, ln_M_max, NULL, GSL_WORKSPACE_SIZE); 
}

double binned_p_cc_incl_halo_exclusion_nointerp(const int nl1, const int nl2, const int na, 
const int nk, const int init_static_vars_only)
{
  const double N_l = Cluster.N200_Nbin;
  if (nl1 < 0 || nl1 > N_l - 1 || nl2 < 0 || nl2 > N_l - 1)
  {
    log_fatal("invalid bin input (nl1, nl2) = (%d, %d)", nl1, nl2);
    exit(1);
  }  
  const int N_k = Ntable.N_k_halo_exclusion;
  if (nk < 0 || nk > N_k - 1)
  {
    log_fatal("invalid bin input nk = %d", nk);
    exit(1);
  }
  const int N_a = Ntable.N_a_halo_exclusion;
  if (na < 0 || na > N_a - 1)
  {
    log_fatal("invalid bin input na = %d", na);
    exit(1);
  }
  const double amin = 1.0/(1.0 + tomo.cluster_zmax[tomo.cluster_Nbin - 1]);
  const double amax = 1.0/(1.0 + tomo.cluster_zmin[0]) + 0.01;
  const double da = (amax - amin)/((double) N_a - 1.0);
  const double a = amin + na*da;
  const double z = 1.0/a - 1.0;
  const double ln_M_min = limits.cluster_util_log_M_min/M_LOG10E;
  const double ln_M_max = limits.cluster_util_log_M_max/M_LOG10E;

  double norm1;
  {
    double params[2] = {nl1, z};
    norm1 = (init_static_vars_only == 1) ? 
      dndlnM_times_binned_P_lambda_obs_given_M(ln_M_min, (void*) params) :
      int_gsl_integrate_low_precision(dndlnM_times_binned_P_lambda_obs_given_M, (void*) params,
        ln_M_min, ln_M_max, NULL, GSL_WORKSPACE_SIZE); 
  }
  double norm2;
  { 
    double params[2] = {nl2, z};
    norm2 = (init_static_vars_only == 1) ? 
      dndlnM_times_binned_P_lambda_obs_given_M(ln_M_min, (void*) params) :
      int_gsl_integrate_low_precision(dndlnM_times_binned_P_lambda_obs_given_M, (void*) params,
        ln_M_min, ln_M_max, NULL, GSL_WORKSPACE_SIZE);  
  }

  double params[5] = {nl1, nl2, nk, na, init_static_vars_only};
  
  const double res = (init_static_vars_only == 1) ? 
    int_for_binned_p_cc_incl_halo_exclusion(ln_M_min, (void*) params) :
    int_gsl_integrate_low_precision(int_for_binned_p_cc_incl_halo_exclusion, (void*) params,
      ln_M_min, ln_M_max, NULL, GSL_WORKSPACE_SIZE);  
  
  return (norm1*norm2 == 0) ? 0.0 : res/(norm1*norm2);
}

double binned_p_cc_incl_halo_exclusion(const double k, const double a, const int nl1, const int nl2)
{
  static cosmopara C;
  static nuisancepara N;
  static double**** table = 0;

  const int N_l = Cluster.N200_Nbin;

  const int N_k = Ntable.N_k_halo_exclusion;
  const double ln_k_min = log(limits.halo_exclusion_k_min); 
  const double ln_k_max = log(limits.halo_exclusion_k_max);
  const double dlnk = (ln_k_max - ln_k_min)/((double) N_k - 1.0);
  
  const int N_a = Ntable.N_a_halo_exclusion;
  const double amin = 1.0/(1.0 + tomo.cluster_zmax[tomo.cluster_Nbin - 1]);
  const double amax = 1.0/(1.0 + tomo.cluster_zmin[0]) + 0.01;
  const double da = (amax - amin)/((double) N_a - 1.0);

  if (table == 0)
  {
    table = (double****) malloc(sizeof(double***)*N_l);
    for(int i=0; i<N_l; i++)
    {          
      table[i] = (double***) malloc(sizeof(double**)*N_l);
      for(int j=i; j<N_l; j++)
      {
        table[i][j] = (double**) malloc(sizeof(double*)*N_a);
        for (int l=0; l<N_a; l++) 
        { 
          table[i][j][l] = (double*) malloc(sizeof(double)*N_k);
        }
      }
    }
  }
  if (recompute_clusters(C, N))
  { 
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init_static_vars_only = binned_p_cc_incl_halo_exclusion_nointerp(0, 0, 0, 0, 1); 
    }
    #pragma GCC diagnostic pop
    #pragma omp parallel for collapse(3)
    for(int i=0; i<N_l; i++)
    {
      for (int p=0; p<N_k; p++) 
      {
        for (int l=0; l<N_a; l++) 
        { 
          for(int j=i; j<N_l; j++)
          {
            const double aa = amin + l*da;
            const double kk = exp(ln_k_min + p*dlnk);
            table[i][j][l][p] = binned_p_cc_incl_halo_exclusion_nointerp(i, j, l, p, 0);   
            // for interpolation
            table[i][j][l][p] = log(kk*kk*kk*sqrt(aa)*table[i][j][l][p] + 1.E8); 
            table[j][i][l][p] = table[i][j][l][p];
          }
        }
      }
    }
  }
  if (nl1 < 0 || nl1 > N_l - 1 || nl2 < 0 || nl2 > N_l - 1)
  {
    log_fatal("invalid bin input (nl1, nl2) = (%d, %d)", nl1, nl2);
    exit(1);
  }
  if (a < amin || a > amax)
  {
    log_fatal("a = %e outside look-up table range [%e,%e]", a, amin, amax);
    exit(1);
  } 
  const double lnk = log(k);
  if (lnk < ln_k_min || lnk > ln_k_max)
  {
    log_fatal("k = %e outside look-up table range [%e,%e]", k, exp(ln_k_min), exp(ln_k_max));
    exit(1);
  } 
  const double val = interpol2d(table[nl1][nl2], N_a, amin, amax, da, a, N_k, ln_k_min, ln_k_max, 
    dlnk, lnk, 1.0, 1.0);
  return (exp(val) - 1.E8)/(k*k*k*sqrt(a));
}

double binned_p_cc_incl_halo_exclusion_with_constant_lambd(const double k, const double a, 
  const int nl1, const int nl2)
{
  static cosmopara C;
  static nuisancepara N;
  static double**** table = 0;

  const int N_l = Cluster.N200_Nbin;

  const int N_k = Ntable.N_k_halo_exclusion;
  const double ln_k_min = log(limits.halo_exclusion_k_min); 
  const double ln_k_max = log(limits.halo_exclusion_k_max);
  const double dlnk = (ln_k_max - ln_k_min)/((double) N_k - 1.0);
  
  const int N_a = Ntable.N_a_halo_exclusion;
  const double amin = 1.0/(1.0 + tomo.cluster_zmax[tomo.cluster_Nbin - 1]);
  const double amax = 1.0/(1.0 + tomo.cluster_zmin[0]) + 0.01;
  const double da = (amax - amin)/((double) N_a - 1.0);

  const int N_k_hankel = Ntable.N_k_hankel_halo_exclusion;
  const double ln_k_min_hankel = log(limits.halo_exclusion_k_min_hankel);
  const double ln_k_max_hankel = log(limits.halo_exclusion_k_max_hankel);
  const double dlnk_hankel = (ln_k_max_hankel - ln_k_min_hankel)/((double) N_k_hankel - 1.0);

  typedef fftw_complex fftwZ;
  const double CO = 1.0; // cutoff = CO
  double arg[2];
  arg[0] = 0;     // bias
  arg[1] = 0.5;   // order of Bessel function  

  if (table == 0)
  {
    table = (double****) malloc(sizeof(double***)*N_l);
    for(int i=0; i<N_l; i++)
    {          
      table[i] = (double***) malloc(sizeof(double**)*N_l);
      for(int j=0; j<N_l; j++)
      {
        table[i][j] = (double**) malloc(sizeof(double*)*N_a);
        for (int l=0; l<N_a; l++) 
        { 
          table[i][j][l] = (double*) malloc(sizeof(double)*N_k);
        }
      }
    }
  }
  if (recompute_clusters(C, N))
  { 
    // ---------------------------------------------------------------------------------------------
    // Step 1 - memory allocation (except for FFTW plans)
    // --------------------------------------------------------------------------------------------- 
    fftw_plan*** plan   = (fftw_plan***) malloc(sizeof(fftw_plan**)*N_l);
    fftw_plan*** plan1  = (fftw_plan***) malloc(sizeof(fftw_plan**)*N_l);
    fftwZ**** flP      = (fftwZ****)    malloc(sizeof(fftwZ***)*N_l);
    fftwZ**** conv      = (fftwZ****)    malloc(sizeof(fftwZ***)*N_l);
    double**** lP       = (double****)   malloc(sizeof(double***)*N_l);
    for(int i=0; i<N_l; i++)
    {
      plan[i]  = (fftw_plan**) malloc(sizeof(fftw_plan*)*N_l);
      plan1[i] = (fftw_plan**) malloc(sizeof(fftw_plan*)*N_l);    
      flP[i]  = (fftwZ***)    malloc(sizeof(fftwZ**)*N_l);
      conv[i]  = (fftwZ***)    malloc(sizeof(fftwZ**)*N_l);
      lP[i]    = (double***)   malloc(sizeof(double**)*N_l);
      for(int j=0; j<N_l; j++)
      {
        plan[i][j]  = (fftw_plan*) malloc(sizeof(fftw_plan)*N_a);
        plan1[i][j] = (fftw_plan*) malloc(sizeof(fftw_plan)*N_a);
        flP[i][j]  = (fftwZ**)    malloc(sizeof(fftwZ*)*N_a);
        conv[i][j]  = (fftwZ**)    malloc(sizeof(fftwZ*)*N_a);
        lP[i][j]    = (double**)   malloc(sizeof(double*)*N_a);
        for (int l=0; l<N_a; l++) 
        { 
          flP[i][j][l] = fftw_malloc((N_k_hankel/2 + 1)*sizeof(fftw_complex));
          conv[i][j][l] = fftw_malloc((N_k_hankel/2 + 1)*sizeof(fftw_complex));
          lP[i][j][l]   = fftw_malloc(N_k_hankel*sizeof(double));
        }
      }
    }
    // ---------------------------------------------------------------------------------------------
    // Step 2: - pk_to_xi
    // ---------------------------------------------------------------------------------------------
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init_static_vars_only  = weighted_B1(0, 1.0/amin - 1.0);     
      double init_static_vars_only2 = Pdelta(exp(ln_k_min_hankel), amin); 
    }
    #pragma GCC diagnostic pop
    #pragma omp parallel for collapse(3)
    for(int i=0; i<N_l; i++)
    {
      for (int l=0; l<N_a; l++) 
      { 
        for (int p=0; p<N_k_hankel; p++) 
        {
          const double aa = amin + l*da;
          const double zz = 1.0/aa - 1.0;
          const double kk = exp(ln_k_min_hankel + p*dlnk_hankel);
          const double B1_1 = weighted_B1(i, zz);
          for(int j=i; j<N_l; j++)
          {
            const double B1_2 = (i == j) ? B1_1 : weighted_B1(j, zz);
            lP[i][j][l][p] = kk*sqrt(kk)*Pdelta(kk, aa)*B1_1*B1_2;
          }
        }
      }
    }
    for(int i=0; i<N_l; i++)
    {
      for(int j=i; j<N_l; j++)
      {
        for (int l=0; l<N_a; l++) 
        { 
          plan[i][j][l] = 
            fftw_plan_dft_r2c_1d(N_k_hankel, lP[i][j][l], flP[i][j][l], FFTW_ESTIMATE);
        }
      }
    }
    #pragma omp parallel for collapse(2)
    for(int i=0; i<N_l; i++)
    {
      for (int l=0; l<N_a; l++) 
      {
        for(int j=i; j<N_l; j++)
        {
          fftw_execute(plan[i][j][l]);
        }
      }
    }
    for(int i=0; i<N_l; i++)
    {
      for(int j=i; j<N_l; j++)
      {
        for (int l=0; l<N_a; l++) 
        { 
          fftw_destroy_plan(plan[i][j][l]);
          plan1[i][j][l] = 
            fftw_plan_dft_c2r_1d(N_k_hankel, conv[i][j][l], lP[i][j][l], FFTW_ESTIMATE);
        }
      }
    }
    #pragma omp parallel for collapse(2)
    for(int i=0; i<N_l; i++)
    {
      for (int l=0; l<N_a; l++) 
      {
        for(int j=i; j<N_l; j++)
        { 
          for(int p=0; p<N_k_hankel/2+1; p++) 
          {            
            fftw_complex  kernel;
            
            hankel_kernel_FT_3D(2.0*M_PI*p/((double) N_k_hankel*dlnk_hankel), &kernel, arg, 2);
            
            conv[i][j][l][p][0] = flP[i][j][l][p][0]*kernel[0] - flP[i][j][l][p][1]*kernel[1];
            conv[i][j][l][p][1] = flP[i][j][l][p][1]*kernel[0] + flP[i][j][l][p][0]*kernel[1];
          }

          conv[i][j][l][0][1] = 0;            
          conv[i][j][l][N_k_hankel/2][1] = 0;

          fftw_execute(plan1[i][j][l]);

          for (int p=0; p<N_k_hankel; p++) 
          {
            const double kk = exp(ln_k_min_hankel + p*dlnk_hankel);    
            const double rr = 1.0/kk;
            lP[i][j][l][p] = pow(2.0*M_PI*rr, -1.5)*lP[i][j][l][p]/((double) N_k_hankel); // XI
          }
        }
      }
    }
    for(int i=0; i<N_l; i++)
    {
      for(int j=i; j<N_l; j++)
      {
        for (int l=0; l<N_a; l++) 
        { 
          fftw_destroy_plan(plan1[i][j][l]);
        }
      }
    }
    // ---------------------------------------------------------------------------------------------
    // STEP 3: exclusion_filter
    // ---------------------------------------------------------------------------------------------
    #pragma omp parallel for collapse(3)
    for(int i=0; i<N_l; i++)
    {
      for (int p=0; p<N_k_hankel; p++) 
      {
        for (int l=0; l<N_a; l++) 
        {
          for(int j=i; j<N_l; j++)
          {
            const double aa = amin + l*da;
            const double kk = exp(ln_k_min_hankel + p*dlnk_hankel);     
            const double rr = 1.0/kk;
            const double R = 1.5*pow(0.25*(Cluster.N_min[i] + Cluster.N_min[j] + 
              Cluster.N_max[i] + Cluster.N_max[j])/100., 0.2)/cosmology.coverH0/aa;
            
            if (Cluster.halo_exclusion_model == 0)
            {
              if (rr < R)
              {
                lP[i][j][l][p] = -1; 
              } 
            }
            else if (Cluster.halo_exclusion_model == 1)
            { // Baldauf 2013 eq. 43
              if (R > 0) 
              {
                gsl_sf_result tmp;
                int status = gsl_sf_erf_e(log(rr/R)/(M_SQRT2*0.0387), &tmp);
                if (status)
                {
                  log_fatal(gsl_strerror(status));
                  exit(1);
                }
                lP[i][j][l][p] = 0.5*(1.0 + tmp.val)*(lP[i][j][l][p] + 1.0) - 1.0;
              }
            }
            else
            {  
              log_fatal("halo exclusion model: %d  not implemented", Cluster.halo_exclusion_model);
              exit(1);
            }
          }
        }
      }
    }
    // ---------------------------------------------------------------------------------------------
    // STEP 4: - xi_to_pk
    // ---------------------------------------------------------------------------------------------
    #pragma omp parallel for collapse(3)
    for(int i=0; i<N_l; i++)
    {
      for (int p=0; p<N_k_hankel; p++) 
      {    
        for (int l=0; l<N_a; l++) 
        { 
          const double kk = exp(ln_k_min_hankel + p*dlnk_hankel);    
          const double rr = 1.0/kk;
          for(int j=i; j<N_l; j++)
          {
            lP[i][j][l][p] = rr*sqrt(rr)*lP[i][j][l][p];
          }
        }
      }
    }
    for(int i=0; i<N_l; i++)
    {
      for(int j=i; j<N_l; j++)
      {
        for (int l=0; l<N_a; l++) 
        { 
          plan[i][j][l] = fftw_plan_dft_r2c_1d(N_k_hankel, lP[i][j][l], flP[i][j][l], FFTW_ESTIMATE);
        }
      }
    }
    #pragma omp parallel for collapse(2)
    for(int i=0; i<N_l; i++)
    {
      for (int l=0; l<N_a; l++) 
      { 
        for(int j=i; j<N_l; j++)
        { 
          fftw_execute(plan[i][j][l]);  
        }
      }
    }
    for(int i=0; i<N_l; i++)
    {
      for(int j=i; j<N_l; j++)
      {
        for (int l=0; l<N_a; l++) 
        { 
          fftw_destroy_plan(plan[i][j][l]);
          plan1[i][j][l] = fftw_plan_dft_c2r_1d(N_k_hankel, conv[i][j][l], lP[i][j][l], FFTW_ESTIMATE);
        }
      }
    }
    #pragma omp parallel for collapse(2)
    for(int i=0; i<N_l; i++)
    {
      for (int l=0; l<N_a; l++) 
      { 
        for(int j=i; j<N_l; j++)
        {           
          for(int p=0; p<N_k_hankel/2+1; p++) 
          { // TODO: ok for dR < 0 in the denom? (original code is neg)
            const double rr = 2.0*M_PI*p/((double) N_k_hankel*(-1.0*dlnk_hankel)); 

            fftw_complex  kernel;
            hankel_kernel_FT_3D(rr, &kernel, arg, 2);
            
            conv[i][j][l][p][0] = flP[i][j][l][p][0]*kernel[0] - flP[i][j][l][p][1]*kernel[1];
            conv[i][j][l][p][1] = flP[i][j][l][p][1]*kernel[0] + flP[i][j][l][p][0]*kernel[1];
          }

          conv[i][j][l][0][1] = 0;            
          conv[i][j][l][N_k_hankel/2][1] = 0;

          fftw_execute(plan1[i][j][l]);

          for (int p=0; p<N_k_hankel; p++) 
          {
            const double kk = exp(ln_k_min_hankel + p*dlnk_hankel); 
            lP[i][j][l][p] = pow(2*M_PI*kk, -1.5)*lP[i][j][l][p]/((double) N_k_hankel);
            
            // normalization 
            lP[i][j][l][p] *= 2.0*2.0*2.0*M_PI*M_PI*M_PI;

            // for interpolation
            lP[i][j][l][p] = log(kk*kk*kk*lP[i][j][l][p] + 1.0E5); // 1E5 avoid negs inside ln()
          }
        }
      }
    }
    for(int i=0; i<N_l; i++)
    {
      for(int j=i; j<N_l; j++)
      {
        for (int l=0; l<N_a; l++) 
        { 
          fftw_destroy_plan(plan1[i][j][l]);
        }
      }
    }
    // ---------------------------------------------------------------------------------------------
    // Step 5 - Evaluate binned_p_cc_incl_halo_exclusion in the interpolation table
    // ---------------------------------------------------------------------------------------------
    #pragma omp parallel for collapse(3)
    for(int i=0; i<N_l; i++)
    {
      for (int p=0; p<N_k; p++) 
      {
        for (int l=0; l<N_a; l++) 
        {
          for(int j=i; j<N_l; j++)
          {
            const double aa = amin + l*da;
            const double zz = 1.0/aa - 1.0;
            const double kk = exp(ln_k_min + p*dlnk);

            const double B1_1 = weighted_B1(i, zz); 
            const double B1_2 = (i == j) ? B1_1 : weighted_B1(j, zz);  
            const double B1B2 = B1_1*B1_2;
            const double PKB1B2 = Pdelta(kk, aa)*B1B2;

            // RedMaPPer exclusion radius (comoving) (Rykoff et al. 2014  eq4) in units coverH0
            const double R = 1.5*pow(0.25*(Cluster.N_min[i] + Cluster.N_min[j] + 
                Cluster.N_max[i] + Cluster.N_max[j])/100.0, 0.2)/cosmology.coverH0/aa;

            const double xx = kk*R;

            const double V_EXCL = 4.0*M_PI*(sin(xx) - xx*cos(xx))/(kk*kk*kk);

            if (xx > CO)
            {
              const double kk_CO = CO/R;
              
              const double ln_kk_CO = log(kk_CO);
              if (ln_kk_CO < ln_k_min || ln_kk_CO > ln_k_max)
              {
                log_fatal(
                "k = %e outside look-up table range [%e,%e]", kk_CO, exp(ln_k_min), exp(ln_k_max));
                exit(1);
              } 

              const double PKB1B2_CO =  Pdelta(kk_CO, aa)*B1B2;
              
              const double V_EXCL_CO =  4.0*M_PI*(sin(CO) - CO*cos(CO))/(kk_CO*kk_CO*kk_CO);

              const double tmp = interpol(lP[i][j][l], N_k_hankel, ln_k_min_hankel, 
                ln_k_max_hankel, dlnk_hankel, log(kk_CO), 1.0, 1.0);

              const double PK_HALO_EXCL_CO = (exp(tmp) - 1.0E5)/(kk_CO*kk_CO*kk_CO);

              const double P_EXCL_CO = (PK_HALO_EXCL_CO + V_EXCL_CO)*B1B2 - V_EXCL_CO;

              const double ENV = pow(kk/kk_CO, -0.7);

              table[i][j][l][p] = 
                (PKB1B2 - V_EXCL*(1 + (PKB1B2_CO - V_EXCL_CO - P_EXCL_CO)/V_EXCL_CO))*ENV;
            }
            else if (R == 0.0)
            {
              table[i][j][l][p] = PKB1B2;
            }
            else
            {
              const double tmp = interpol(lP[i][j][l], N_k_hankel, ln_k_min_hankel, ln_k_max_hankel,
                dlnk_hankel, log(kk), 1.0, 1.0);
              
              const double PK_HALO_EXCL = (exp(tmp) - 1.0E5)/(kk*kk*kk);

              table[i][j][l][p] = (PK_HALO_EXCL + V_EXCL)*B1B2 - V_EXCL;
            }
            
            // for interpolation
            table[i][j][l][p] =  log(kk*kk*kk*sqrt(aa)*table[i][j][l][p] + 1.E8);
            table[j][i][l][p] =  table[i][j][l][p];
          }
        }
      }
    }
    // ---------------------------------------------------------------------------------------------
    // Step 6 - clean up
    // ---------------------------------------------------------------------------------------------
    for(int i=0; i<N_l; i++)
    {
      for(int j=0; j<N_l; j++)
      {
        for (int l=0; l<N_a; l++) 
        { 
          fftw_free(flP[i][j][l]); 
          fftw_free(conv[i][j][l]);  
          fftw_free(lP[i][j][l]);   
        }
        free(plan[i][j]);
        free(plan1[i][j]);
        free(flP[i][j]);
        free(conv[i][j]);
        free(lP[i][j]);
      }
      free(plan[i]);
      free(plan1[i]);
      free(flP[i]);
      free(conv[i]);
      free(lP[i]);
    }
    free(plan);
    free(plan1);
    free(flP);
    free(conv);
    free(lP);
    // ---------------------------------------------------------------------------------------------
    // Step 7 Update cosmo and nuisance param classes
    // ---------------------------------------------------------------------------------------------    
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  if (nl1 < 0 || nl1 > N_l - 1 || nl2 < 0 || nl2 > N_l - 1)
  {
    log_fatal("invalid bin input (nl1, nl2) = (%d, %d)", nl1, nl2);
    exit(1);
  }
  if (a < amin || a > amax)
  {
    log_fatal("a = %e outside look-up table range [%e,%e]", a, amin, amax);
    exit(1);
  } 
  const double lnk = log(k);
  if (lnk < ln_k_min || lnk > ln_k_max)
  {
    log_fatal("k = %e outside look-up table range [%e,%e]", k, exp(ln_k_min), exp(ln_k_max));
    exit(1);
  } 
  const double val = interpol2d(table[nl1][nl2], N_a, amin, amax, da, a, N_k, ln_k_min, ln_k_max, 
    dlnk, lnk, 1.0, 1.0);
  return (exp(val) - 1.E8)/(k*k*k*sqrt(a));
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// CLUSTER CROSS SPECTRUM WITH GALAXIES 
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

// nl = lambda_obs bin, nj = galaxy redshift bin
double binned_p_cg(const double k, const double a, const int nl, const int nj, 
  const int use_linear_ps)
{ // binend in lambda_obs
  if (!(a>0) || !(a<1))
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  const double z = 1.0/a - 1.0;
  const double cluster_bias = weighted_B1(nl, z);
  const double PK = use_linear_ps == 1 ? p_lin(k, a) : Pdelta(k, a);
  const double PCG = cluster_bias * gbias.b1_function(z, nj) * PK;

  double PCG_1L = 0.;
  if (has_b2_galaxies() && Cluster.nonlinear_bias == 1 && use_linear_ps != 1)
  {
    const double g = growfac(a);
    const double g4 = g*g*g*g;
    
    const double b1c = weighted_B1(nl, z);
    const double b2c = weighted_B2(nl, z);
    const double bs2c = bs2_from_b1(b1c);

    const double b1g = gbias.b1_function(z, nj);
    const double b2g = (gbias.b2[nj]  == 0) ? b2_from_b1(gbias.b[nj]) : gbias.b2[nj];
    const double bs2g = gbias.bs2[nj];
    
    PCG_1L = 0.5*(b1c*b2g+b2c*b1g)*PT_d1d2(k) + 0.25*b2g*b2c*PT_d2d2(k) 
      + 0.5*(b1c*bs2g+b1g*bs2c)*PT_d1s2(k) + 0.25*(b2c*bs2g+b2g*bs2c)*PT_d2s2(k) 
      + 0.25*(bs2g*bs2c)*PT_s2s2(k) + 0.5*(b1c*b3nl_from_b1(b1g) 
      + b1g*b3nl_from_b1(b1c))*PT_d1d3(k);
    PCG_1L *= g4;
  }
  // Cosmolike: this makes the code much much faster like 1000 times faster
  return (PCG + PCG_1L) < 0 ? 0.0 : (PCG + PCG_1L);
}

// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
// cluster number counts
// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------
// nl = lambda_obs bin, ni = cluster redshift bin

double binned_Ndensity_nointerp(const int nl, const double z, const int init_static_vars_only)
{
  const double ln_M_min = limits.cluster_util_log_M_min/M_LOG10E;
  const double ln_M_max = limits.cluster_util_log_M_max/M_LOG10E;
  double params[2] = {(double) nl, z};
  return (init_static_vars_only == 1) ?
    dndlnM_times_binned_P_lambda_obs_given_M(ln_M_min, (void*) params) :
    int_gsl_integrate_low_precision(dndlnM_times_binned_P_lambda_obs_given_M, (void*) params,
      ln_M_min, ln_M_max, NULL, GSL_WORKSPACE_SIZE);
}

double binned_Ndensity(const int nl, const double z)
{
  static cosmopara C;
  static nuisancepara N;
  static double** table;

  const int N_l = Cluster.N200_Nbin;
  const int N_a = Ntable.N_a;
  const double zmin = fmax(tomo.cluster_zmin[0] - 0.05, 0.01);
  const double zmax = tomo.cluster_zmax[tomo.cluster_Nbin - 1] + 0.05;
  const double amin = 1.0/(1.0 + zmax);
  const double amax = 1.0/(1.0 + zmin);
  const double da = (amax - amin)/((double) N_a - 1.0);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*N_l);
    for (int i=0; i<N_l; i++)
    {
      table[i] = (double*) malloc(sizeof(double)*N_a);
    }
  }
  if (recompute_clusters(C, N))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init_static_vars_only = binned_Ndensity_nointerp(0, 1.0/amin - 1.0, 1);
    }
    #pragma GCC diagnostic pop
    #pragma omp parallel for collapse(2)
    for (int i=0; i<N_l; i++)
    {
      for (int j=0; j<N_a; j++)
      {
        const double aa = amin + j*da;
        table[i][j] = binned_Ndensity_nointerp(i, 1.0/aa - 1.0, 0);
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
  if (z < zmin || z > zmax)
  {
    log_fatal("z = %e outside look-up table range [%e,%e]", z, zmin, zmax);
    exit(1);
  }
  return interpol(table[nl], N_a, amin, amax, da, 1.0/(z + 1.0), 0., 0.);
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// CLUSTER CROSS SPECTRUM WITH MATTER)
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

double int_for_binned_p_cm(double lnM, void* params)
{ // binned in lambda_obs (nl = lambda_obs bin)
  double* ar = (double*) params;
  const int nl = (int) ar[0];
  const double z = ar[1];
  const double k = ar[2];
  const int include_1h_term = (int) ar[3];
  const int use_linear_ps = (int) ar[4];
  
  const double M = exp(lnM); 
  const double a = 1.0/(1.0 + z);
  if (!(a>0) || !(a<1))
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }

  const double PCM_1H = (include_1h_term == 1) ? 
    M/(cosmology.rho_crit*cosmology.Omega_m)*u_nfw_c(conc(M, a), k, M, a) : 0.0; 
   
  const double B1BS = B1_x_BSF(M, a);
  const double PCM_2H = use_linear_ps == 1 ? B1BS*p_lin(k, a) : B1BS*Pdelta(k, a);
  
  double params_in[2] = {(double) nl, z};
  return (PCM_1H + PCM_2H)*dndlnM_times_binned_P_lambda_obs_given_M(lnM, (void*) params_in);
}


// nl = lambda_obs bin
double binned_p_cm_nointerp(const double k, const double a, const int nl, const int include_1h_term,
const int use_linear_ps, const int init_static_vars_only)
{
  if (!(a>0) || !(a<1))
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  const double z = 1.0/a - 1.0;
  
  double params[5] = {(double) nl, z, k, (double) include_1h_term, (double) use_linear_ps};
  
  const double ln_M_min = limits.cluster_util_log_M_min/M_LOG10E;
  const double ln_M_max = limits.cluster_util_log_M_max/M_LOG10E;
  
  const double norm = binned_Ndensity_nointerp(nl, z, init_static_vars_only);

  return (init_static_vars_only == 1) ? int_for_binned_p_cm(ln_M_min, (void*) params) :
    int_gsl_integrate_low_precision(int_for_binned_p_cm, (void*) params, ln_M_min, ln_M_max, NULL, 
      GSL_WORKSPACE_SIZE)/norm; 
}

// nl = lambda_obs bin, ni = cluster redshift bin
double binned_p_cm(const double k, const double a, const int nl, const int ni, 
const int include_1h_term, const int use_linear_ps)
{
  static cosmopara C;
  static nuisancepara N;
  static double**** table = 0;
  static double* amin = 0;
  static double* amax = 0;
  static double* da = 0;

  const int N_l = Cluster.N200_Nbin;
  const int N_z = tomo.cluster_Nbin;
  const int N_k = Ntable.N_ell;
  const double ln_k_min = log(limits.k_min_cH0);
  const double ln_k_max = log(limits.k_max_cH0);
  const double dlnk = (ln_k_max - ln_k_min)/((double) N_k - 1.0);
  const int N_a = Ntable.binned_p_cm_size_a_table;

  if (table == 0)
  {
    table = (double****) malloc(sizeof(double***)*N_l);
    for (int i=0; i<N_l; i++)
    {
      table[i] = (double***) malloc(sizeof(double**)*N_z);
      for (int j=0; j<N_z; j++)
      {
        table[i][j] = (double**) malloc(sizeof(double*)*N_k);
        for (int l=0; l<N_k; l++)
        {
          table[i][j][l] = (double*) malloc(sizeof(double)*N_a);
        }
      }
    }
    amin = (double*) malloc(sizeof(double)*N_z);
    amax = (double*) malloc(sizeof(double)*N_z);
    da = (double*) malloc(sizeof(double)*N_z);
    for(int j=0; j<N_z; j++)
    {
      const double zmin = tomo.cluster_zmin[j];
      const double zmax = tomo.cluster_zmax[j];
      amin[j] = 1.0/(1.0 + zmax); 
      amax[j] = 1.0/(1.0 + zmin);
      da[j] = (amax - amin)/((double) N_a - 1.0);
    }
  }
  if (recompute_clusters(C, N))
  { 
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {    
      double init_static_vars_only = binned_p_cm_nointerp(exp(ln_k_min), amin[0], 0, 
        include_1h_term, use_linear_ps, 1); 
    }
    #pragma GCC diagnostic pop
    #pragma omp parallel for collapse(4)
    for(int i=0; i<N_l; i++)
    {
      for(int j=0; j<N_z; j++)
      {
        for(int l=0; l<N_k; l++)
        {
          for(int p=0; p<N_a; p++)
          { 
            table[i][j][l][p] = binned_p_cm_nointerp(exp(ln_k_min + l*dlnk), amin[j] + p*da[j], i,
              include_1h_term, use_linear_ps, 0);
          }
        }
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);   
  }
  if (ni < 0 || ni > N_z - 1 || nl < 0 || nl > N_l - 1)
  {
    log_fatal("invalid bin input (nl, ni) = (%d, %d)", nl, ni);
    exit(1);
  }
  if (a < amin[ni] || a > amax[ni])
  {
    log_fatal("a = %e outside look-up table range [%e,%e]", a, amin[ni], amax[ni]);
    exit(1);
  }
  const double lnk = log(k);
  if (lnk < ln_k_min || lnk > ln_k_max)
  {
    log_fatal("k = %e outside look-up table range [%e,%e]", k, exp(ln_k_min), exp(ln_k_max));
    exit(1);
  } 
  return interpol2d(table[nl][ni], N_k, ln_k_min, ln_k_max, dlnk, lnk, N_a, amin[ni], amax[ni], 
    da[ni], a, 1.0, 1.0); 
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// Area
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

// Cocoa: we try to avoid reading of files in the cosmolike_core code 
// Cocoa: (reading is done in the C++/python interface)
void setup_get_area(int* io_nz, double** io_z, double** io_A, int io) 
{
  static int nz;
  static double* z = NULL;
  static double* A = NULL;

  if (io == 1)
  { // IO == 1 IMPLES THAT IO_A(Z) WILL COPIED TO LOCAL A(Z)
    if (z != NULL)
    {
      free(z);
    }
    if (A != NULL)
    {
      free(A);
    }
    nz = (*io_nz);
    if (!(nz > 5))
    {
      log_fatal("array to small");
      exit(1);
    }
    z = (double*) malloc(nz*sizeof(double));
    A = (double*) malloc(nz*sizeof(double));
    if (z == NULL || A == NULL)
    {
      log_fatal("fail allocation");
      exit(1);
    }
    for (int i=0; i<nz; i++)
    {
      z[i] = (*io_z)[i];
      A[i] = (*io_A)[i];
    }
  }
  else
  { // IO != 1 IMPLES THAT LOCAL A(Z) WILL BE COPIED TO IO_A(Z)
    if (z == NULL || A == NULL)
    {
      log_fatal("array/pointer not allocated");
      exit(1);
    }
    if ((*io_z) != NULL)
    {
      free((*io_z));
      (*io_z) = NULL;
    }
    else if ((*io_A) != NULL)
    {
      free((*io_A));
      (*io_A) = NULL;
    }
    (*io_z) = (double*) malloc(nz*sizeof(double));
    (*io_A) = (double*) malloc(nz*sizeof(double));
    if ((*io_z) == NULL || (*io_A) == NULL)
    {
      log_fatal("fail allocation");
      exit(1);
    }
    for (int i=0; i<nz; i++)
    {
      (*io_z)[i] = z[i];
      (*io_A)[i] = A[i];
    }
  }
}

double get_area(const double zz, const int interpolate_survey_area)
{
  static gsl_spline* fA = NULL;

  if (interpolate_survey_area == 1)
  {
    if (fA == NULL)
    {
      int nz;
      double** z = (double**) malloc(1*sizeof(double*));
      double** A = (double**) malloc(1*sizeof(double*));
         
      if (z == NULL || A == NULL)
      {
        log_fatal("fail allocation");
        exit(1);
      }
      else 
      {
        (*z) = NULL;
        (*A) = NULL;
      }

      setup_get_area(&nz, z, A, 0);

      const gsl_interp_type* T = gsl_interp_linear;
      fA = gsl_spline_alloc(T, nz);
      if (fA == NULL)
      {
        log_fatal("fail allocation");
        exit(1);
      }

      int status = 0;
      status = gsl_spline_init(fA, (*z), (*A), nz);
      if (status) 
      {
        log_fatal(gsl_strerror(status));
        exit(1);
      }
      free(z[0]); // spline makes a copy of the data
      free(A[0]); // spline makes a copy of the data
      free(z);    // spline makes a copy of the data
      free(A);    // spline makes a copy of the data
    }
    double res;
    int status = gsl_spline_eval_e(fA, zz, NULL, &res);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
    return res;
  }
  else
  {
    return survey.area;
  }
}