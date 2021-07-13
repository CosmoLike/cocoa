#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>

#include "basics.h"
#include "bias.h"
#include "cosmo3D.h"
#include "cluster_util.h"
#include "halo.h"
#include "recompute.h"
#include "radial_weights.h"
#include "redshift_spline.h"
#include "structs.h"

#include "log.c/src/log.h"

static int GSL_WORKSPACE_SIZE = 250;
static double log_M_min = 12.0;
static double log_M_max = 15.9; 

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// BUZZARD
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

double buzzard_P_lambda_obs_given_M_nointerp(double obs_lambda, void* params)
{      
  double* array = (double*) params;
  
  const double mass = array[0];
  const double cluster_z = array[1];
  const double lnlm0 = nuisance.cluster_MOR[0]; 
  const double Alm = nuisance.cluster_MOR[1];
  const double sigma_lnlm_intrinsic = nuisance.cluster_MOR[2];
  const double Blm = nuisance.cluster_MOR[3];
  const double Mpiv = 5E14; //Msun/h
  const double lnlm = lnlm0 + (Alm)*log(mass/Mpiv) + Blm*log((1+cluster_z)/1.45);
  const double sigma_total = (lnlm>0) ? 
    sqrt(sigma_lnlm_intrinsic*sigma_lnlm_intrinsic+(exp(lnlm)-1.)/exp(2*lnlm)) : 
    sqrt(sigma_lnlm_intrinsic*sigma_lnlm_intrinsic);
  const double x = 1.0/2.0*(log(obs_lambda)-lnlm)*(log(obs_lambda)-lnlm)/pow(sigma_total,2.0);
  
  return exp(-x)/M_SQRTPI/M_SQRT2/sigma_total/obs_lambda;
}

// \int_(bin_lambda_obs_min)^(bin_lambda_obs_max) \dlambda_obs P(\lambda_obs|M)
// (see for example https://arxiv.org/pdf/1810.09456.pdf - eq3 qnd 6) 
double buzzard_binned_P_lambda_obs_given_M_nointerp(int nl, double M, double z)
{
  double params[2] = {M, z};

  const double bin_lambda_obs_min = Cluster.N_min[nl];
  const double bin_lambda_obs_max = Cluster.N_max[nl];

  return int_gsl_integrate_medium_precision(buzzard_P_lambda_obs_given_M_nointerp,
    (void*) params, bin_lambda_obs_min, bin_lambda_obs_max, NULL, 1000);
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// SDSS
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

// Cocoa: we try to avoid reading of files in the cosmolike_core code 
// Cocoa: (reading is done in the C++/python interface)
void setup_SDSS_P_true_lambda_given_mass(int* io_nintrinsic_sigma,  double** io_intrinsic_sigma, 
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
    if(intrinsic_sigma == NULL || atsrgm == NULL || alpha == NULL || sigma == NULL) 
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
    if((*io_intrinsic_sigma) != NULL)
    {
      free((*io_intrinsic_sigma));
      (*io_intrinsic_sigma) = NULL;
    }
    else if((*io_atsrgm) != NULL)
    {
      free((*io_atsrgm));
      (*io_atsrgm) = NULL;
    }
    else if((*io_alpha) != NULL)
    {
      free((*io_alpha));
      (*io_alpha) = NULL;
    }
    else if((*io_sigma) != NULL)
    {
      free((*io_sigma));
      (*io_sigma) = NULL;
    }

    (*io_intrinsic_sigma) = (double*) malloc(nintrinsic_sigma*sizeof(double));
    (*io_atsrgm) = (double*) malloc(natsrgm*sizeof(double));
    (*io_alpha) = (double*) malloc(nintrinsic_sigma*natsrgm*sizeof(double));
    (*io_sigma) = (double*) malloc(nintrinsic_sigma*natsrgm*sizeof(double));
    if(*io_intrinsic_sigma == NULL || *io_atsrgm == NULL || *io_alpha == NULL || *io_sigma == NULL) 
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

// SKEW-NORMAL APPROXIMATION eq B1 of https://arxiv.org/pdf/1810.09456.pdf
double SDSS_P_true_lambda_given_mass(const double true_lambda, const double mass, const double z)
{
  static int first = 0;
  static gsl_spline2d* falpha = NULL; // skewness of the skew-normal distribution
  static gsl_spline2d* fsigma = NULL; // variance of the skew-normal distribution

  if(first == 0) 
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

    // we don't want to guess the appropriate GSL z array ordering in z = f(x, y)
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
  // intrisic scatter of mass-richness relation
  const double intrinsic_sigma = nuisance.cluster_MOR[3];
  
  // average true satellite richness given mass
  const double tmp = (mass - mass_min)/(mass_M1 - mass_min);
  double atsrgm = pow(tmp, intrinsic_alpha)*pow(((1+z)/1.45), nuisance.cluster_MOR[4]); 
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
  double sigma = 0.0;
  status = gsl_spline2d_eval_e(falpha, intrinsic_sigma, atsrgm, NULL, NULL, &alpha);
  if (status) 
  {
    log_fatal(gsl_strerror(status));
    exit(1);
  }
  status = gsl_spline2d_eval_e(fsigma, intrinsic_sigma, atsrgm, NULL, NULL, &sigma);
  if (status) 
  {
    log_fatal(gsl_strerror(status));
    exit(1);
  }

  const double x = 1.0/(M_SQRT2*abs(sigma));
  const double result = exp(-(true_lambda-atsrgm)*(true_lambda-atsrgm)*x*x)*x/M_SQRTPI;
  return result*gsl_sf_erfc(-alpha*(true_lambda-atsrgm)*x);
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

    z = (double*) malloc(nz*sizeof(double));
    lambda = (double*) malloc(nlambda*sizeof(double));
    tau = (double*) malloc(nz*nlambda*sizeof(double));
    mu = (double*) malloc(nz*nlambda*sizeof(double));
    sigma = (double*) malloc(nz*nlambda*sizeof(double));
    fmask = (double*) malloc(nz*nlambda*sizeof(double));
    fprj = (double*) malloc(nz*nlambda*sizeof(double));
    if(z == NULL || lambda == NULL || tau == NULL || mu == NULL || sigma == NULL || 
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
    if((*io_z) != NULL)
    {
      free((*io_z));
      (*io_z) = NULL;
    }
    else if((*io_lambda) != NULL)
    {
      free((*io_lambda));
      (*io_lambda) = NULL;
    }
    else if((*io_tau) != NULL)
    {
      free((*io_tau));
      (*io_tau) = NULL;
    }
    else if((*io_mu) != NULL)
    {
      free((*io_mu));
      (*io_mu) = NULL;
    }
    else if((*io_sigma) != NULL)
    {
      free((*io_sigma));
      (*io_sigma) = NULL;
    }
    else if((*io_fmask) != NULL)
    {
      free((*io_fmask));
      (*io_fmask) = NULL;
    }   
    else if((*io_fprj) != NULL)
    {
      free((*io_fprj));
      (*io_fprj) = NULL;
    }
   
    (*io_z) = (double*) malloc(nz*sizeof(double));
    (*io_lambda) = (double*) malloc(nlambda*sizeof(double));
    (*io_tau) = (double*) malloc(nz*nlambda*sizeof(double));
    (*io_mu) = (double*) malloc(nz*nlambda*sizeof(double));
    (*io_sigma) = (double*) malloc(nz*nlambda*sizeof(double));
    (*io_fmask) = (double*) malloc(nz*nlambda*sizeof(double));
    (*io_fprj) = (double*) malloc(nz*nlambda*sizeof(double));
    if((*io_z) == NULL || (*io_lambda) == NULL || (*io_tau) == NULL || (*io_mu) == NULL ||
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

double SDSS_P_lambda_obs_given_true_lambda(double observed_lambda, double true_lambda, double zz) 
{
  static int first = 0;
  static gsl_spline2d* ftau;
  static gsl_spline2d* fmu;
  static gsl_spline2d* fsigma;
  static gsl_spline2d* ffmask;
  static gsl_spline2d* ffprj;

  if(first == 0) 
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

    z = (double**) malloc(1*sizeof(double*));
    lambda = (double**) malloc(1*sizeof(double*));
    tmp_tau = (double**) malloc(1*sizeof(double*));
    tmp_mu = (double**) malloc(1*sizeof(double*));
    tmp_sigma = (double**) malloc(1*sizeof(double*));
    tmp_fmask = (double**) malloc(1*sizeof(double*));
    tmp_fprj = (double**) malloc(1*sizeof(double*));
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
    tau = (double*) malloc(nz*nlambda*sizeof(double));
    mu = (double*) malloc(nz*nlambda*sizeof(double));
    sigma = (double*) malloc(nz*nlambda*sizeof(double));
    fmask = (double*) malloc(nz*nlambda*sizeof(double));
    fprj = (double*) malloc(nz*nlambda*sizeof(double));
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

double SDSS_int_P_lambda_obs_given_M_nointerp(double true_lambda, void* params)
{
  double* ar = (double*) params;
  
  const double M = ar[0];
  const double z = ar[1];
  const double obs_lambda = ar[2];

  const double r1 = SDSS_P_lambda_obs_given_true_lambda(obs_lambda, true_lambda, z);
  const double r2 = SDSS_P_true_lambda_given_mass(true_lambda, M, z);

  return r1*r2;
}

double SDSS_P_lambda_obs_given_M_nointerp(double obs_lambda, void* params)
{
  double* ar = (double*) params; // mass , redshift 
  
  const double M = ar[0];
  const double z = ar[1];
  
  double params_in[3] = {M, z, obs_lambda}; 
  
  const double true_lambda_min =  3.;
  const double true_lambda_max =  160.;

  return int_gsl_integrate_medium_precision(SDSS_int_P_lambda_obs_given_M_nointerp, 
    (void*) params_in, true_lambda_min, true_lambda_max, NULL, GSL_WORKSPACE_SIZE);
}

// \int_(bin_lambda_obs_min)^(bin_lambda_obs_max) \dlambda_obs P(\lambda_obs|M)
// (see for example https://arxiv.org/pdf/1810.09456.pdf - eq3 qnd 6) 
double SDSS_binned_P_lambda_obs_given_M_nointerp(int nl, double M, double z)
{
  double params[2] = {M, z};
  
  const double bin_lambda_obs_min = Cluster.N_min[nl];
  const double bin_lambda_obs_max = Cluster.N_max[nl];
  
  return int_gsl_integrate_medium_precision(SDSS_P_lambda_obs_given_M_nointerp, (void*) params, 
    bin_lambda_obs_min, bin_lambda_obs_max, NULL, GSL_WORKSPACE_SIZE);
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// INTERFACE
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

// \int_(bin_lambda_obs_min)^(bin_lambda_obs_max) \dlambda_obs P(\lambda_obs|M)
// (see for example https://arxiv.org/pdf/1810.09456.pdf - eq3 qnd 6) 
double binned_P_lambda_obs_given_M_nointerp(int nl, double M, double z)
{
  if (strcmp(Cluster.model, "SDSS") == 0) 
  {
    return SDSS_binned_P_lambda_obs_given_M_nointerp(nl, M, z);
  }
  else if (strcmp(Cluster.model, "BUZZARD") == 0)
  {
    return buzzard_binned_P_lambda_obs_given_M_nointerp(nl, M, z);
  }
  else
  {
    log_fatal("option not implemented");
    exit(1);
  }
}

// Threading this function is not allowed (unless nl is fixed)
double binned_P_lambda_obs_given_M(int nl, double M, double z) 
{
  static int BIN_OBS_LAMBDA = -1; // lambda = bin number of observed lambda
  static cosmopara C;
  static nuisancepara N;
  static double** table_P1 = 0;

  const double logM = log10(M);
  const double logmmin = log_M_min;
  const double logmmax = log_M_max; 
  const double zmin = 0.20; 
  const double zmax = 0.80;
  const int nz = 10;
  const int nm = 50;
  const double dm = (logmmax-logmmin)/(nm-1.);
  const double dz = (zmax-zmin)/(nz-1.);
  const double mmin = pow(10.0, logmmin);
  const double mmax = pow(10.0, logmmax);

  if (table_P1 == 0)
  {
    table_P1 = create_double_matrix(0, nz-1, 0, nm-1);
  }
  if ((BIN_OBS_LAMBDA != nl) || recompute_clusters(C, N))
  { 
    {
      const int i = 0;
      {
        const int j = 0;
        const double M = pow(10.0, logmmin + j*dm);
        table_P1[i][j] = binned_P_lambda_obs_given_M_nointerp(nl, M, z);
      } 
      #pragma omp parallel for     
      for (int j=1; j<nm; j++) 
      {
        const double M = pow(10.0, logmmin + j*dm);
        table_P1[i][j] = binned_P_lambda_obs_given_M_nointerp(nl, M, z);
      }
    }
    #pragma omp parallel for
    for (int i=1; i<nz; i++) 
    {
      for (int j=0; j<nm; j++) 
      {
        const double M = pow(10.0, logmmin + j*dm);
        table_P1[i][j] = binned_P_lambda_obs_given_M_nointerp(nl, M, z);
      }
    }
    
    update_cosmopara(&C);
    update_nuisance(&N);
    BIN_OBS_LAMBDA = nl;
  }
  if (z > zmin && z < zmax && M > mmin && M < mmax)
  {
    return interpol2d(table_P1, nz, zmin, zmax, dz, z, nm, logmmin, logmmax, dm, logM, 1, 1);
  }
  return 0.;
}

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// CLUSTER BIAS
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

double int_dndlogM_times_binned_P_lambda_obs_given_M(double logM, void* params)
{
   double *ar = (double*) params; //n_lambda, z
   
   const int nl = (int) ar[0];
   const double z = ar[1];

   const double M = exp(logM);
   const double a = 1./(1. + z);
   
   return massfunc(M, a)*M*binned_P_lambda_obs_given_M(nl, M , z);
}

double int_weighted_bias_nointerp(double logM, void* params)
{
  double *ar = (double*) params; //nl, z
  
  const double z = ar[1];
  
  const double M = exp(logM);  
  const double a = 1./(1. + z);
  
  return B1(M, a)*int_dndlogM_times_binned_P_lambda_obs_given_M(logM, params)*
    nuisance.cluster_selection[0]; 
}

double weighted_bias_nointerp(int nl, double z)
{
  double param[2] = {(double) nl, z};
  
  const double mmin = log_M_min/0.4342944819; // log(pow(10.,12.)) --- 0.4342944819 = log10(e)
  const double mmax = log_M_max/0.4342944819; // log(pow(10.,15.9)) --- 0.4342944819 = log10(e)
  
  const double r1 = int_gsl_integrate_low_precision(int_weighted_bias_nointerp, (void*) param, mmin, 
    mmax, NULL, GSL_WORKSPACE_SIZE);
  
  const double r2 = int_gsl_integrate_low_precision(int_dndlogM_times_binned_P_lambda_obs_given_M, 
    (void*) param, mmin, mmax, NULL, GSL_WORKSPACE_SIZE);

  return r1/r2; 
}

double weighted_bias(int nl, double z)
{
  static cosmopara C;
  static nuisancepara N;
  static double** table;

  const int nasize = Ntable.N_a;
  const int nlsize = Cluster.N200_Nbin;
  const double zmin = fmax(tomo.cluster_zmin[0]-0.05, 0.01); 
  const double zmax = tomo.cluster_zmax[tomo.cluster_Nbin-1] + 0.05;
  const double amin = 1./(1.+zmax); 
  const double amax = 1./(1.+zmin);
  const double da = (amax-amin)/(nasize-1.);

  if (table == 0) 
  {
    table = create_double_matrix(0, nlsize - 1, 0, nasize);
  }

  if (recompute_clusters(C, N))
  {
    {
      const int n = 0;
      {
        const int i = 0;
        const double a = amin + i*da;
        table[n][i] = weighted_bias_nointerp(n, 1./a - 1.);
      }
      #pragma omp parallel for
      for (int i=1; i<nasize; i++)
      {
        const double a = amin + i*da;
        table[n][i] = weighted_bias_nointerp(n, 1./a - 1.);
      }
    }
    #pragma omp parallel for
    for (int n=1; n<nlsize; n++)
    {
      for (int i=0; i<nasize; i++)
      {
        const double a = amin + i*da;
        table[n][i] = weighted_bias_nointerp(n, 1./a - 1.);
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }

  const double a = 1./(z+1.);
  return interpol(table[nl], nasize, amin, amax, da, a, 0., 0.);
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

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// CLUSTER POWER SPECTRUM AND CROSS SPECTRUM (WITH GALAXIES AND MATTER)
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

// nl{1,2} = lambda_obs bins, n{i,j} = cluster redshift bins
double binned_p_cc(double k, double a, int nl1, int nl2, int use_linear_ps)
{ // binned in lambda_obs (nl1, nl2 = lambda_obs bin (cluster 1 and 2))
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  const double z = 1./a - 1;
  const double cluster_bias1 = weighted_bias(nl1, z); 
  const double PK = use_linear_ps == 1 ? p_lin(k, a) : Pdelta(k, a);

  if (nl1 == nl2)
  {
    return cluster_bias1*cluster_bias1*PK;
  }
  else
  {
    const double cluster_bias2 = weighted_bias(nl2, z); 
    return cluster_bias1*cluster_bias2*PK;
  }
}

// nl = lambda_obs bin, nj = galaxy redshift bin
double binned_p_cg(double k, double a, int nl, int nj, int use_linear_ps)
{ // binend in lambda_obs
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  const double z = 1.0/a - 1.0;
  const double cluster_bias = weighted_bias(nl, z);
  const double PK = use_linear_ps == 1 ? p_lin(k, a) : Pdelta(k, a);
  return cluster_bias * gbias.b1_function(z, nj)* PK;
}

// nl = lambda_obs bin, ni = cluster redshift bin, nj = source redshift bin
double int_binned_p_cm(double logM, void* params)
{ // binned in lambda_obs (nl = lambda_obs bin)
  double* ar = (double*) params;
   
  const int nl = (int) ar[0];
  const double z = ar[1];
  const double k = ar[2];
  const int use_linear_ps = (int) ar[3];
  
  const double M = exp(logM); 
  const double a = 1./(1 + z);

  const double PCM_1h = M/(cosmology.rho_crit*cosmology.Omega_m)*u_nfw_c(conc(M, a), k, M, a); 
   
  const double PCM_2h = use_linear_ps == 1 ? B1(M, a)*p_lin(k, a) : B1(M, a)*Pdelta(k, a);
   
  return M*massfunc(M, a)*(PCM_1h + PCM_2h)*binned_P_lambda_obs_given_M(nl, M , z);
}

double binned_p_cm(double k, double a, int nl, int use_linear_ps)
{
  if(!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  const double z = 1.0/a - 1;
  double params[4] = {(double) nl, z, k, (double) use_linear_ps};
  const double mmin = log_M_min/0.4342944819; // log(pow(10.,12.)) --- 0.4342944819 = log10(e)
  const double mmax = log_M_max/0.4342944819; // log(pow(10.,15.9)) --- 0.4342944819 = log10(e)  
  const double ncounts = binned_average_number_counts(nl, z);

  return int_gsl_integrate_low_precision(int_binned_p_cm, params, mmin, mmax, NULL, 
    GSL_WORKSPACE_SIZE)/ncounts; 
}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// Area
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

double get_area(double z)
{
  static double *areatab =0;
  static double*  z_tab_info =0;
  static int nz_tab = 850;
  if (areatab==0)
  {
     double* ztab;
     FILE *Farea;
     char filename[400];
     sprintf(filename,"/home/users/chto/code/lighthouse/tables/area_z.txt");
     Farea = fopen(filename,"r");
     ztab = create_double_vector(0, nz_tab);
     areatab = create_double_vector(0, nz_tab);
     z_tab_info = create_double_vector(0,3);
     for(int i=0; i<nz_tab; ++i){
        fscanf(Farea, "%lf", &ztab[i]);
        fscanf(Farea, "%lf", &areatab[i]);
     }
     fclose(Farea);
     z_tab_info[0] = ztab[0];
     z_tab_info[1] = ztab[nz_tab-1];
     z_tab_info[2] = ztab[1]-ztab[0];
     free_double_vector(ztab,0,nz_tab-1);
  }
  return interpol(areatab, nz_tab, z_tab_info[0], z_tab_info[1], z_tab_info[2], z, 1.0, 1.0);
}