#include <assert.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_integration.h>

#include "halo.h"
#include "basics.h"
#include "cosmo3D.h"
#include "redshift_spline.h"
#include "structs.h"
#include "recompute.h"

#include "log.c/src/log.h"

static int GSL_WORKSPACE_SIZE = 250;

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// halo model routines to convert scales/radii to halo masses, and vice versa
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double delta_c(double a __attribute__((unused))) 
{ // set to 1.686 for consistency with Tinker mass & bias definition
  return 1.686; 
}

double delta_Delta(double a __attribute__((unused))) 
{ // using Tinker et al mass & bias functions with \Delta = 200 rho_{mean}
  return 200.0;
}

double rho_Delta(const double a) // virial density in solar masses/h/(H0/c)^3 
{ // using Tinker et al mass & bias functions with \Delta = 200 rho_{mean}
  return (delta_Delta(a) * cosmology.rho_crit * cosmology.Omega_m);
}

double m_Delta(const double r, const double a) 
{
  return (4.0*M_PI /3.0)*r*r*r*rho_Delta(a);
}

double r_Delta(const double m, const double a) 
{ // calculate r_Delta in c/H0 given m in (solar masses/h)
  return pow(3./(4.0*M_PI)*(m/rho_Delta(a)), 1./3.);
}

double r_s(const double m, const double a) 
{ 
  return r_Delta(m, a) / conc(m, a); 
}

double radius(const double m) 
{
  return pow(3.0/4.0*m/(M_PI*cosmology.rho_crit*cosmology.Omega_m), 1.0/3.0);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Variance of the density field
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_sigma2(double x, void* params) // inner integral
{ // refactored FT of spherical top-hat to avoid numerica divergence of 1/x
  double* ar = (double*) params;
  const double k = x/ar[0];
  const double PK = p_lin(k, 1.0);
  
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

double sigma2_nointerp(const double M, const int init_static_vars_only) 
{
  double ar[1] = {radius(M)};
  const double xmin = 0;
  const double xmax = 14.1;

  return (init_static_vars_only == 1) ? int_for_sigma2((xmin+xmax)/2.0, (void*) ar) :
    Ntable.high_def_integration > 0 ?
    int_gsl_integrate_high_precision(int_for_sigma2, (void*) ar, xmin, xmax, 
      NULL, GSL_WORKSPACE_SIZE) :
    int_gsl_integrate_medium_precision(int_for_sigma2, (void*) ar, xmin, xmax, 
      NULL, GSL_WORKSPACE_SIZE);
}

double sigma2(const double M) 
{
  static cosmopara C;
  static double* table;
  
  const int nlnM = Ntable.N_M;
  const double lnMmin = log(limits.M_min);
  const double lnMmax = log(limits.M_max);
  const double dlnM = (lnMmax - lnMmin)/((double) nlnM - 1.0);
  
  if (table == 0) 
  {
    table = (double*) malloc(sizeof(double)*nlnM);
  } 
  if (recompute_cosmo3D(C)) 
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init = sigma2_nointerp(exp(lnMmin), 1);
    }
    #pragma GCC diagnostic pop
    #pragma omp parallel for
    for (int i=0; i<nlnM; i++) 
    {
      table[i] = log(sigma2_nointerp(exp(lnMmin + i*dlnM), 0));
    }
    update_cosmopara(&C);
  }

  const double lnM = log(M);
  if (lnM < lnMmin)
  {
    log_warn("M = %e < l_min = %e. Extrapolation adopted", M, exp(lnMmin));
  }
  if (lnM > lnMmax)
  {
    log_warn("M = %e > l_max = %e. Extrapolation adopted", M, exp(lnMmax));
  }
  return exp(interpol(table, nlnM, lnMmin, lnMmax, dlnM, lnM, 1.0, 1.0));
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// \nu
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double nu(const double M, const double a) 
{
  static int init = 1;
  if ((cosmology.Omega_nu > 0) && init) 
  {
    log_fatal("halo.c does not support cosmologies with massive neutrinos");
    exit(1);
    init = 0;
  }
  return delta_c(a)/(sqrt(sigma2(M))*growfac(a));
}

double nuv2(const double M, const double a, const double growfac) 
{
  static int init = 1;
  if ((cosmology.Omega_nu > 0) && init) 
  {
    log_fatal("halo.c does not support cosmologies with massive neutrinos");
    exit(1);
    init = 0;
  }
  return delta_c(a)/(sqrt(sigma2(M))*growfac);
}

double int_for_dlognudlogm(double lnM, void* params __attribute__((unused))) 
{
  return log(nu(exp(lnM), 1.0));
}

double dlognudlogm(const double M) 
{
  static cosmopara C;
  static double* table;

  const int nlnM = Ntable.N_M;
  const double lnMmin = log(limits.M_min);
  const double lnMmax = log(limits.M_max);
  const double dlnM = (lnMmax - lnMmin)/((double) nlnM - 1.0);

  if (table == 0) 
  {
    table = (double*) malloc(sizeof(double) * nlnM);
  }
  if (recompute_cosmo3D(C)) 
  {
    {
      const int i = 0;
      double result, abserr;
      gsl_function F;
      F.function = &int_for_dlognudlogm;
      F.params = 0;

      int status = gsl_deriv_central(&F, lnMmin + i*dlnM, 0.1*(lnMmin + i*dlnM), &result, &abserr);
      if (status) 
      {
        log_fatal(gsl_strerror(status));
        exit(1);
      }
      table[i] = result;
    }
    #pragma omp parallel for
    for (int i=1; i <nlnM; i++) 
    {
      double result, abserr;
      gsl_function F;
      F.function = &int_for_dlognudlogm;
      F.params = 0;

      int status = gsl_deriv_central(&F, lnMmin + i*dlnM, 0.1*(lnMmin + i*dlnM), &result, &abserr);
      if (status) 
      {
        log_fatal(gsl_strerror(status));
        exit(1);
      }
      table[i] = result;
    }
    update_cosmopara(&C);
  }
  
  const double lnM = log(M);
  if (lnM < lnMmin)
  {
    log_warn("M = %e < M_min = %e. Extrapolation adopted", M, exp(lnMmin));
  }
  if (lnM > lnMmax)
  {
    log_warn("M = %e > M_max = %e. Extrapolation adopted", M, exp(lnMmax));
  }  
  return interpol(table, nlnM, lnMmin, lnMmax, dlnM, lnM, 1.0, 1.0);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// mass function & halo bias (Tinker et al.)
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double f_tinker(const double n, const double a_in) 
{ 
  // Eqs. (8-12) + Table 4 from Tinker et al. 2010
  // aa = alpha, b = beta, c = gamma, f = phi, e = eta
  double a = fmax(0.25, a_in); // limit fit range of mass function evolution to z <= 3, as
                               // discussed after Eq. 12 of http://arxiv.org/pdf/1001.3162v2.pdf
  const double aa = 0.368;
  const double b = 0.589 * pow(a, -0.2);
  const double c = 0.864 * pow(a, 0.01);
  const double f = -0.729 * pow(a, .08);
  const double e = -0.243 * pow(a, -0.27);
  return aa * (1.0 + pow(b * n, -2.0 * f)) * pow(n, 2.0 * e) * exp(-c * n * n / 2.0);
}

double fnu_tinker(const double n, const double a) 
{ 
  return f_tinker(n, a) * n; 
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// norm for redshift evolution of < function (Eq.8 in Tinker 2010) + enforcing <M> to be unbiased
double int_for_mass_norm(double x, void* params) 
{ // inner integral
  double* ar = (double*) params;
  
  const double a = ar[0];
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  
  return f_tinker(x, a);
}

double mass_norm_nointerp(const double a, const int init_static_vars_only) 
{
  double ar[1] = {a}; 
  const double xmin = 0.0;
  const double xmax = 10.0;

  return (init_static_vars_only == 1) ? int_for_mass_norm((xmin+xmax)/2.0, (void*) ar) :
    Ntable.high_def_integration > 0 ?
    int_gsl_integrate_high_precision(int_for_mass_norm, (void*) ar, xmin, xmax, 
      NULL, GSL_WORKSPACE_SIZE) :
    int_gsl_integrate_medium_precision(int_for_mass_norm, (void*) ar, xmin, xmax, 
      NULL, GSL_WORKSPACE_SIZE);
}

double mass_norm(const double a) 
{
  static cosmopara C;
  static double* table;
  
  const int na = Ntable.N_a;
  const double amin = limits.a_min;
  const double amax = 0.9999999;
  const double da = (amax - amin) / ((double) na - 1.0);
  
  if (table == 0) 
  {
    table = (double*) malloc(sizeof(double)*na);
  }
  if (recompute_cosmo3D(C)) 
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init = mass_norm_nointerp(amin, 1);
    }
    #pragma GCC diagnostic pop
    #pragma omp parallel for
    for (int i=0; i<na; i++)
    {   
      table[i] = mass_norm_nointerp(amin + i*da, 0);
    }
    update_cosmopara(&C);
  }
  if (a < amin)
  {
    log_warn("a = %e < a_min = %e. Extrapolation adopted", a, amin);
  }
  return interpol(table, na, amin, amax, da, fmin(a, amax-da), 1.0, 1.0);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double B1_nu(const double n, const double a) 
{
  // Eq (6) + Table 2 from Tinker et al. 2010
  const double y = log10(200.0);
  const double A = 1.0 + 0.24 * y * exp(-pow(4.0 / y, 4.0));
  const double aa = 0.44 * y - 0.88;
  const double B = 0.183;
  const double b = 1.5;
  const double C = 0.019 + 0.107 * y + 0.19 * exp(-pow(4.0 / y, 4.0));
  const double c = 2.4;
  return 1.0 - A * pow(n, aa) / (pow(n, aa) + pow(delta_c(a), aa)) +
         B * pow(n, b) + C * pow(n, c);
}

double int_for_bias_norm(double x, void* params) 
{ // correct bias for to halo mass cuts (so large-scale 2-h matches PT results at all zs)
  double* ar = (double*) params;
  
  double a = ar[0];
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  
  return B1_nu(x, a) * f_tinker(x, a);
}

double bias_norm_nointerp(const double a, const int init_static_vars_only)
{
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }

  double ar[1] = {a};
  const double numin = nu(limits.M_min, a);
  const double numax = nu(limits.M_max, a);

  return (init_static_vars_only == 1) ? int_for_bias_norm((numin+numax)/2.0, (void*) ar) :
    Ntable.high_def_integration > 0 ?
    int_gsl_integrate_high_precision(int_for_bias_norm, (void*) ar, numin, numax, 
      NULL, GSL_WORKSPACE_SIZE) :
    int_gsl_integrate_medium_precision(int_for_bias_norm, (void*) ar, numin, numax, 
      NULL, GSL_WORKSPACE_SIZE);
}

double bias_norm(const double a) 
{
  static cosmopara C;
  static double* table;

  const int na = Ntable.N_a;
  const double amin = limits.a_min; 
  const double amax = 0.9999999;
  const double da = (amax - amin) / ((double) na - 1.0);

  if (table == 0) 
  {
    table = (double*) malloc(sizeof(double)*na);
  }
  if (recompute_cosmo3D(C)) 
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init = bias_norm_nointerp(amin, 1); 
    }
    #pragma GCC diagnostic pop
    #pragma omp parallel for
    for (int i=0; i<na; i++)
    {     
      table[i] = bias_norm_nointerp(amin + i*da, 0);
    }
    table[na - 1] = 1.0;
    update_cosmopara(&C);
  }
  return interpol(table, na, amin, amax, da, fmin(a, amax - da), 1.0, 1.0);
}

double massfunc(const double m, const double a) 
{
  return fnu_tinker(nu(m, a), a) * cosmology.rho_crit * cosmology.Omega_m/m/m * dlognudlogm(m);
}

double B1(const double m, const double a) 
{ // b(m,a) based on redshift evolution fits in Tinker et al. paper, no additional normalization
  return B1_nu(nu(m, a), a); 
}

double B1_normalized(const double m, const double a) 
{ // divide by bias norm only in matter spectra, not in HOD modeling/cluster analyses
  return B1_nu(nu(m, a), a)/bias_norm(a);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// halo profiles
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double conc(const double m, const double a) 
{ // mass-concentration relation: Bhattacharya et al. 2013, Delta = 200 rho_{mean} (Table 2)
  const double g = growfac(a);
  return 9.0*pow(nuv2(m, a, g), -0.29)*pow(g, 1.15); 
}

double u_nfw_c(const double c, const double k, const double m, const double a) 
{ // analytic FT of NFW profile, from Cooray & Sheth 01
  const double x = k * r_Delta(m, a) / c;
  const double xu = (1. + c) * x;

  gsl_sf_result SI_XU;
  int status = gsl_sf_Si_e(xu, &SI_XU);
  if (status) 
  {
    log_fatal(gsl_strerror(status));
    exit(1);
  }

  gsl_sf_result SI_X;
  {
    int status = gsl_sf_Si_e(x, &SI_X);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
  }

  gsl_sf_result CI_XU;
  {
    int status = gsl_sf_Ci_e(xu, &CI_XU);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
  }

  gsl_sf_result CI_X;
  {
    int status = gsl_sf_Ci_e(x, &CI_X);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
  }

  return (sin(x)*(SI_XU.val - SI_X.val) - sinl(c*x)/xu + cos(x)*(CI_XU.val - CI_X.val))/
    (log(1. + c) - c/(1. + c));
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// matter power spectrum
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_I02(double logm, void* params) 
{
  double* ar = (double*) params;
  const double a = ar[0];
  const double k1 = ar[1];
  const double k2 = ar[2];
  const double m = exp(logm);
  const double c = conc(m, a);
  const double u = u_nfw_c(c, k1, m, a)*u_nfw_c(c, k2, m, a);
  const double tmp = m/(cosmology.rho_crit * cosmology.Omega_m);
  return massfunc(m, a) * m * u * tmp * tmp;
}

double I02_nointerp(const double k1, const double k2, const double a, const int init_static_vars_only) 
{
  double ar[3] = {a, k1, k2};
  const double lnMmin = log(limits.M_min);
  const double lnMmax = log(limits.M_max);

  return init_static_vars_only == 1 ? int_for_I02(lnMmin, (void*) ar) :
    Ntable.high_def_integration > 0 ?
    int_gsl_integrate_high_precision(int_for_I02, (void*) ar, lnMmin, lnMmax, 
      NULL, GSL_WORKSPACE_SIZE) :
    int_gsl_integrate_medium_precision(int_for_I02, (void*) ar, lnMmin, lnMmax, 
      NULL, GSL_WORKSPACE_SIZE);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_I11(double logm, void* params) 
{
  const double* ar = (double*) params;  
  const double a = ar[0];
  const double k = ar[1];
  const double m = exp(logm);
  const double c = conc(m, a);
  const double u = u_nfw_c(c, k, m, a);
  return massfunc(m, a) * m * m * u * B1_normalized(m, a)/(cosmology.rho_crit * cosmology.Omega_m);
}

double I11_nointerp(const double k, const double a, const int init_static_vars_only) 
{ 
  double ar[2] = {a, k};
  const double lnMmin = log(limits.M_min);
  const double lnMmax = log(limits.M_max);
  
  return init_static_vars_only == 1 ? int_for_I11(lnMmin, (void*) ar) :
    Ntable.high_def_integration > 0 ?
    int_gsl_integrate_high_precision(int_for_I11, (void*) ar, lnMmin, lnMmax, 
      NULL, GSL_WORKSPACE_SIZE) :
    int_gsl_integrate_medium_precision(int_for_I11, (void*) ar, lnMmin, lnMmax, 
      NULL, GSL_WORKSPACE_SIZE);
}

double I11(const double k, const double a) 
{
  static cosmopara C;
  static double** table = 0;

  const int na = Ntable.N_a;
  const double amin = limits.a_min;
  const double amax = 0.999999;
  const double da = (amax - amin) / ((double) na - 1.0);

  const int nk = Ntable.N_k_nlin;
  const double lnkmin = log(limits.k_min_cH0);
  const double lnkmax = log(limits.k_max_cH0);
  const double dlnk = (lnkmax - lnkmin) / ((double) nk - 1.0);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*Ntable.N_a);
    for(int i=0; i<Ntable.N_a;i++)
    {
      table[i] = (double*) malloc(sizeof(double)*Ntable.N_k_nlin);
    }
  }
  if (recompute_cosmo3D(C)) 
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init = I11_nointerp(exp(lnkmin), amin, 1);
    }
    #pragma GCC diagnostic pop
    #pragma omp parallel for collapse(2)
    for (int i=0; i<na; i++)
    {
      for (int j=0; j<nk; j++)
      {
        table[i][j] = I11_nointerp(exp(lnkmin + j*dlnk), amin + i*da, 0);
      }
    }
    update_cosmopara(&C);
  }

  const double lnk = log(k);
  if (lnk < lnkmin)
  {
    log_warn("k = %e < k_min = %e. Extrapolation adopted", k, exp(lnkmin));
  }
  if (lnk > lnkmax)
  {
    log_warn("k = %e > k_max = %e. Extrapolation adopted", k, exp(lnkmax));
  }
  if (a < amin)
  {
    log_warn("a = %e < a_min = %e. Extrapolation adopted", a, amin);
  }
  if (a > amax)
  {
    log_warn("a = %e > a_max = %e. Extrapolation adopted", a, amax);
  }
  return exp(interpol2d(table, na, amin, amax, da, a, nk, lnkmin, lnkmax, dlnk, lnk, 1.0, 1.0));
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------



// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// simplistic abundance matching to get approximate b(z) of source galaxies
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double n_s_cmv(double a) 
{
  double dV_dz = pow(f_K(chi(a)), 2.0) / hoverh0(a); // comoving dV/dz(z = 1./a-1) per radian^2
  // dN/dz/radian^2/(dV/dz/radian^2)
  return zdistr_photoz(1.0/a - 1., -1) * survey.n_gal * survey.n_gal_conversion_factor / dV_dz; 
}

double int_for_ngmatched(double logm, void* params) 
{
  double* ar = (double*) params;
  const double a = ar[0];
  const double m = exp(logm);
  return massfunc(m, a) * m;
}

double int_for_b_ngmatched(double logm, void* params) 
{
  double* ar = (double*) params;
  const double a = ar[0];
  const double m = exp(logm);
  return massfunc(m, a) * m * B1(exp(logm), a);
}

double b_ngmatched(double a, double n_cmv) 
{
  double ar[1] = {a};
  const double dm = 0.1;
  
  double mlim = log(5.e+14);
  
  double n = int_gsl_integrate_medium_precision(int_for_ngmatched, (void*) ar, mlim, 
    log(limits.M_max), NULL, GSL_WORKSPACE_SIZE);
  
  while (n < n_cmv) 
  {
    mlim -= dm;
    n += int_gsl_integrate_medium_precision(int_for_ngmatched, (void*) ar, mlim, mlim + dm, 
      NULL, GSL_WORKSPACE_SIZE);
  }

  return int_gsl_integrate_medium_precision(int_for_b_ngmatched, (void*) ar, mlim, 
    log(limits.M_max), NULL, GSL_WORKSPACE_SIZE) /n;
}

double b_source(double a) // lookup table for b1 of source galaxies
{ 
  static cosmopara C;
  static double* table;

  const double amin = limits.a_min;
  const double amax = 1.0 / (1.0 + redshift.shear_zdist_zmin_all);
  const double da = (amax - amin) / ((double) Ntable.N_a - 1.0);

  if (table == 0)
  {
    table = (double*) malloc(sizeof(double)*Ntable.N_a);
  }
  if (recompute_cosmo3D(C))
  {
    {
      const int i = 0;
      table[i] = b_ngmatched(amin + i*da, n_s_cmv(amin + i*da));
    }
    #pragma omp parallel for
    for (int i=1; i<Ntable.N_a; i++) 
    {
      table[i] = b_ngmatched(amin + i*da, n_s_cmv(amin + i*da));
    }
    update_cosmopara(&C);
  }
  return interpol(table, Ntable.N_a, amin, amax, da, a, 1.0, 1.0);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// functions important for y correlation functions
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// FT of bound gas K-S Profile (sinc integral part)
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_F_KS(double x, void* params)
{
  double* ar = (double*) params;
  const double y = ar[0];
  return (x*sinl(y*x)/y)*pow(log(1.0 + x)/x, ynuisance.gas_Gamma_KS/(ynuisance.gas_Gamma_KS - 1.0));
}

double F_KS_nointerp(double c, double krs, const int init_static_vars_only) 
{
  double ar[1] ={krs};
  const double xmin = 0.0;
  const double xmax = c;
  return (init_static_vars_only == 1) ? int_F_KS((xmin + xmax)/2.0, (void*) ar) :
    Ntable.high_def_integration > 0 ?
    int_gsl_integrate_high_precision(int_F_KS, (void*) ar, xmin, xmax, NULL, GSL_WORKSPACE_SIZE) :
    int_gsl_integrate_medium_precision(int_F_KS, (void*) ar, xmin, xmax, NULL, GSL_WORKSPACE_SIZE);
}

double int_F0_KS(double x, void* params __attribute__((unused)))
{
  return x*x*pow(log(1.0 + x)/x, 1.0/(ynuisance.gas_Gamma_KS - 1.0));
}

double F0_KS_nointerp(double c, const int init_static_vars_only)
{
  double ar[1] = {0.0};
  const double xmin = 0.0;
  const double xmax = c;
  return (init_static_vars_only == 1) ? int_F0_KS((xmin + xmax)/2.0, (void*) ar) :
    Ntable.high_def_integration > 0 ?
    int_gsl_integrate_high_precision(int_F0_KS, (void*) ar, xmin, xmax, NULL, GSL_WORKSPACE_SIZE):
    int_gsl_integrate_medium_precision(int_F0_KS, (void*) ar, xmin, xmax, NULL, GSL_WORKSPACE_SIZE);
}

double u_KS(double c, double k, double rv)
{ // TODO: IMPLEMENT FFT VERSION WITH OPENMP
  static ynuisancepara N;
  static double** table = 0;
  static double* norm = 0;
 
  const int N_c = Ntable.halo_uKS_nc;
  const double cmin = limits.halo_uKS_cmin; 
  const double cmax = limits.halo_uKS_cmax;
  const double dc = (cmax - cmin)/((double) N_c - 1.0);

  const int N_x = Ntable.halo_uks_nx;
  const double xmin = limits.halo_uKS_xmin; 
  const double xmax = limits.halo_uKS_xmax; // full range of possible k*R_200/c in the code
  const double lnxmin = log(xmin);
  const double lnxmax = log(xmax); 
  const double dx = (lnxmax - lnxmin)/((double) N_x - 1.0); 
 
  if (table == 0) 
  {    
    table = (double**) malloc(sizeof(double*)*N_c);
    for (int i=0; i<N_c; i++) 
    {
      table[i] = (double*) malloc(sizeof(double)*N_x);
    }
    norm = (double*) malloc(sizeof(double)*N_c);
  }
  if (recompute_yhalo(N)) // Only ynuisance.gas_Gamma is used in creating the table
  { 
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init = F0_KS_nointerp(cmin, 1);
      double init2 = F_KS_nointerp(cmin, exp(lnxmin), 1);
    }
    #pragma GCC diagnostic pop
    #pragma omp parallel for
    for (int i=0; i<N_c; i++) 
    {
      norm[i] = F0_KS_nointerp(cmin + i*dc, 0);
    }
    #pragma omp parallel for collapse(2)
    for (int i=0; i<N_c; i++)
    {
      for (int j=0; j<N_x; j++) 
      {
        table[i][j] = F_KS_nointerp(cmin + i*dc, exp(lnxmin + j*dx), 0)/norm[i];
      }
    }
    update_ynuisance(&N);
  }
  const double x = k * rv/c;
  const double lnx = log(x);
  if (c < cmin) 
  {
    log_warn("c = %e < c_min = %e. Extrapolation adopted", c, cmin);
    return 0.;
  }
  if (c > cmax) 
  {
    log_warn("c = %e > c_max = %e. Extrapolation adopted", c, cmax);
    return 0.;
  }
  if (lnx < lnxmin) 
  {
    log_warn("x = %e < x_min = %e. Extrapolation adopted", x, exp(lnxmin));
    return 0.;
  }
  if (lnx > lnxmax) 
  {
    log_warn("x = %e > x_max = %e. Extrapolation adopted", x, exp(lnxmax));
    return 0.;
  }
  return interpol2d(table, N_c, cmin, cmax, dc, c, N_x, lnxmin, lnxmax, dx, lnx, 0.0, 0.0);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double frac_bnd(double M)
{
  const double M0 = pow(10.0, ynuisance.gas_lgM0);
  return cosmology.Omega_b/(cosmology.Omega_m*(1.0+ pow(M0/M, ynuisance.gas_beta)));
}

double frac_ejc(double M)
{
  const double logM = log10(M);
  const double delta = (logM - ynuisance.gas_lgM_star)/ynuisance.gas_sigma_star;
  const double tmp = ynuisance.gas_A_star * exp(-0.5*delta*delta);  
  const double frac_star = ((logM > ynuisance.gas_lgM0) && (tmp < ynuisance.gas_A_star/3.0)) ?
    ynuisance.gas_A_star/3.0 : tmp; 
  return frac_bnd(M) - frac_star;
}

double u_y_bnd(double c, double k, double m, double a)
{ //unit: [G(M_solar/h)^2 / (c/H0)]
  const double rv = r_Delta(m,a);
  const double mu_p = 4.0/(3.0 + 5*ynuisance.gas_f_H);
  const double mu_e = 2.0/(1.0 + ynuisance.gas_f_H);
  return (2.0*ynuisance.gas_alpha/(3.0*a))*(mu_p/mu_e)*frac_bnd(m)*m*(m/rv)*u_KS(c, k, rv);
}

double u_y_ejc(double m)
{ // [m] = [Msun/h]
  const double num_p = 1.1892e57; // proton number in 1 solar mass, unit [1/Msun]
  // convert ejected gas temperature (K) to energy eV then to unit [G (Msun/h)^2 / (c/H0) * h]
  // double E_w = pow(10,ynuisance.gas_lgT_w) * 8.6173e-5 * 1.81945e-68;
  const double E_w = pow(10,ynuisance.gas_lgT_w) * 8.6173e-5 * 5.616e-44;
  const double mu_e = 2./(1.+ynuisance.gas_f_H);
  return (num_p * m * frac_ejc(m) / mu_e) * E_w; // final unit in [G(Msun/h)^2 / (c/H0)]
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// matter-y and yy power spectra
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_I11_y(double logm, void* params) 
{
  const double* ar = (double*) params;  
  const double a = ar[0];
  const double k = ar[1];
  const double m = exp(logm);
  const double c = conc(m,a); 
  const double u = u_y_bnd(c, k, m, a) + u_y_ejc(m);
  const double tmp = m/(cosmology.rho_crit*cosmology.Omega_m);
  return massfunc(m, a) * m * u * tmp * B1_normalized(m, a);
}

double I11_y_nointerp(const double k, const double a, const int init_static_vars_only) 
{ 
  double ar[2] = {a, k};
  const double lnMmin = log(limits.M_min);
  const double lnMmax = log(limits.M_max);
  
  return init_static_vars_only == 1 ? int_for_I11_y(lnMmin, (void*) ar): 
    Ntable.high_def_integration > 0 ?
    int_gsl_integrate_high_precision(int_for_I11_y, (void*) ar, lnMmin, lnMmax, 
      NULL, GSL_WORKSPACE_SIZE) :
    int_gsl_integrate_medium_precision(int_for_I11_y, (void*) ar, lnMmin, lnMmax, 
      NULL, GSL_WORKSPACE_SIZE);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_I02_yy(double logm, void* params)
{
  double* ar = (double*) params;
  const double a = ar[0];
  const double k1 = ar[1];
  const double k2 = ar[2];
  const double m = exp(logm);
  const double c = conc(m,a); 
  const double u = u_y_bnd(c, k1, m, a) * u_y_bnd(c, k2, m, a);
  const double tmp = m/(cosmology.rho_crit*cosmology.Omega_m);
  return massfunc(m, a) * m * u * tmp * tmp;
}

double I02_yy_nointerp(const double k1, const double k2, const double a, 
const int init_static_vars_only) 
{
  double ar[3] = {a, k1, k2};
  const double lnMmin = log(limits.M_min);
  const double lnMmax = log(limits.M_max);

  return init_static_vars_only == 1 ? int_for_I02_yy(lnMmin, (void*) ar) :
    Ntable.high_def_integration > 0 ?
    int_gsl_integrate_high_precision(int_for_I02_yy, (void*) ar, lnMmin, lnMmax, 
      NULL, GSL_WORKSPACE_SIZE) :
    int_gsl_integrate_medium_precision(int_for_I02_yy, (void*) ar, lnMmin, lnMmax, 
      NULL, GSL_WORKSPACE_SIZE);
}

double p_yy_1h_nointerp(const double k, const double a, const int init_static_vars_only) 
{
  const double I02YY = I02_yy_nointerp(k, k, a, init_static_vars_only);
  // convert to code unit, Table 2, 2009.01858
  const double ks = 0.05618/pow(cosmology.sigma_8*a, 1.013)*cosmology.coverH0; 
  const double x = ks*ks*ks*ks;
  return I02YY*(1.0/(x + 1.0)); // suppress lowk, Eq17, 2009.01858
}

double p_yy_2h_nointerp(const double k, const double a, const int init_static_vars_only)
{
  const double I11YY = I11_y_nointerp(k, a, init_static_vars_only);
  return I11YY*I11YY*p_lin(k,a);
}



// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_I02_my(double logm, void* params)
{
  double* ar = (double*) params;
  const double a = ar[0];
  const double k1 = ar[1];
  const double k2 = ar[2];
  const double m = exp(logm);
  const double c = conc(m, a);
  const double u = u_y_bnd(c, k1, m, a) * u_nfw_c(c, k2, m, a);
  const double tmp = m/(cosmology.rho_crit*cosmology.Omega_m);
  return massfunc(m, a) * m * u * tmp * tmp;
}

double I02_my_nointerp(const double k1, const double k2, const double a, 
const int init_static_vars_only) 
{
  double ar[3] = {a, k1, k2};
  const double lnMmin = log(limits.M_min);
  const double lnMmax = log(limits.M_max);

  return init_static_vars_only == 1 ? int_for_I02_my(lnMmin, (void*) ar) :
    Ntable.high_def_integration > 0 ?
    int_gsl_integrate_high_precision(int_for_I02_my, (void*) ar, lnMmin, lnMmax, 
      NULL, GSL_WORKSPACE_SIZE) :
    int_gsl_integrate_medium_precision(int_for_I02_my, (void*) ar, lnMmin, lnMmax, 
      NULL, GSL_WORKSPACE_SIZE);
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// HOD
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double n_c(const double mh, const double a, const int ni)
{
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  { // avoid segfault for accessing wrong array index
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }
  if (nuisance.hod[ni][0] < 10 || nuisance.hod[ni][0] > 16)
  {
    log_fatal("HOD parameters in redshift bin %d not set", ni);
    exit(1);
  }

  const double x = (log10(mh) - nuisance.hod[ni][0])/nuisance.hod[ni][1];
  gsl_sf_result ERF;
  {
    int status = gsl_sf_erf_e(x, &ERF);
    if (status) 
    {
      log_fatal(gsl_strerror(status));
      exit(1);
    }
  }
  
  return 0.5*(1.0 + ERF.val);
}

double n_s(const double mh, const double a, const int ni)
{
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  { // avoid segfault for accessing wrong array index
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }
  const double x = (mh - pow(10.0,nuisance.hod[ni][3]))/pow(10.0, nuisance.hod[ni][2]);
  const double ns = n_c(mh, a, ni)*pow(x, nuisance.hod[ni][4]);
  if (ns > 0)
  {
    return ns;
  }
  else 
  {
    return 1.e-15;
  }
}

double f_c(const int ni)
{
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  { // avoid segfault for accessing wrong array index
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }
  if (nuisance.hod[ni][5])
  {
    return nuisance.hod[ni][5];
  }
  return 1.0;
}

double u_g(const double k, const double m, const double a, const int ni)
{ // Fourier transformed normalized galaxy density profile, NFW with rescaled concentraion
  // analytic FT of NFW profile, from Cooray & Sheth 01
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  { // avoid segfault for accessing wrong array index
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }
  if (!(nuisance.gc[ni] > 0))
  {
    log_fatal("galaxy concentration parameter ill-defined in z-bin %d", ni); 
    exit(1);
  }
  return u_nfw_c(conc(m, a)*nuisance.gc[ni], k, m, a);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_ngal(double lnM, void *params)
{
  double* ar = (double*) params;
  
  const double a = ar[0];
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  const int ni = (int) ar[1];
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }

  const double m = exp(lnM);
  return massfunc(m,a)*m*(f_c(ni)*n_c(m, a, ni) + n_s(ni, a, ni));
}

double ngal_nointerp(
    const int ni, 
    const double a, 
    const int init_static_vars_only
  )
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
    const size_t szint = 1000 + 500 * (Ntable.high_def_integration);
    if (w != NULL)  gsl_integration_glfixed_table_free(w);
    w = malloc_gslint_glfixed(szint);
    cache[0] = Ntable.random;
  }

  double ar[2] = {a, (double) ni};
  const double lnMmin = log(10.0)*(nuisance.hod[ni][0] - 2.);
  const double lnMmax = log(limits.M_max);

  double res = 0.0;
  
  if (init_static_vars_only == 1)
    res = int_ngal((lnMmin + lnMmax)/2.0, (void*) ar);
  else
  {
    gsl_function F;
    F.params = (void*) ar;
    F.function = int_ngal;
    res = gsl_integration_glfixed(&F, lnMmin, lnMmax, w);
  }
  
  return res;
}

double ngal(const int ni, const double a)
{
  static double cache[MAX_SIZE_ARRAYS];
  static double** table = NULL;
  static double amin;
  static double amax;
  static double da;

  if (table == NULL || 
      fdiff(cache[1], Ntable.random) ||
      fdiff(cache[3], redshift.random_clustering)) 
  {  
    amin = 1.0/(redshift.clustering_zdist_zmax_all + 1.0);
    amax = 1.0/(redshift.clustering_zdist_zmin_all + 1.0);
    da = (amax - amin)/((double) Ntable.N_a - 1.0);

    if (table != NULL) free(table);
    table = (double**) malloc2d(redshift.clustering_nbin, Ntable.N_a);
  }
  
  if (fdiff(cache[0], cosmology.random) || 
      fdiff(cache[1], Ntable.random)    ||
      fdiff(cache[2], nuisance.random_galaxy_bias) ||
      fdiff(cache[3], redshift.random_clustering))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init = ngal_nointerp(0, amin, 1);
    }
    #pragma GCC diagnostic pop
    
    #pragma omp parallel for collapse(2)
    for (int i=0; i<redshift.clustering_nbin; i++) 
      for (int j=0; j<Ntable.N_a; j++) 
        table[j][i] = ngal_nointerp(i, amin + j*da, 0);
      
    cache[0] = cosmology.random;
    cache[1] = Ntable.random;
    cache[2] = nuisance.random_galaxy_bias;
    cache[3] = redshift.random_clustering;
  }

  if ((a < amin) || (a > amax))
    return 0.0;
  return interpol(table[ni], Ntable.N_a, amin, amax, da, a, 1.0, 1.0);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_m_mean(double lnM, void* params)
{
  double* ar  = (double*) params;
  
  const double a = ar[0];
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  const int ni = (int) ar[1];
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }

  const double m = exp(lnM);
  
  return massfunc(m, a)*m*m*(f_c(ni)*n_c(m, a, ni) + n_s(m, a, ni));
}

double m_mean_nointerp(
    const int ni, 
    const double a, 
    const int init_static_vars_only
  )
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
    const size_t szint = 1000 + 500 * (Ntable.high_def_integration);
    if (w != NULL)  gsl_integration_glfixed_table_free(w);
    w = malloc_gslint_glfixed(szint);
    cache[0] = Ntable.random;
  }

  double ar[2] = {a, (double) ni};

  const double lnMmin = log(10.0)*(nuisance.hod[ni][0] - 2.);
  const double lnMmax = log(limits.M_max);
  const double ngalaxies = ngal(ni, a);

  double res = 0.0;
  
  if (init_static_vars_only == 1)
    res = int_m_mean((lnMmin + lnMmax)/2.0, (void*) ar);
  else
  {
    gsl_function F;
    F.params = (void*) ar;
    F.function = int_m_mean;
    res = gsl_integration_glfixed(&F, lnMmin, lnMmax, w)/ngalaxies;
  }
  
  return res;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_bgal(double lnM, void* params)
{
  double* ar  = (double*) params;
  
  const double a = ar[0];
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  const int ni = (int) ar[1];
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }

  const double m = exp(lnM);
  return massfunc(m, a)*m*B1(m,a)*(f_c(ni)*n_c(m, a, ni) + n_s(m, a, ni));
}

double bgal_nointerp(
    const int ni, 
    const double a, 
    const int init_static_vars_only
  )
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
    const size_t szint = 1000 + 500 * (Ntable.high_def_integration);
    if (w != NULL)  gsl_integration_glfixed_table_free(w);
    w = malloc_gslint_glfixed(szint);
    cache[0] = Ntable.random;
  }

  double ar[2] = {a, (double) ni};
  const double lnMmin = log(10.0)*(nuisance.hod[ni][0] - 2.0);
  const double lnMmax = log(limits.M_max);

  double res = 0.0;
  
  if (init_static_vars_only == 1)
    res = int_bgal((lnMmin + lnMmax)/2.0, (void*) ar);
  else
  {
    gsl_function F;
    F.params = (void*) ar;
    F.function = int_bgal;
    res = gsl_integration_glfixed(&F, lnMmin, lnMmax, w);
  }
  
  return res;
}

double bgal(const int ni, const double a)
{
  static double cache[MAX_SIZE_ARRAYS];
  static double** table = NULL;
  static double amin;
  static double amax;
  static double da;

  if (table == NULL || 
      fdiff(cache[1], Ntable.random) ||
      fdiff(cache[3], redshift.random_clustering))  
  {  
    if (table != NULL) free(table); 
    table = (double**) malloc2d(redshift.clustering_nbin, Ntable.N_a);

    amin = 1.0/(redshift.clustering_zdist_zmax_all + 1.0);
    amax = 1.0/(redshift.clustering_zdist_zmin_all + 1.0);
    da = (amax - amin)/((double) Ntable.N_a - 1.0);
  }

  if (fdiff(cache[0], cosmology.random) || 
      fdiff(cache[1], Ntable.random)    ||
      fdiff(cache[2], nuisance.random_galaxy_bias) ||
      fdiff(cache[3], redshift.random_clustering))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init = bgal_nointerp(0, amin, 1);
    }
    #pragma GCC diagnostic pop
    
    #pragma omp parallel for collapse(2)
    for (int i=0; i<redshift.clustering_nbin; i++) 
      for (int j=0; j<Ntable.N_a; j++) 
        table[j][i] = bgal_nointerp(i, amin + j*da, 0);

    cache[0] = cosmology.random;
    cache[1] = Ntable.random;
    cache[2] = nuisance.random_galaxy_bias;
    cache[3] = redshift.random_clustering;
  }  

  if ((a < amin) || (a > amax))
    return 0.0;
  return interpol(table[ni], na, amin, amax, da, a, 1.0, 1.0);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_fsat(double lnM, void *params)
{
  double* ar  = (double*) params;
  
  double a = ar[0];
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  const int ni = (int) ar[1];
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }

  const double m = exp(lnM);
  return massfunc(m,a)*m*n_s(m, a, ni);
}

double fsat_nointerp(int ni, double a, const int init_static_vars_only)
{
  static double cache[MAX_SIZE_ARRAYS];
  static gsl_integration_glfixed_table* w = NULL;

  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  { // avoid segfault for accessing wrong array index
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }

  if (w == NULL || fdiff(cache[0], Ntable.random))
  {
    const size_t szint = 1000 + 500 * (Ntable.high_def_integration);
    if (w != NULL)  gsl_integration_glfixed_table_free(w);
    w = malloc_gslint_glfixed(szint);
    cache[0] = Ntable.random;
  }

  double ar[2] = {a, (double) ni};
  const double lnMmin = log(10.0)*(nuisance.hod[ni][0] - 2.);
  const double lnMmax = log(limits.M_max);
  const double ngalaxies = ngal(ni, a);

  if (init_static_vars_only == 1)
    res = int_fsat((lnMmin + lnMmax)/2.0, (void*) ar);
  else
  {
    gsl_function F;
    F.params = (void*) ar;
    F.function = int_fsat;
    res = gsl_integration_glfixed(&F, lnMmin, lnMmax, w)/ngalaxies;
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_G02(double logm, void* param)
{
  double* ar = (double*) param;
  
  const double k = ar[0];
  const double a = ar[1];
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  const int ni = (int) ar[2];
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }

  const double m = exp(logm);
  const double u = u_g(k, m, a, ni);
  const double ns = n_s(m, a, ni);
  
  return massfunc(m, a)*m*(u*u*ns*ns + 2.0*u*ns*n_c(m, a, ni)*f_c(ni));
}

double int_GM02(double logm, void* params)
{ // 1-halo galaxy-matter spectrum
  double* ar = (double*) params;
  
  const double k = ar[0];
  const double a = ar[1];
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  const int ni = (int) ar[2];
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }

  const double m = exp(logm);
  
  return (massfunc(m,a)*m*m/(cosmology.rho_crit*cosmology.Omega_m))*
          u_nfw_c(conc(m, a), k, m, a)*(u_g(k, m, a, ni)*n_s(m, a, ni) + 
          n_c(m, a, ni)*f_c(ni));
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double G02_nointerp(
    double k, 
    double a, 
    int ni, 
    const int init_static_vars_only
  )
{ //needs to be divided by ngal^2
  static double cache[MAX_SIZE_ARRAYS];
  static gsl_integration_glfixed_table* w = NULL;

  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  { // avoid segfault for accessing wrong array index
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }

  if (w == NULL || fdiff(cache[0], Ntable.random))
  {
    const size_t szint = 1000 + 500 * (Ntable.high_def_integration);
    if (w != NULL)  gsl_integration_glfixed_table_free(w);
    w = malloc_gslint_glfixed(szint);
    cache[0] = Ntable.random;
  }

  double ar[3] = {k, a, (double) ni};
  const double lnMmin = log(1.0e+8);
  const double lnMmax = log(limits.M_max);

  if (init_static_vars_only == 1)
    res = int_for_G02((lnMmin + lnMmax)/2.0, (void*) ar);
  else
  {
    gsl_function F;
    F.params = (void*) ar;
    F.function = int_for_G02;
    res = gsl_integration_glfixed(&F, lnMmin, lnMmax, w);
  }
}

double GM02_nointerp(
    double k, 
    double a, 
    int ni, 
    const int init_static_vars_only
  )
{ // needs to be divided by ngal
  static double cache[MAX_SIZE_ARRAYS];
  static gsl_integration_glfixed_table* w = NULL;

  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  { // avoid segfault for accessing wrong array index
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }

  if (w == NULL || fdiff(cache[0], Ntable.random))
  {
    const size_t szint = 1000 + 500 * (Ntable.high_def_integration);
    if (w != NULL)  gsl_integration_glfixed_table_free(w);
    w = malloc_gslint_glfixed(szint);
    cache[0] = Ntable.random;
  }

  double ar[3] = {k, a, (double) ni};
  const double lnMmin = log(10.)*(nuisance.hod[ni][0] - 1.0);
  const double lnMmax = log(limits.M_max);

  if (init_static_vars_only == 1)
    res = int_GM02((lnMmin + lnMmax)/2.0, (void*) ar);
  else
  {
    gsl_function F;
    F.params = (void*) ar;
    F.function = int_GM02;
    res = gsl_integration_glfixed(&F, lnMmin, lnMmax, w);
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/*
double int_G11(double logM, void* params)
{ // 2-halo term
  double* ar = (double*) params;
  
  const double k = ar[0];
  const double a = ar[1];
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  const int ni = (int) ar[2];
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }

  const double m = exp(logM);
  const double u = u_g(k, m, a, ni);
  
  return m*massfunc(m,a)*(u*n_s(m, a, ni) + n_c(m, a, ni)*f_c(ni))*B1(m, a);
}

double G11_nointerp(double k, double a, int ni, const int init_static_vars_only)
{ // needs to be devided by ngal
  double ar[3]  ={k, a, (double) ni};

  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  { // avoid segfault for accessing wrong array index
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }
  const double min = log(10.0)*(nuisance.hod[ni][0] - 1.0);
  const double max = log(limits.M_max);
  
  return (init_static_vars_only == 1) ? int_G11((min+max)/2.0, (void*) ar) :
    Ntable.high_def_integration > 0 ?
    int_gsl_integrate_high_precision(int_G11, (void*) ar, min, max, NULL, GSL_WORKSPACE_SIZE) : 
    int_gsl_integrate_medium_precision(int_G11, (void*) ar, min, max, NULL, GSL_WORKSPACE_SIZE);
}
*/

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// P_XY_NOINTERP (non-interpolated version of P_XY)
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double p_mm_1h_nointerp(
    const double k, 
    const double a, 
    const int init_static_vars_only
  ) 
{
  const double I02 = I02_nointerp(k, k, a, init_static_vars_only);
  return I02;
}

double p_mm_2h_nointerp(
    const double k, 
    const double a, 
    const int 
    init_static_vars_only
  ) 
{
  const double I11 = I11_nointerp(k, a, init_static_vars_only);
  return I11 * I11 * p_lin(k, a); 
}

double p_mm_nointerp(
    const double k, 
    const double a, 
    const int init_static_vars_only
  ) 
{ 
  const double PMM1H = p_mm_1h_nointerp(k, a, init_static_vars_only);
  const double PMM2H = p_mm_2h_nointerp(k, a, init_static_vars_only);
  return PMM1H + PMM2H; 
}

double p_gm_nointerp(
    const double k, 
    const double a, 
    const int ni, 
    const int use_2h_and_linear_only, 
    const int init_static_vars_only
  )
{
  const double bg = bgal(ni, a);
  
  if(use_2h_and_linear_only == 1)
    return p_lin(k, a)*bg;
  else
  {
    const double ng = ngal(ni, a);
    return Pdelta(k, a)*bg + GM02_nointerp(k, a, ni, init_static_vars_only)/ng;
  }
}

double p_gg_nointerp(
    const double k, 
    const double a, 
    const int ni, 
    const int nj,
    const int use_2h_and_linear_only, 
    const int init_static_vars_only
  )
{
  if (ni != nj)
  {
    log_fatal("cross-tomography (ni,nj) = (%d,%d) bins not supported", ni, nj);
    exit(1);
  }

  const double bg = bgal(ni, a);
  
  double res;
  
  if(use_2h_and_linear_only == 1)
    res = p_lin(k, a)*bg*bg;
  else
  {
    const double ng = ngal(ni, a);
    res = Pdelta(k, a)*bg*bg + 
            G02_nointerp(k, a, ni, init_static_vars_only)/(ng*ng);
  }

  return res;
}

double p_my_1h_nointerp(
    const double k, 
    const double a, 
    const int init_static_vars_only
  ) 
{
  const double I02MY = I02_my_nointerp(k, k, a, init_static_vars_only);
  
  // convert to code unit, Table 2, 2009.01858
  const double ks = 0.05618/pow(cosmology.sigma_8*a, 1.013)*cosmology.coverH0; 
  const double x = ks*ks*ks*ks;
  
  return I02MY*(1.0/(x + 1.0)); // suppress lowk, Eq17, 2009.01858
}

double p_my_2h_nointerp(
    const double k, 
    const double a, 
    const int init_static_vars_only
  )
{
  const double I11 = I11_nointerp(k, a, init_static_vars_only);
  
  const double I11YY = I11_y_nointerp(k, a, init_static_vars_only);
  
  return I11*I11YY*p_lin(k,a);
}

double p_my_nointerp(
    const double k, 
    const double a, 
    const int use_2h_only,
    const int init_static_vars_only
  )
{
  const double PM1H = 
    (use_2h_only == 1) ? 0.0 : p_my_1h_nointerp(k, a, init_static_vars_only);
  
  const double PM2H = p_my_2h_nointerp(k, a, init_static_vars_only);
  
  return PM1H + PM2H;
}

double p_yy_nointerp(
    const double k, 
    const double a, 
    const int use_2h_only, 
    const int init_static_vars_only
  )
{
  const double PYY1H = 
    (use_2h_only == 1) ? 0.0 : p_yy_1h_nointerp(k, a, init_static_vars_only);
  
  const double PYY2H = p_yy_2h_nointerp(k, a, init_static_vars_only);
  
  return PYY1H + PYY2H;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// P_XY (interpolated version of P_XY)
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double p_mm(const double k, const double a)
{ 
  static double cache[MAX_SIZE_ARRAYS];
  static double** table = NULL;

  const int na = Ntable.N_a;
  const double amin = limits.a_min;
  const double amax = 0.9999999;
  const double da = (amax - amin) / ((double) na - 1.0);
  
  const int nk = Ntable.N_k_nlin;
  const double lnkmin = log(limits.k_min_cH0);
  const double lnkmax = log(limits.k_max_cH0);
  const double dlnk = (lnkmax - lnkmin) / ((double) nk - 1.0);

  if (table == NULL || fdiff(cache[1], Ntable.random))
  {
    if (table != NULL) free(table);
    table = (double**) malloc2d(na, nk);
  }

  if (fdiff(cache[0], cosmology.random) || 
      fdiff(cache[1], Ntable.random))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      const int init_static_vars_only = 1; // TRUE
      double trash = p_mm_nointerp(exp(lnkmin), amin, init_static_vars_only);
    }
    #pragma GCC diagnostic pop
    
    #pragma omp parallel for collapse(2)
    for (int i=1; i<na; i++) 
      for (int j=0; j<nk; j++) 
        table[i][j] = log(p_mm_nointerp(exp(lnkmin + j*dk), amin + i*da, 0));

    cache[0] = cosmology.random;
    cache[1] = Ntable.random;
  }

/*
if (a < limits.a_min_hm){return Pdelta_halo(k,limits.a_min_hm)*pow(growfac(a)/growfac(limits.a_min_hm),2);}
  if (a > 0.999){return Pdelta_halo(k,0.999)*pow(growfac(a)/growfac(0.999),2);}
*/

  return exp(interpol2d(table[ni], na, amin, amax, da, a, Ntable.N_k_nlin, 
    lnkmin, lnkmax, dlnk, log(k), 0.0, 0.0));
}

double p_gm(
    const double k, 
    const double a, 
    const int ni, 
    const int use_2h_and_linear_only
  )
{
  static double cache[MAX_SIZE_ARRAYS];
  static double*** table = NULL;
  static double** alim = NULL; // range[0, :] = amin;
                               // range[1, :] = amax; 
                               // range[2, :] = da; 
  static double lnklim[3];  // lnklim[0] = lnkmin; 
                            // lnklim[1] = lnkmax; 
                            // lnklim[2] = dlnk; 

  const int na = (int) Ntable.N_a/5.0; // range is the (\delta a) of a single bin
  
  if (table == NULL || fdiff(cache[1], Ntable.random))
  {
    lnklim[0] = log(limits.k_min_cH0);
    lnklim[1] = log(limits.k_max_cH0);
    lnklim[2] = (lnklim[1] - lnklim[0])/((double) Ntable.N_k_nlin - 1.0);

    if (table != NULL) free(table);
    table = (double***) malloc3d(redshift.clustering_nbin, na, Ntable.N_k_nlin);

    if (alim != NULL) free(alim);
    alim = (double*) malloc2d(3, redshift.clustering_nbin);

    for (int l=0; l<redshift.clustering_nbin; l++) 
    {
      alim[0][l] = amin_lens(l);
      alim[1][l] = amax_lens(l);
      alim[2][l] = (alim[1][l] - alim[0][l])/((double) na - 1.0);
    }
  }

  if (fdiff(cache[0], cosmology.random) || 
      fdiff(cache[1], Ntable.random)    ||
      fdiff(cache[2], nuisance.random_galaxy_bias) ||
      fdiff(cache[3], redshift.random_clustering))
  { 
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      const int init_static_vars_only = 1; // TRUE
      double trash = p_gm_nointerp(exp(lnklim[0]), alim[0][0], 0, 
        use_2h_and_linear_only, init_static_vars_only);
    }
    #pragma GCC diagnostic pop
    
    #pragma omp parallel for collapse(3)
    for (int l=0; l<redshift.clustering_nbin; l++) 
      for (int i=0; i<na; i++) 
        for (int j=0; j<nk; j++) 
          table[l][i][j] = log(p_gm_nointerp(exp(lnklim[0] + j*lnklim[2]), 
            alim[0][l] + i*da[l], l, use_2h_and_linear_only, 0));

    cache[0] = cosmology.random;
    cache[1] = Ntable.random;
    cache[2] = nuisance.random_galaxy_bias;
    cache[3] = redshift.random_clustering;
  }

  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  { // avoid segfault for accessing wrong array index
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }

  return (a < alim[0][ni] || a > alim[1][ni]) ? 0.0 : exp(
      interpol2d(table[ni], na, alim[0][ni], alim[1][ni], alim[2][ni], 
      a, Ntable.N_k_nlin, lnklim[0], lnklim[1], lnklim[2], log(k), 0.0, 0.0));
}

double p_gg(
    const double k, 
    const double a, 
    const int ni, 
    const int nj, 
    const int use_2h_and_linear_only
  )
{
  static double cache[MAX_SIZE_ARRAYS];
  static double*** table = NULL;
  static double** alim = NULL; // range[0, :] = amin;
                               // range[1, :] = amax; 
                               // range[2, :] = da; 
  static double lnklim[3];  // lnklim[0] = lnkmin; 
                            // lnklim[1] = lnkmax; 
                            // lnklim[2] = dlnk; 
  
  const int na = (int) Ntable.N_a/5.0;

  if (table == NULL || fdiff(cache[1], Ntable.random))
  {
    lnklim[0] = log(limits.k_min_cH0);
    lnklim[1] = log(limits.k_max_cH0);
    lnklim[2] = (lnklim[1] - lnklim[0])/((double) Ntable.N_k_nlin - 1.0);

    if (table != NULL) free(table);
    table = (double***) malloc3d(redshift.clustering_nbin, na, Ntable.N_k_nlin);

    if (alim != NULL) free(alim);
    alim = (double*) malloc2d(3, redshift.clustering_nbin);

    for (int l=0; l<redshift.clustering_nbin; l++) 
    {
      alim[0][l] = amin_lens(l);
      alim[1][l] = amax_lens(l);
      alim[2][l] = (alim[1][l] - alim[0][l])/((double) na - 1.0);
    }
  }

  if (fdiff(cache[0], cosmology.random) || 
      fdiff(cache[1], Ntable.random)    ||
      fdiff(cache[2], nuisance.random_galaxy_bias) ||
      fdiff(cache[3], redshift.random_clustering))
  { 
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      const int init_static_vars_only = 1; // TRUE
      double trash = p_gg_nointerp(exp(lnklim[0]), alim[0][0], 0, 0, 
        use_2h_and_linear_only, init_static_vars_only);
    }
    #pragma GCC diagnostic pop
    
    #pragma omp parallel for collapse(3)
    for (int l=0; l<nbin; l++) 
      for (int i=0; i<na; i++) 
        for (int j=0; j<nk; j++) 
          table[l][i][j] = log(p_gg_nointerp(exp(lnklim[0] + j*lnklim[2]), 
            alim[0][l] + i*da[l], l, l, use_2h_and_linear_only, 0));
        
    cache[0] = cosmology.random;
    cache[1] = Ntable.random;
    cache[2] = nuisance.random_galaxy_bias;
    cache[3] = redshift.random_clustering;
  }
  
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }
  if (ni != nj)
  {
    log_fatal("cross-tomography (ni,nj) = (%d,%d) bins not supported", ni, nj);
    exit(1);
  }  

  return (a < alim[0][ni] || a > alim[1][ni]) ? 0.0 : exp(
      interpol2d(table[ni], na, alim[0][ni], alim[1][ni], alim[2][ni], 
      a, Ntable.N_k_nlin, lnklim[0], lnklim[1], lnklim[2], log(k), 0.0, 0.0));
}

double p_yy(
    const double k, 
    const double a, 
    const int use_2h_only
  )
{ 
  static double cache[MAX_SIZE_ARRAYS];
  static double** table = 0;

  const int na = Ntable.N_a;
  const double amin = limits.a_min;
  const double amax = 0.999999;
  const double da = (amax - amin) / ((double) na - 1.0);
  
  const int nk = Ntable.N_k_nlin;
  const double lnkmin = log(limits.k_min_cH0);
  const double lnkmax = log(limits.k_max_cH0);
  const double dlnk = (lnkmax - lnkmin) / ((double) nk - 1.0);

  if (table == NULL || fdiff(cache[1], Ntable.random))
  {
    if (table != NULL) free(table);
    table = (double**) malloc2d(na, nk);
  }

  if (fdiff(cache[0], cosmology.random) || 
      fdiff(cache[1], Ntable.random) ||
      fdiff(cache[2], Ntable.random_yparams))
  { 
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      const int init_static_vars_only = 1; // TRUE
      double init = 
        p_yy_nointerp(exp(lnkmin), amin, use_2h_only, init_static_vars_only);
    }
    #pragma GCC diagnostic pop   
    
    #pragma omp parallel for collapse(2)
    for (int i=0; i<na; i++) 
      for (int j=0; j<nk; j++) 
        table[i][j] = log(p_yy_nointerp(exp(lnkmin + j*dlnk), 
          amin + i*da, use_2h_only, 0));  
    
    cache[0] = cosmology.random;
    cache[1] = Ntable.random;
    cache[2] = nuisance.random_yparams;
  }

  return exp(interpol2d(table, na, amin, amax, da, a, nk, lnkmin, 
    lnkmax, dlnk, log(k), 1.0, 1.0));
}

double p_my(const double k, const double a, const int use_2h_only)
{ 
  static double cache[MAX_SIZE_ARRAYS];
  static double** table = 0;

  const int na = Ntable.N_a;
  const double amin = limits.a_min;
  const double amax = 0.9999999;
  const double da = (amax - amin) / ((double) na - 1.0);
  
  const int nk = Ntable.N_k_nlin;
  const double lnkmin = log(limits.k_min_cH0);
  const double lnkmax = log(limits.k_max_cH0);
  const double dlnk = (lnkmax - lnkmin) / ((double) nk - 1.0);

  if (table == 0)
  {
    table = (double**) malloc(sizeof(double*)*na);
    for(int i=0; i<na;i++)
    {
      table[i] = (double*) malloc(sizeof(double)*nk);
    } 
  }

  if (fdiff(cache[0], cosmology.random) || 
      fdiff(cache[1], Ntable.random) ||
      fdiff(cache[2], Ntable.random_yparams))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      const int init_static_vars_only = 1; // TRUE
      double init = 
        p_my_nointerp(exp(lnkmin), amin, use_2h_only, init_static_vars_only);
    }
    #pragma GCC diagnostic pop 
    
    #pragma omp parallel for collapse(2)
    for (int i=0; i<na; i++) 
      for (int j=0; j<nk; j++) 
        table[i][j] = log(p_my_nointerp(exp(lnkmin + j*dlnk), 
          amin + i*da, use_2h_only, 0));

    cache[0] = cosmology.random;
    cache[1] = Ntable.random;
    cache[2] = nuisance.random_yparams;
  }

  return exp(interpol2d(table, na, amin, amax, da, a, nk, lnkmin, 
    lnkmax, dlnk, log(k), 1.0, 1.0));
}

void set_HOD(int n)
{ 
  // n >= 0: set HOD parameters in redshift bin n; 
  // n = -1: unset HOD parameters (code then uses linear bias + non-linear matter power spectrum 
  // instead of halo model)

  const double z = zmean(n);
  const double a = 1.0/(z + 1.0);
  
  // set HOD parameters, parameterization of Zehavi et al.
  // these example values from Coupon et al. 2012 for red galaxies with M_r < -21.8 (Table B.2)
  // hod[zi][] = {lg(M_min), sigma_{lg M}, lg M_1, lg M_0, alpha, f_c}
  // gbias.cg[] = {f_g} (shift of concentration parameter: c_g(M) = f_g c(M))
  
  switch (n)
  {
    case -1:
    {
      for (int i = 0; i<MAX_SIZE_ARRAYS; i++)
      {
        nuisance.hod[i][0] = 0.;
        nuisance.hod[i][1] = 0.;
        nuisance.hod[i][2] = 0.;
        nuisance.hod[i][3] = 0.;
        nuisance.hod[i][4] = 0.;
        nuisance.hod[i][5] = 0.;
      }
      
      log_debug("unset HOD parameters; code uses non-linear matter power spectrum + linear bias");
      
      break;
    }
    case 0:
    {
      nuisance.hod[0][0] = 13.17;
      nuisance.hod[0][1] = 0.39;
      nuisance.hod[0][2] = 14.53;
      nuisance.hod[0][3] = 11.09;
      nuisance.hod[0][4] = 1.27;
      nuisance.hod[0][5] = 1.00;
      nuisance.gb[0][n] = bgal(n, a);
      
      log_debug("HOD quantities in bin %d at <z> = %.2f: <n_g> = %e(h/Mpc)^3", 
        n, z, ngal_nointerp(n, a, 0)*pow(cosmology.coverH0, -3.0));
      
      log_debug("HOD quantities in bin %d at <z> = %.2f: <b_g> = %.2f", 
        n, z, nuisance.gb[0][n]);
      
      log_debug("HOD quantities in bin %d at <z> = %.2f: lg(<M_h> h/M_sun) = %.3f", 
        n, z, log10(m_mean_nointerp(n, a, 0)));
      
      log_debug("HOD quantities in bin %d at <z> = %.2f: f_sat = %.3f", 
        n, z, fsat_nointerp(n, a, 0));
      
      break;
    }
    case 1:
    {
      nuisance.hod[1][0] = 13.18;
      nuisance.hod[1][1] = 0.30;
      nuisance.hod[1][2] = 14.47;
      nuisance.hod[1][3] = 10.93;
      nuisance.hod[1][4] = 1.36;
      nuisance.hod[1][5] = 1.00;
      
      nuisance.gb[0][n] = bgal(n, a);
      
      log_debug("HOD quantities in bin %d at <z> = %.2f: <n_g> = %e(h/Mpc)^3", n, z, 
        ngal_nointerp(n, a, 0)*pow(cosmology.coverH0, -3.0));
      
      log_debug("HOD quantities in bin %d at <z> = %.2f: <b_g> = %.2f", n, z, nuisance.gb[0][n]);
      
      log_debug("HOD quantities in bin %d at <z> = %.2f: lg(<M_h> h/M_sun) = %.3f", n, z, 
        log10(m_mean_nointerp(n, a, 0)));
      
      log_debug("HOD quantities in bin %d at <z> = %.2f: f_sat = %.3f",n, z, fsat_nointerp(n, a, 0));

      break;
    }
    case 2:
    {
      nuisance.hod[2][0] = 12.96;
      nuisance.hod[2][1] = 0.38;
      nuisance.hod[2][2] = 14.10;
      nuisance.hod[2][3] = 12.47;
      nuisance.hod[2][4] = 1.28;
      nuisance.hod[2][5] = 1.00;
      
      nuisance.gb[0][n] = bgal(n, a);

      log_debug("HOD quantities in bin %d at <z> = %.2f: <n_g> = %e(h/Mpc)^3", n, z, 
        ngal_nointerp(n, a, 0)*pow(cosmology.coverH0, -3.0));
      
      log_debug("HOD quantities in bin %d at <z> = %.2f: <b_g> = %.2f", n, z, nuisance.gb[0][n]);
      
      log_debug("HOD quantities in bin %d at <z> = %.2f: lg(<M_h> h/M_sun) = %.3f", n, z, 
        log10(m_mean_nointerp(n, a, 0)));
      
      log_debug("HOD quantities in bin %d at <z> = %.2f: f_sat = %.3f",n, z, fsat_nointerp(n, a, 0));

      break;
    }
    case 3:
    {
      nuisance.hod[3][0] = 12.80;
      nuisance.hod[3][1] = 0.35;
      nuisance.hod[3][2] = 13.94;
      nuisance.hod[3][3] = 12.15;
      nuisance.hod[3][4] = 1.52;
      nuisance.hod[3][5] = 1.00;
      
      nuisance.gb[0][n] = bgal(n, a);

      log_debug("HOD quantities in bin %d at <z> = %.2f: <n_g> = %e(h/Mpc)^3", n, z, 
        ngal_nointerp(n, a, 0)*pow(cosmology.coverH0, -3.0));
      
      log_debug("HOD quantities in bin %d at <z> = %.2f: <b_g> = %.2f", n, z, nuisance.gb[0][n]);
      
      log_debug("HOD quantities in bin %d at <z> = %.2f: lg(<M_h> h/M_sun) = %.3f", n, z, 
        log10(m_mean_nointerp(n, a, 0)));
      
      log_debug("HOD quantities in bin %d at <z> = %.2f: f_sat = %.3f",n, z, fsat_nointerp(n, a, 0));

      break;
    }
    case 4:
    {
      // no information for higher redshift populations - copy 1<z<1.2 values
      nuisance.hod[4][0] = 12.80;
      nuisance.hod[4][1] = 0.35;
      nuisance.hod[4][2] = 13.94;
      nuisance.hod[4][3] = 12.15;
      nuisance.hod[4][4] = 1.52;
      nuisance.hod[4][5] = 1.00;
      
      nuisance.gb[0][n] = bgal(n, a);

      log_debug("HOD quantities in bin %d at <z> = %.2f: <n_g> = %e(h/Mpc)^3", n, z, 
        ngal_nointerp(n, a, 0)*pow(cosmology.coverH0, -3.0));
      
      log_debug("HOD quantities in bin %d at <z> = %.2f: <b_g> = %.2f", n, z, nuisance.gb[0][n]);
      
      log_debug("HOD quantities in bin %d at <z> = %.2f: lg(<M_h> h/M_sun) = %.3f", n, z, 
        log10(m_mean_nointerp(n, a, 0)));
      
      log_debug("HOD quantities in bin %d at <z> = %.2f: f_sat = %.3f",n, z, fsat_nointerp(n, a, 0));

      break;
    }
    default:
    {
      printf("no HOD parameters specified to initialize bin %d\n", n);
    }
  }

  // for sparse lens populations, update angular galaxy density
  // TODO: WHAT?
  survey.n_lens = 2.0;

  log_debug("Using HOD model for lens population with n_lens = %.2f/arcmin^2", survey.n_lens);
}

