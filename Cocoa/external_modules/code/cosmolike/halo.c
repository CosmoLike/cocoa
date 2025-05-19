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

#include "log.c/src/log.h"

#define DEFAULT_INT_PREC 1000
#define delta_c 1.686
#define Delta 200

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// BASIC PEAK BACKGROUND SPLIT ROUTINES
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double hb1nu(const double nu, const double a)
{ // Halo bias based on peak-background split

  int ans;
  switch(like.halo_model[1])
  {
    case HALO_BIAS_TINKER_2010:
    {
      const double y = log10(200.0);
      
      const double ALPHA    = 1.0 + 0.24 * y * exp(-pow(4.0 / y, 4.0));
      const double nu_alpha = pow(nu, 0.44 * y - 0.88);
      
      const double BETA = 0.183;
      const double nu_beta = pow(nu, 1.5);
      
      const double GAMMA = 0.019 + 0.107 * y + 0.19 * exp(-pow(4.0 / y, 4.0));
      const double nu_gamma = pow(nu, 2.4);
      
      ans = 1.0 - ALPHA * nu_alpha / (nu_alpha + pow(delta_c, 0.44 * y - 0.88)) 
               + BETA * nu_beta + GAMMA * nu_gamma;
      break;
    }
    default:
    {
      log_fatal("like.halo_model[1] = %d not supported", like.halo_model[1]);
      exit(1);  
    }
  }
  return ans;
}

double fnu(const double nu, const double a)
{ // Halo bias based on peak-background split 
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }

  int ans;
  switch(like.halo_model[0])
  {
    case HMF_TINKER_2010:
    { // Eqs. (8-12) + Table 4 from Tinker et al. 2010
      const double aa = fmax(0.25, a); // limit fit range of mass function evolution to
                                       // z <= 3 (discussed after Eq. 12 of 1001.3162)
      const double alpha = 0.368;
      const double beta = 0.589 * pow(aa, -0.2);
      const double gamma = 0.864 * pow(aa, 0.01);
      const double phi = -0.729 * pow(aa, .08);
      const double eta = -0.243 * pow(aa, -0.27);

      ans = alpha*(1. + pow(beta*nu,-2*phi))*pow(nu,2*eta)*exp(-gamma*nu*nu/2.);
      break;
    }
    default:
    {
      log_fatal("like.halo_model[0] = %d not supported", like.halo_model[0]);
      exit(1);  
    }
  }
  return ans;
}

double conc(const double m, const double growfac_a) 
{
  int ans;
  switch(like.halo_model[2])
  {
    case CONCENTRATION_BHATTACHARYA_2013:
    { // Bhattacharya et al. 2013, Delta = 200 rho_{mean} (Table 2)
      const double nu = delta_c/(sqrt(sigma2(m))*growfac_a);
      ans =  9.0*pow(nu, -0.29)*pow(growfac_a, 1.15); 
      break;
    }
    default:
    {
      log_fatal("like.halo_model[2] = %d not supported", like.halo_model[2]);
      exit(1);  
    }
  }
  return ans;
}

double int_for_bias_norm(double nu, void* params) 
{ // correction for halo mass cuts so large-scale 2h matches PT at all redshifts 
  double* ar = (double*) params;
  const double a = ar[0];
  return hb1nu(nu, a) * fnu(nu, a);
}

double bias_norm_nointerp(
    const double a, 
    const int init
  )
{
  static double cache[MAX_SIZE_ARRAYS];
  static gsl_integration_glfixed_table* w = NULL;

  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }

  if (w == NULL || fdiff(cache[0], Ntable.random))
  {
    const size_t szint = DEFAULT_INT_PREC + 500*Ntable.high_def_integration;
    if (w != NULL)  gsl_integration_glfixed_table_free(w);
    w = malloc_gslint_glfixed(szint);
    cache[0] = Ntable.random;
  }

  const double growfac_a = growfac(a);
  const double nu_min = delta_c/(sqrt(sigma2(limits.M_min))*growfac_a);
  const double nu_max = delta_c/(sqrt(sigma2(limits.M_max))*growfac_a);

  double ar[2] = {a, growfac_a};

  double res;
  if (init == 1)
    res = int_for_bias_norm((nu_min+nu_max)/2.0, (void*) ar);
  else
  {
    gsl_function F;
    F.params = (void*) ar;
    F.function = int_for_bias_norm;
    res = gsl_integration_glfixed(&F, nu_min, nu_max, w);
  }
  return res;
}

double bias_norm(const double a) 
{
  static double cache[MAX_SIZE_ARRAYS];
  static double* table = NULL;
  static double lim[3];
  
  if (table == NULL || fdiff(cache[1], Ntable.random)) 
  {
    if (table != NULL) free(table);
    table = (double*) malloc(sizeof(double)*Ntable.N_a);
    lim[0] = limits.a_min; 
    lim[1] = 0.9999999;
    lim[2] = (lim[1] - lim[0]) / ((double) Ntable.N_a - 1.0);
  }
  if (fdiff(cache[0], cosmology.random) || 
      fdiff(cache[1], Ntable.random))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init = bias_norm_nointerp(lim[0], 1); 
    }
    #pragma GCC diagnostic pop

    #pragma omp parallel for
    for (int i=0; i<Ntable.N_a; i++) 
      table[i] = bias_norm_nointerp(lim[0] + i*lim[2], 0);
    table[Ntable.N_a-1] = 1.0;
    
    cache[0] = cosmology.random;
    cache[1] = Ntable.random;
  }
  return interpol1d(table, Ntable.N_a, lim[0], lim[1], lim[2], fmin(a,lim[1]-lim[2]));
}

double lognu0_gsl(double lnM, void* params __attribute__((unused))) 
{ 
  return log(delta_c/sqrt(sigma2(exp(lnM)))); 
}

double dlognudlogm(const double m) 
{ // if sigma(z) \propto to D(z), then d\ln \nu/dlnM independent of z
  static double cache[MAX_SIZE_ARRAYS];
  static double* table;
  static double lim[3];

  if (table == 0 || fdiff(cache[1], Ntable.random))
  {
    if (table != NULL) free(table);
    table = (double*) malloc(sizeof(double) * Ntable.N_M);
    lim[0] = log(limits.M_min);
    lim[1] = log(limits.M_max);
    lim[2] = (lim[1] - lim[0])/((double) Ntable.N_M - 1.0);
  }

  if (fdiff(cache[0], cosmology.random) || fdiff(cache[1], Ntable.random))
  {
    {
      const int i = 0;
      double result, abserr;
      double ar[1] = {1.0};

      gsl_function F;
      F.function = &lognu0_gsl;
      F.params = (void*) ar;

      int status = gsl_deriv_central(&F, lim[0] + i*lim[2], 
        0.1*(lim[0] + i*lim[2]), &result, &abserr);
      if (status) 
      {
        log_fatal(gsl_strerror(status));
        exit(1);
      }
      table[i] = result;
    }

    #pragma omp parallel for
    for (int i=1; i <Ntable.N_M; i++) 
    {
      double result, abserr;
      double ar[1] = {1.0};

      gsl_function F;
      F.function = &lognu0_gsl;
      F.params = (void*) ar;

      int status = gsl_deriv_central(&F, lim[0] + i*lim[2], 
        0.1*(lim[0] + i*lim[2]), &result, &abserr);
      if (status) 
      {
        log_fatal(gsl_strerror(status));
        exit(1);
      }
      table[i] = result;
    }

    cache[0] = cosmology.random;
    cache[1] = Ntable.random;
  }  
  return interpol1d(table, Ntable.N_M, lim[0], lim[1], lim[2], log(m));
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// HALO PROFILES
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double u_nfw_c(
    const double c, 
    const double k, 
    const double m, 
    const double a
  ) 
{ // analytic FT of NFW profile, from Cooray & Sheth 01 
  const double rho_delta = Delta * cosmology.rho_crit * cosmology.Omega_m;
  const double r_delta = pow(3./(4.0*M_PI)*(m/rho_delta), 1./3.);
  const double x = k * r_delta / c;
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

  return (sin(x)*(SI_XU.val - SI_X.val) 
          - sinl(c*x)/xu 
          + cos(x)*(CI_XU.val - CI_X.val))/(log(1. + c) - c/(1. + c));
}

double u_c(
    const double c, 
    const double k, 
    const double m, 
    const double a
  ) 
{
  int ans;
  switch(like.halo_model[3])
  {
    case HALO_PROFILE_NFW:
    {
      ans = u_nfw_c(c, k, m, a);
      break;
    }
    default:
    {
      log_fatal("like.halo_model[3] = %d not supported", like.halo_model[3]);
      exit(1);  
    }
  }
  return ans;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// GALAXY PROFILES
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


double u_g(
    const double c, 
    const double k, 
    const double m, 
    const double a,
    const int ni
  ) 
{
  return u_nfw_c(c*nuisance.gc[ni], k, m, a);

}


double HOD_nc(const double m, const double a, const int ni)
{
  if (!(a>0) || !(a<1)) 
  {
    log_fatal("a>0 and a<1 not true");
    exit(1);
  }
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  { 
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }
  if (nuisance.hod[ni][0] < 10 || nuisance.hod[ni][0] > 16)
  {
    log_fatal("HOD parameters in redshift bin %d not set", ni);
    exit(1);
  }

  const double x = (log10(m) - nuisance.hod[ni][0])/nuisance.hod[ni][1];
  
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

double HOD_ns(
    const double m, 
    const double a, 
    const int ni
  )
{
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  { 
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }
  
  const double x = (m - pow(10.0,nuisance.hod[ni][3]))/pow(10.0, nuisance.hod[ni][2]);
  const double ns = HOD_nc(m, a, ni)*pow(x, nuisance.hod[ni][4]);
  return (ns > 0) ? ns : 1.e-15;
}

double HOD_fc(const int ni)
{
  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }
  return (nuisance.hod[ni][5]) ? nuisance.hod[ni][5] : 1.0;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// GAS PROFILES
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_F0_KS(double x, void* params __attribute__((unused)))
{
  return x*x*pow(log(1.0 + x)/x, 1.0/(nuisance.gas[0] - 1.0));
}

double F0_KS_nointerp(double c, const int init)
{
  static double cache[MAX_SIZE_ARRAYS];
  static gsl_integration_glfixed_table* w = NULL;

  if (w == NULL || fdiff(cache[0], Ntable.random))
  {
    const size_t szint = DEFAULT_INT_PREC + 500*Ntable.high_def_integration;
    if (w != NULL)  gsl_integration_glfixed_table_free(w);
    w = malloc_gslint_glfixed(szint);
    cache[0] = Ntable.random;
  }

  double ar[1] = {0.0};
  const double xmin = 0.0;
  const double xmax = c;
  
  double res = 0.0;
  
  if (init == 1)
    res = int_F0_KS((xmin + xmax)/2.0, (void*) ar);
  else
  {
    gsl_function F;
    F.params = (void*) ar;
    F.function = int_F0_KS;
    res = gsl_integration_glfixed(&F, xmin, xmax, w);
  }

  return res;
}

double int_F_KS(double x, void* params)
{
  double* ar = (double*) params;
  const double y = ar[0];
  
  return (x*sinl(y*x)/y)*pow(log(1.0 + x)/x, 
      nuisance.gas[0]/(nuisance.gas[0] - 1.0));
}

double F_KS_nointerp(double c, double krs, const int init) 
{
  static double cache[MAX_SIZE_ARRAYS];
  static gsl_integration_glfixed_table* w = NULL;

  if (w == NULL || fdiff(cache[0], Ntable.random))
  {
    const size_t szint = DEFAULT_INT_PREC + 500*Ntable.high_def_integration;
    if (w != NULL)  gsl_integration_glfixed_table_free(w);
    w = malloc_gslint_glfixed(szint);
    cache[0] = Ntable.random;
  }

  double ar[1] = {krs};
  const double cmin = 0.0;
  const double cmax = c;

  double res = 0.0;
  
  if (init == 1)
    res = int_F_KS((cmin + cmax)/2.0, (void*) ar);
  else
  {
    gsl_function F;
    F.params = (void*) ar;
    F.function = int_F_KS;
    res = gsl_integration_glfixed(&F, cmin, cmax, w);
  }

  return res;
}

double u_KS(double c, double k, const double rv)
{
  static double cache[MAX_SIZE_ARRAYS];
  static double** table = 0;
  static double* norm = 0;
  static double  lim[2][3]; // lim[0][0] = cmin;  lim[1][0] = lnxmin;
                            // lim[0][1] = cmax;  lim[1][1] = lnxmax;
                            // lim[0][2] = dc;    lim[1][2] = dlnx; 
 
  if (table == NULL || fdiff(cache[1], Ntable.random))
  {   
    if (table != NULL) free(table); 
    table = (double**) malloc2d(Ntable.halo_uKS_nc, Ntable.halo_uks_nx);
    if (norm != NULL) free(norm); 
    norm = (double*) malloc(sizeof(double)*Ntable.halo_uKS_nc);

    lim[0][0] = limits.halo_uKS_cmin; 
    lim[0][1] = limits.halo_uKS_cmax;
    lim[0][2] = (lim[0][1] - lim[0][0])/((double) Ntable.halo_uKS_nc - 1.0);
    lim[1][0] = log(limits.halo_uKS_xmin); // full range of possible k*R_200/c in the code
    lim[1][1] = log(limits.halo_uKS_xmax); 
    lim[1][2] = (lim[1][1] - lim[1][0])/((double) Ntable.halo_uks_nx - 1.0); 
  }

  if (fdiff(cache[0], nuisance.random_gas) || fdiff(cache[1], Ntable.random))
  { 
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init  = F0_KS_nointerp(lim[0][0], 1);
      init = F_KS_nointerp(lim[0][0], exp(lim[1][0]), 1);
    }
    #pragma GCC diagnostic pop
    
    #pragma omp parallel for
    for (int i=0; i<Ntable.halo_uKS_nc; i++) 
      norm[i] = F0_KS_nointerp(lim[0][0] + i*lim[0][2], 0);

    #pragma omp parallel for collapse(2)
    for (int i=0; i<Ntable.halo_uKS_nc; i++)
      for (int j=0; j<Ntable.halo_uks_nx; j++) 
        table[i][j] = F_KS_nointerp(lim[0][0] + i*lim[0][2], 
                                    exp(lim[1][0] + j*lim[1][2]), 0)/norm[i];

    cache[0] = nuisance.random_gas; 
    cache[1] = Ntable.random;
  }
  return interpol2d(table, 
    Ntable.halo_uKS_nc, lim[0][0], lim[0][1], lim[0][2], c, 
    Ntable.halo_uks_nx, lim[1][0], lim[1][1], lim[1][2], log(k * rv/c));
}

double frac_bnd(double M)
{
  const double M0 = pow(10.0, nuisance.gas[2]);
  return cosmology.Omega_b/(cosmology.Omega_m*(1.0+ pow(M0/M, nuisance.gas[1])));
}

double frac_ejc(double M)
{
  const double logM = log10(M);
  const double delta = (logM - nuisance.gas[7])/nuisance.gas[8];
  
  const double tmp = nuisance.gas[6] * exp(-0.5*delta*delta);  
  const double frac_star = ((logM > nuisance.gas[2]) && 
                           (tmp < nuisance.gas[6]/3.0)) ? nuisance.gas[6]/3.0 : tmp; 
  
  return frac_bnd(M) - frac_star;
}

double u_y_bnd(double c, double k, double m, double a)
{ //unit: [G(M_solar/h)^2 / (c/H0)]
  
  const double rho_delta = Delta * cosmology.rho_crit * cosmology.Omega_m;
  const double r_delta = pow(3./(4.0*M_PI)*(m/rho_delta), 1./3.);
  const double rv = r_delta;

  const double mu_p = 4.0/(3.0 + 5*nuisance.gas[10]);
  const double mu_e = 2.0/(1.0 + nuisance.gas[10]);
  
  return (2.0*nuisance.gas[5]/(3.0*a))*(mu_p/mu_e)*frac_bnd(m)*m*(m/rv)*u_KS(c, k, rv);
}

double u_y_ejc(double m)
{ // [m] = [Msun/h]
  const double num_p = 1.1892e57; // proton number in 1 solar mass, unit [1/Msun]
  
  // convert ejected gas T (K) to E (eV) then to [G (Msun/h)^2 / (c/H0) * h]
  const double E_w = pow(10,nuisance.gas[9]) * 8.6173e-5 * 5.616e-44;
  const double mu_e = 2./(1.+nuisance.gas[10]);
  
  return (num_p * m * frac_ejc(m) / mu_e) * E_w; // final unit in [G(Msun/h)^2 / (c/H0)]
}

double n_s_cmv(double a) 
{ 
  double dV_dz = pow(f_K(chi(a)), 2.0) / hoverh0(a); // comoving dV/dz per radian^2
  return zdistr_photoz(1.0/a - 1., -1) * survey.n_gal * 
    survey.n_gal_conversion_factor / dV_dz; // dN/dz/radian^2/(dV/dz/radian^2)
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// HALO MODEL ROUTINES
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_hm_funcs(double lnM, void* params)
{ // 0 = ngal, 1 = m_mean, 2 = fsat, 3 = bgal 
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
  const int func = (int) ar[2];
  const double growfac_a = (double) ar[3];

  const double m = exp(lnM);
  
  const double nu = delta_c/(sqrt(sigma2(m))*growfac_a);
  const double gnu  = fnu(nu, a) * nu; 
  const double rhom = cosmology.rho_crit * cosmology.Omega_m;
  const double dNdlnM = gnu * (rhom/m) * dlognudlogm(m);
  
  const double nc = HOD_fc(ni)*HOD_nc(m, a, ni);
  const double ns = HOD_ns(ni, a, ni);
  
  double res;
  switch(func)
  {
    case 0:
    { // N_gal = \int dM n(M)*(nc + ns) = \int dlnM M n(M)*(nc + ns) 
      res = dNdlnM*(nc + ns);
      break;
    }
    case 1:
    { // <M> = \int dM M*n(M)*(nc + ns) = \int dlnM M^2 n(M)*(nc + ns) 
      res = m*(dNdlnM*(nc + ns));
      break;
    }
    case 2:
    {
      res = dNdlnM*ns;
      break;
    }
    case 3:
    {
      res = hb1nu(nu, a)*(dNdlnM*(nc + ns));
      break;
    }
    default:
    {
      log_fatal("option not supported");
      exit(1);
    }
  }
  return res;
}

double ngal_nointerp(
    const int ni, 
    const double a, 
    const int init
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
    const size_t szint = DEFAULT_INT_PREC + 500*Ntable.high_def_integration;
    if (w != NULL)  gsl_integration_glfixed_table_free(w);
    w = malloc_gslint_glfixed(szint);
    cache[0] = Ntable.random;
  }

  double ar[4] = {a, (double) ni, (double) 0, growfac(a)};
  const double lnMmin = log(10.0)*(nuisance.hod[ni][0] - 2.);
  const double lnMmax = log(limits.M_max);

  double res = 0.0;
  if (init == 1)
    res = int_hm_funcs((lnMmin + lnMmax)/2.0, (void*) ar);
  else
  {
    gsl_function F;
    F.params = (void*) ar;
    F.function = int_hm_funcs;
    res = gsl_integration_glfixed(&F, lnMmin, lnMmax, w);
  }
  return res;
}

double ngal(const int ni, const double a)
{
  static double cache[MAX_SIZE_ARRAYS];
  static double** table = NULL;
  static double lim[3]; // [0] = amin; [1] = amax; [2] = da

  if (table == NULL || 
      fdiff(cache[1], Ntable.random) ||
      fdiff(cache[3], redshift.random_clustering)) 
  { 
    if (table != NULL) free(table);
    table = (double**) malloc2d(redshift.clustering_nbin, Ntable.N_a);

    lim[0] = 1.0/(redshift.clustering_zdist_zmax_all + 1.0);
    lim[1] = 1.0/(redshift.clustering_zdist_zmin_all + 1.0);
    lim[2] = (lim[1] - lim[0])/((double) Ntable.N_a - 1.0);
  }
  
  if (fdiff(cache[0], cosmology.random) || 
      fdiff(cache[1], Ntable.random)    ||
      fdiff(cache[2], nuisance.random_galaxy_bias) ||
      fdiff(cache[3], redshift.random_clustering))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init = ngal_nointerp(0, lim[0], 1);
    }
    #pragma GCC diagnostic pop
    
    #pragma omp parallel for collapse(2)
    for (int i=0; i<redshift.clustering_nbin; i++) 
      for (int j=0; j<Ntable.N_a; j++) 
        table[j][i] = ngal_nointerp(i, lim[0] + j*lim[2], 0);

    cache[0] = cosmology.random;
    cache[1] = Ntable.random;
    cache[2] = nuisance.random_galaxy_bias;
    cache[3] = redshift.random_clustering;
  }

  if ((a < lim[0]) || (a > lim[1]))
    return 0.0;
  return interpol1d(table[ni], Ntable.N_a, lim[0], lim[1], lim[2], a);
}

double hm_funcs_nointerp(
    const int ni, 
    const double a, 
    const int func,
    const int init
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
    const size_t szint = DEFAULT_INT_PREC + 500*Ntable.high_def_integration;
    if (w != NULL)  gsl_integration_glfixed_table_free(w);
    w = malloc_gslint_glfixed(szint);
    cache[0] = Ntable.random;
  }

  double ar[4] = {a, (double) ni, (double) func, growfac(a)}; 
  const double lnMmin = log(10.0)*(nuisance.hod[ni][0] - 2.);
  const double lnMmax = log(limits.M_max);

  double res = 0.0;
  if (init == 1)
    res = int_hm_funcs((lnMmin + lnMmax)/2.0, (void*) ar);
  else
  {
    gsl_function F;
    F.params = (void*) ar;
    F.function = int_hm_funcs;
    res = gsl_integration_glfixed(&F, lnMmin, lnMmax, w);
  }
  return (func == 1 || func == 2) ? res/ngal(ni, a) : res;
}

double mmean_nointerp(
    const int ni, 
    const double a, 
    const int init
  )
{
  return hm_funcs_nointerp(ni, a, 1, init);
}

double fsat_nointerp(
    const int ni, 
    const double a, 
    const int init
  )
{
  return hm_funcs_nointerp(ni, a, 2, init);
} 

double bgal_nointerp(
    const int ni, 
    const double a, 
    const int init
  )
{
  return hm_funcs_nointerp(ni, a, 3, init);
}

double bgal(const int ni, const double a)
{
  static double cache[MAX_SIZE_ARRAYS];
  static double** table = NULL;
  static double lim[3]; // [0] = amin; [1] = amax; [2] = da

  if (table == NULL || 
      fdiff(cache[1], Ntable.random) ||
      fdiff(cache[3], redshift.random_clustering))  
  {  
    if (table != NULL) free(table); 
    table = (double**) malloc2d(redshift.clustering_nbin, Ntable.N_a);
    
    lim[0] = 1.0/(redshift.clustering_zdist_zmax_all + 1.0);
    lim[1] = 1.0/(redshift.clustering_zdist_zmin_all + 1.0);
    lim[2] = (lim[1] - lim[0])/((double) Ntable.N_a - 1.0);
  }
  if (fdiff(cache[0], cosmology.random) || 
      fdiff(cache[1], Ntable.random)    ||
      fdiff(cache[2], nuisance.random_galaxy_bias) ||
      fdiff(cache[3], redshift.random_clustering))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init = bgal_nointerp(0, lim[0], 1);
    }
    #pragma GCC diagnostic pop
    
    #pragma omp parallel for collapse(2)
    for (int i=0; i<redshift.clustering_nbin; i++) 
      for (int j=0; j<Ntable.N_a; j++) 
        table[j][i] = bgal_nointerp(i, lim[0] + j*lim[2], 0);

    cache[0] = cosmology.random;
    cache[1] = Ntable.random;
    cache[2] = nuisance.random_galaxy_bias;
    cache[3] = redshift.random_clustering;
  }  
  return (a < lim[0]) || (a > lim[1]) ? 0.0 : 
    interpol1d(table[ni], Ntable.N_a, lim[0], lim[1], lim[2], a);
}

double int_for_I02_XY(double lnM, void* params) 
{
  double* ar = (double*) params;
  const double a = ar[0];
  const double k1 = ar[1];
  const double k2 = ar[2];
  const int XY = (int) ar[3];
  const double growfac_a = ar[4];

  const double m = exp(lnM);
  
  const double nu = delta_c/(sqrt(sigma2(m))*growfac_a);
  const double gnu  = fnu(nu, a) * nu; 
  const double rhom = cosmology.rho_crit * cosmology.Omega_m;
  const double dNdlnM = gnu * (rhom/m) * dlognudlogm(m); // mass function

  const double c = conc(m, growfac_a);

  double u;
  switch(XY)
  {
    case 0:
    { // matter-matter
      u = u_c(c, k1, m, a) * u_c(c, k2, m, a);
      break;
    }
    case 1:
    { // matter-y
      u = u_y_bnd(c, k1, m, a) * u_c(c, k2, m, a);
      break;
    }
    case 2:
    { // y-y 
      u = u_y_bnd(c, k1, m, a) * u_y_bnd(c, k2, m, a);
      break;
    }
    default:
    {
      log_fatal("option not supported");
      exit(1);
    }
  }
  return dNdlnM * u * (m/rhom) * (m/rhom);
}  

double I02_XY_nointerp(
    const double k1, 
    const double k2, 
    const double a,
    const int func, // 1 = MM, 2 = MY, 3 = YY
    const int init
  ) 
{
  static double cache[MAX_SIZE_ARRAYS];
  static gsl_integration_glfixed_table* w = NULL;

  if (w == NULL || fdiff(cache[0], Ntable.random))
  {
    const size_t szint = DEFAULT_INT_PREC + 500*Ntable.high_def_integration;
    if (w != NULL)  gsl_integration_glfixed_table_free(w);
    w = malloc_gslint_glfixed(szint);
    cache[0] = Ntable.random;
  }

  double ar[5] = {a, k1, k2, func, growfac(a)};
  const double lnMmin = log(limits.M_min);
  const double lnMmax = log(limits.M_max);

  double res;
  if (init == 1)
    res = int_for_I02_XY((lnMmin + lnMmax)/2.0, (void*) ar);
  else
  {
    gsl_function F;
    F.params = (void*) ar;
    F.function = int_for_I02_XY;
    res = gsl_integration_glfixed(&F, lnMmin, lnMmax, w);
  }
  return res;
}

double int_for_I11_X(double lnM, void* params) 
{
  const double* ar = (double*) params;  
  const double a = ar[0];
  const double k = ar[1];
  const int func = (int) ar[2];
  const double growfac_a = ar[3];

  const double m = exp(lnM);
  
  const double nu = delta_c/(sqrt(sigma2(m))*growfac_a);
  const double gnu = fnu(nu, a) * nu; 
  const double rhom = cosmology.rho_crit * cosmology.Omega_m;
  const double dNdlnM = gnu * (rhom/m) * dlognudlogm(m);

  const double c = conc(m, growfac_a);

  double u;
  switch(func)
  {
    case 0:
    { // matter
      u = u_c(c, k, m, a);
      break;
    }
    case 1:
    { // y
      u = u_y_bnd(c, k, m, a) + u_y_ejc(m);
      break;
    }
    default:
    {
      log_fatal("option not supported");
      exit(1);
    }
  }
  return dNdlnM * u * (m/rhom) * hb1nu(nu, a)/bias_norm(a);
}

double I11_X_nointerp(
    const double k, 
    const double a,
    const int func, 
    const int init
  ) 
{ 
  static double cache[MAX_SIZE_ARRAYS];
  static gsl_integration_glfixed_table* w = NULL;

  if (w == NULL || fdiff(cache[0], Ntable.random))
  {
    const size_t szint = DEFAULT_INT_PREC + 500*Ntable.high_def_integration;
    if (w != NULL)  gsl_integration_glfixed_table_free(w);
    w = malloc_gslint_glfixed(szint);
    cache[0] = Ntable.random;
  }

  double ar[4] = {a, k, func, growfac(a)};
  const double lnMmin = log(limits.M_min);
  const double lnMmax = log(limits.M_max);
  
  double res;
  if (init == 1)
    res = int_for_I11_X((lnMmin + lnMmax)/2.0, (void*) ar);
  else
  {
    gsl_function F;
    F.params = (void*) ar;
    F.function = int_for_I11_X;
    res = gsl_integration_glfixed(&F, lnMmin, lnMmax, w);
  }
  return res;
}

double int_for_G02(double lnM, void* param)
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
  const double growfac_a = ar[3];

  const double m  = exp(lnM);

  const double nu = delta_c/(sqrt(sigma2(m))*growfac_a);
  const double gnu = fnu(nu, a) * nu; 
  const double rhom = cosmology.rho_crit * cosmology.Omega_m;
  const double dNdlnM = gnu * (rhom/m) * dlognudlogm(m);

  const double c  = conc(m, growfac_a);
  const double u  = u_g(c, k, m, a, ni);
  const double ns = HOD_ns(m, a, ni);
  const double nc = HOD_nc(m, a, ni);
  const double fc = HOD_fc(ni);



  return dNdlnM*(u*u*ns*ns + 2.0*u*ns*nc*fc);
}

double G02_nointerp(
    double k, 
    double a, 
    int ni, 
    const int init
  )
{ //needs to be divided by ngal^2
  static double cache[MAX_SIZE_ARRAYS];
  static gsl_integration_glfixed_table* w = NULL;

  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }

  if (w == NULL || fdiff(cache[0], Ntable.random))
  {
    const size_t szint = DEFAULT_INT_PREC + 500*Ntable.high_def_integration;
    if (w != NULL)  gsl_integration_glfixed_table_free(w);
    w = malloc_gslint_glfixed(szint);
    cache[0] = Ntable.random;
  }

  double ar[4] = {k, a, (double) ni, growfac(a)};
  const double lnMmin = log(1.0e+8);
  const double lnMmax = log(limits.M_max);

  double res;
  if (init == 1)
    res = int_for_G02((lnMmin + lnMmax)/2.0, (void*) ar);
  else
  {
    gsl_function F;
    F.params = (void*) ar;
    F.function = int_for_G02;
    res = gsl_integration_glfixed(&F, lnMmin, lnMmax, w);
  }
  return res;
}

double int_GM02(double lnM, void* params)
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
  const double growfac_a = ar[3];

  const double m = exp(lnM);

  const double nu = delta_c/(sqrt(sigma2(m))*growfac_a);
  const double gnu  = fnu(nu, a) * nu; 
  const double rhom = cosmology.rho_crit * cosmology.Omega_m;
  const double dNdlnM = gnu * (rhom/m) * dlognudlogm(m);

  const double c = conc(m, growfac_a);
  const double ns = HOD_ns(m, a, ni);
  const double nc = HOD_nc(m, a, ni);
  const double fc = HOD_fc(ni);

  return dNdlnM*(m/rhom)*u_c(c,k,m,a)*(u_g(c,k,m,a,ni)*ns + nc*fc);
}

double GM02_nointerp(
    double k, 
    double a, 
    int ni, 
    const int init
  )
{ // needs to be divided by ngal
  static double cache[MAX_SIZE_ARRAYS];
  static gsl_integration_glfixed_table* w = NULL;

  if (ni < 0 || ni > redshift.clustering_nbin - 1)
  {
    log_fatal("error in selecting bin number ni = %d", ni);
    exit(1);
  }

  if (w == NULL || fdiff(cache[0], Ntable.random))
  {
    const size_t szint = DEFAULT_INT_PREC + 500*Ntable.high_def_integration;
    if (w != NULL)  gsl_integration_glfixed_table_free(w);
    w = malloc_gslint_glfixed(szint);
    cache[0] = Ntable.random;
  }

  double ar[4] = {k, a, (double) ni, growfac(a)};
  const double lnMmin = log(10.)*(nuisance.hod[ni][0] - 1.0);
  const double lnMmax = log(limits.M_max);

  double res;
  if (init == 1)
    res = int_GM02((lnMmin + lnMmax)/2.0, (void*) ar);
  else
  {
    gsl_function F;
    F.params = (void*) ar;
    F.function = int_GM02;
    res = gsl_integration_glfixed(&F, lnMmin, lnMmax, w);
  }
  return res;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// HALO MODEL POWER SPECTRA
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double p_xy_nointerp(
    const double k, 
    const double a,
    const int func,
    const int init
  ) 
{
  const double I02 = I02_XY_nointerp(k, k, a, func, init);

  double P1H, I11X, I11Y;

  switch(func)
  {
    case 0:
    { // PMM
      P1H  = I02;
      I11X = I11_X_nointerp(k, a, func, init);
      I11Y = I11X;
      break;
    }
    case 1:
    { // PMY
      // convert to code unit, Table 2, 2009.01858
      const double ks = 0.05618/pow(cosmology.sigma_8*a, 1.013)*cosmology.coverH0; 
      const double x = ks*ks*ks*ks;

      P1H  = I02*(1.0/(x + 1.0)); // suppress lowk (Eq17;2009.01858)
      I11X = I11_X_nointerp(k, a, 0, init);
      I11Y = I11_X_nointerp(k, a, 2, init);
      break;
    }
    case 2:
    { // PYY
      // convert to code unit, Table 2, 2009.01858
      const double ks = 0.05618/pow(cosmology.sigma_8*a, 1.013)*cosmology.coverH0; 
      const double x = ks*ks*ks*ks;
      
      P1H  = I02*(1.0/(x + 1.0)); // suppress lowk (Eq17;2009.01858)
      I11X = I11_X_nointerp(k, a, func, init);
      I11Y = I11X;
      break;
    }
    default:
    {
      log_fatal("option not supported");
      exit(1);
    }
  }

  return P1H + (I11X * I11Y * p_lin(k, a));;
}

double p_mm(
    const double k, 
    const double a
  )
{ 
  static double cache[MAX_SIZE_ARRAYS];
  static double** table = NULL;
  static double lim[2][3]; // lim[0][0] = amin, lim[0][1] = amax, lim[0][2] = da 
                           // lim[1][0] = lnkmin, lim[1][1] = lnkmax, lim[1][2] = dlnk

  if (table == NULL || fdiff(cache[1], Ntable.random))
  {
    if (table != NULL) free(table);
    table = (double**) malloc2d(Ntable.N_a, Ntable.N_k_nlin);
    
    lim[0][0] = limits.a_min;
    lim[0][1] = 0.9999999;
    lim[0][2] = (lim[0][1] - lim[0][0]) / ((double) Ntable.N_a - 1.0);
    lim[1][0] = log(limits.k_min_cH0);
    lim[1][1] = log(limits.k_max_cH0);
    lim[1][2] = (lim[1][1] - lim[1][0]) / ((double) Ntable.N_k_nlin - 1.0);
  }

  if (fdiff(cache[0], cosmology.random) || fdiff(cache[1], Ntable.random))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    { // p_xy_nointerp(k, a, func, init)
      double init = p_xy_nointerp(exp(lim[1][0]), lim[0][0], 0, 1);
    }
    #pragma GCC diagnostic pop
    
    #pragma omp parallel for collapse(2)
    for (int i=1; i<Ntable.N_a; i++) 
      for (int j=0; j<Ntable.N_k_nlin; j++) 
        table[i][j] = log(p_xy_nointerp(exp(lim[1][0] + j*lim[1][2]), 
                                        lim[0][0] + i*lim[0][2], 0, 0));
    cache[0] = cosmology.random;
    cache[1] = Ntable.random;
  }

/*
if (a < limits.a_min_hm){return Pdelta_halo(k,limits.a_min_hm)*pow(growfac(a)/growfac(limits.a_min_hm),2);}
  if (a > 0.999){return Pdelta_halo(k,0.999)*pow(growfac(a)/growfac(0.999),2);}
*/
  return exp(interpol2d(table, 
    Ntable.N_a, lim[0][0], lim[0][1], lim[0][2], a, 
    Ntable.N_k_nlin, lim[1][0], lim[1][1], lim[1][2], log(k)));
}

double p_my(
    const double k, 
    const double a
  )
{ 
  static double cache[MAX_SIZE_ARRAYS];
  static double** table = 0;
  static double lim[2][3]; // lim[0][0] = amin, lim[0][1] = amax, lim[0][2] = da 
                           // lim[1][0] = lnkmin, lim[1][1] = lnkmax, lim[1][2] = dlnk

  if (table == NULL || fdiff(cache[1], Ntable.random))
  {
    if (table != NULL) free(table);
    table = (double**) malloc2d(Ntable.N_a, Ntable.N_k_nlin); 
    
    lim[0][0] = limits.a_min;
    lim[0][1] = 0.9999999;
    lim[0][2] = (lim[0][1] - lim[0][0]) / ((double) Ntable.N_a - 1.0);
    lim[1][0] = log(limits.k_min_cH0);
    lim[1][1] = log(limits.k_max_cH0);
    lim[1][2] = (lim[1][1] - lim[1][0]) / ((double) Ntable.N_k_nlin - 1.0);
  }
  if (fdiff(cache[0], cosmology.random) || 
      fdiff(cache[1], Ntable.random) ||
      fdiff(cache[2], nuisance.random_gas))
  {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init = p_xy_nointerp(exp(lim[1][0]), lim[0][0], 1, 1);
    }
    #pragma GCC diagnostic pop 
    
    #pragma omp parallel for collapse(2)
    for (int i=1; i<Ntable.N_a; i++) 
      for (int j=0; j<Ntable.N_k_nlin; j++) 
        table[i][j] = log(p_xy_nointerp(exp(lim[1][0] + j*lim[1][2]), 
                                            lim[0][0] + i*lim[0][2], 1, 0));

    cache[0] = cosmology.random;
    cache[1] = Ntable.random;
    cache[2] = nuisance.random_gas;
  }
  return exp(interpol2d(table, 
    Ntable.N_a, lim[0][0], lim[0][1], lim[0][2], a, 
    Ntable.N_k_nlin, lim[1][0], lim[1][1], lim[1][2], log(k)));
}

double p_yy(
    const double k, 
    const double a
  )
{ 
  static double cache[MAX_SIZE_ARRAYS];
  static double** table = 0;
  static double lim[2][3]; // lim[0][0] = amin, lim[0][1] = amax, lim[0][2] = da 
                           // lim[1][0] = lnkmin, lim[1][1] = lnkmax, lim[1][2] = dlnk

  if (table == NULL || fdiff(cache[1], Ntable.random))
  {
    if (table != NULL) free(table);
    table = (double**) malloc2d(Ntable.N_a, Ntable.N_k_nlin);
    
    lim[0][0] = limits.a_min;
    lim[0][1] = 0.9999999;
    lim[0][2] = (lim[0][1] - lim[0][0]) / ((double) Ntable.N_a - 1.0);
    lim[1][0] = log(limits.k_min_cH0);
    lim[1][1] = log(limits.k_max_cH0);
    lim[1][2] = (lim[1][1] - lim[1][0]) / ((double) Ntable.N_k_nlin - 1.0);
  }

  if (fdiff(cache[0], cosmology.random) || 
      fdiff(cache[1], Ntable.random) ||
      fdiff(cache[2], nuisance.random_gas))
  { 
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init = p_xy_nointerp(exp(lim[1][0]), lim[0][0], 2, 1);
    }
    #pragma GCC diagnostic pop   

    #pragma omp parallel for collapse(2)
    for (int i=1; i<Ntable.N_a; i++) 
      for (int j=0; j<Ntable.N_k_nlin; j++) 
        table[i][j] = log(p_xy_nointerp(exp(lim[1][0] + j*lim[1][2]), 
                                            lim[0][0] + i*lim[0][2], 2, 0));

    cache[0] = cosmology.random;
    cache[1] = Ntable.random;
    cache[2] = nuisance.random_gas;
  }
  return exp(interpol2d(table, 
    Ntable.N_a, lim[0][0], lim[0][1], lim[0][2], a, 
    Ntable.N_k_nlin, lim[1][0], lim[1][1], lim[1][2], log(k)));
}

double p_gm_nointerp(
    const double k, 
    const double a, 
    const int ni,
    const int init
  )
{
  return Pdelta(k, a)*bgal(ni, a) + 
            GM02_nointerp(k, a, ni, init)/ngal(ni, a);
}

double p_gm(
    const double k, 
    const double a, 
    const int ni
  )
{
  static double cache[MAX_SIZE_ARRAYS];
  static double*** table = NULL;
  static double** lim = NULL; //lim[:,0] = amin; lim[:,1] = amax; lim[:,2] = da; 
                              //lim[redshift.clustering_nbin][0] = lnkmin; 
                              //lim[redshift.clustering_nbin][1] = lnkmax; 
                              //lim[redshift.clustering_nbin][2] = dlnk; 

  const int nbin = redshift.clustering_nbin;
  const int na = (int) Ntable.N_a/5.0; // range is the (\delta a) of a single bin
  
  if (table == NULL || fdiff(cache[1], Ntable.random))
  {
    if (table != NULL) free(table);
    table = (double***) malloc3d(nbin, na, Ntable.N_k_nlin);
    if (lim != NULL) free(lim);
    lim = (double**) malloc2d(nbin+1, 3);

    for (int l=0; l<redshift.clustering_nbin; l++) 
    {
      lim[l][0] = amin_lens(l);
      lim[l][1] = amax_lens(l);
      lim[l][2] = (lim[l][1] - lim[l][0])/((double) na - 1.0);
    }
    lim[nbin][0] = log(limits.k_min_cH0);
    lim[nbin][1] = log(limits.k_max_cH0);
    lim[nbin][2] = (lim[nbin][1]-lim[nbin][0])/((double) Ntable.N_k_nlin - 1.0);
  }

  if (fdiff(cache[0], cosmology.random) || 
      fdiff(cache[1], Ntable.random)    ||
      fdiff(cache[2], nuisance.random_galaxy_bias) ||
      fdiff(cache[3], redshift.random_clustering))
  { 
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init = p_gm_nointerp(exp(lim[nbin][0]), lim[0][0], 0, 1);
    }
    #pragma GCC diagnostic pop
    
    #pragma omp parallel for collapse(3)
    for (int l=0; l<redshift.clustering_nbin; l++) 
      for (int i=0; i<na; i++) 
        for (int j=0; j<Ntable.N_k_nlin; j++) 
          table[l][i][j] = log(p_gm_nointerp(exp(lim[nbin][0] + j*lim[nbin][2]), 
            lim[l][0] + i*lim[l][2], l, 0));

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
  return (a < lim[ni][0] || a > lim[ni][1]) ? 0.0 : exp(interpol2d(table[ni], 
    na, lim[ni][0], lim[ni][1], lim[ni][2], a, 
    Ntable.N_k_nlin, lim[nbin][0], lim[nbin][1], lim[nbin][2], log(k)));
}

double p_gg_nointerp(
    const double k, 
    const double a, 
    const int ni, 
    const int nj,
    const int init
  )
{
  if (ni != nj)
  {
    log_fatal("cross-tomography (ni,nj) = (%d,%d) bins not supported", ni, nj);
    exit(1);
  }
  const double bg = bgal(ni, a);
  const double ng = ngal(ni, a);
  return Pdelta(k, a)*bg*bg + G02_nointerp(k, a, ni, init)/(ng*ng);
}

double p_gg(
    const double k, 
    const double a, 
    const int ni, 
    const int nj
  )
{
  static double cache[MAX_SIZE_ARRAYS];
  static double*** table = NULL;
  static double** lim = NULL; //lim[0,:] = amin; lim[:,1] = amax; lim[:,2] = da; 
                              //lim[redshift.clustering_nbin] = lnkmin; 
                              //lim[redshift.clustering_nbin] = lnkmax; 
                              //lim[redshift.clustering_nbin] = dlnk; 
  
  const int nbin = redshift.clustering_nbin;
  const int na = (int) Ntable.N_a/5.0;

  if (table == NULL || fdiff(cache[1], Ntable.random))
  {
    if (table != NULL) free(table);
    table = (double***) malloc3d(nbin, na, Ntable.N_k_nlin);
    if (lim != NULL) free(lim);
    lim = (double**) malloc2d(nbin+1, 3);

    for (int l=0; l<redshift.clustering_nbin; l++) 
    {
      lim[l][0] = amin_lens(l);
      lim[l][1] = amax_lens(l);
      lim[l][2] = (lim[l][1] - lim[l][0])/((double) na - 1.0);
    }
    lim[nbin][0] = log(limits.k_min_cH0);
    lim[nbin][1] = log(limits.k_max_cH0);
    lim[nbin][2] = (lim[nbin][1]-lim[nbin][0])/((double) Ntable.N_k_nlin - 1.0);
  }

  if (fdiff(cache[0], cosmology.random) || 
      fdiff(cache[1], Ntable.random)    ||
      fdiff(cache[2], nuisance.random_galaxy_bias) ||
      fdiff(cache[3], redshift.random_clustering))
  { 
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    {
      double init = p_gg_nointerp(exp(lim[nbin][0]), lim[0][0], 0, 0, 1);
    }
    #pragma GCC diagnostic pop
    
    #pragma omp parallel for collapse(3)
    for (int l=0; l<nbin; l++) 
      for (int i=0; i<na; i++) 
        for (int j=0; j<Ntable.N_k_nlin; j++) 
          table[l][i][j] = log(p_gg_nointerp(exp(lim[nbin][0]+j*lim[nbin][2]),
            lim[l][0]+i*lim[l][2], l, l, 0));

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
  return (a < lim[ni][0] || a > lim[ni][1]) ? 0.0 : exp(interpol2d(table[ni], 
      na, lim[ni][0], lim[ni][1], lim[ni][2], a, 
      Ntable.N_k_nlin, lim[nbin][0], lim[nbin][1], lim[nbin][2], log(k)));
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// MISCELLANEOUS
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void set_HOD(const int ni)
{ 
  const double z = zmean(ni);
  const double a = 1.0/(z + 1.0);
  
  // Parameterization of Zehavi et al. 
  // hod[zi][] = {lg(M_min), sigma_{lg M}, lg M_1, lg M_0, alpha, f_c}
  // gbias.gc[] = {f_g} (shift of concentration parameter: c_g(M) = f_g c(M))
  
  // Values from Coupon etal. (2012) for red gals with M_r < -21.8 (Table B.2)
  switch (ni)
  {
    case 0:
    {
      nuisance.hod[0][0] = 13.17;
      nuisance.hod[0][1] = 0.39;
      nuisance.hod[0][2] = 14.53;
      nuisance.hod[0][3] = 11.09;
      nuisance.hod[0][4] = 1.27;
      nuisance.hod[0][5] = 1.00;
      nuisance.gb[0][ni] = bgal(ni, a);      
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
      nuisance.gb[0][ni] = hm_funcs_nointerp(ni, a, 3, 0);
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
      nuisance.gb[0][ni] = hm_funcs_nointerp(ni, a, 3, 0);
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
      nuisance.gb[0][ni] = hm_funcs_nointerp(ni, a, 3, 0);
      break;
    }
    case 4:
    { // no information for higher redshift populations - copy 1<z<1.2 values
      nuisance.hod[4][0] = 12.80;
      nuisance.hod[4][1] = 0.35;
      nuisance.hod[4][2] = 13.94;
      nuisance.hod[4][3] = 12.15;
      nuisance.hod[4][4] = 1.52;
      nuisance.hod[4][5] = 1.00;
      nuisance.gb[0][ni] = hm_funcs_nointerp(ni, a, 3, 0);
      break;
    }
    default:
    {
      log_fatal("no HOD parameters specified to initialize bin %d\n", ni);
      exit(1);
    }
  }

  log_debug("HOD: bin %d; <z> %.2f; <n_g> %e(h/Mpc)^3", ni, z, 
    ngal_nointerp(ni, a, 0)*pow(cosmology.coverH0, -3.0));
  
  log_debug("HOD: bin %d; <z> %.2f; <M> h/Msun %.4e", ni, z, mmean_nointerp(ni,a,0));
  
  log_debug("HOD: bin %d; <z> %.2f; f_sat %.3f", ni, z, fsat_nointerp(ni,a,0));
  
  log_debug("HOD: bin %d; <z> %.2f; <b_g> %.2f", ni, z, nuisance.gb[0][ni]);
}



// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// External Power Spectrum
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double external_p_gm(const double k, const double a)
{
  // convert from (x/Mpc/h - dimensioneless) to h/Mpc with x = c/H0 (Mpc)
  const double log10k = log10(k/cosmology.coverH0);
  const double z = 1.0/a-1.0;
  
  int i = 0;
  {
    size_t ilo = 0;
    size_t ihi = hmemu.lnPGM_nk-1;
    while (ihi>ilo+1) 
    {
      size_t ll = (ihi+ilo)/2;
      if(hmemu.lnPGM[ll][hmemu.lnPGM_nz] > log10k)
        ihi = ll;
      else
        ilo = ll;
    }
    i = ilo;
  }

  int j = 0;
  {
    size_t ilo = 0;
    size_t ihi = hmemu.lnPGM_nz-1;
    while (ihi>ilo+1) 
    {
      size_t ll = (ihi+ilo)/2;
      if(hmemu.lnPGM[hmemu.lnPGM_nk][ll] > z)
        ihi = ll;
      else
        ilo = ll;
    }
    j = ilo;
  }

  double dx = (log10k                           - hmemu.lnPGM[i][hmemu.lnPGM_nz])/
              (hmemu.lnPGM[i+1][hmemu.lnPGM_nz] - hmemu.lnPGM[i][hmemu.lnPGM_nz]);

  double dy = (z                                - hmemu.lnPGM[hmemu.lnPGM_nk][j])/
              (hmemu.lnPGM[hmemu.lnPGM_nk][j+1] - hmemu.lnPGM[hmemu.lnPGM_nk][j]);

  const double out_lnP =    (1-dx)*(1-dy)*hmemu.lnPGM[i][j]
                          + (1-dx)*dy*hmemu.lnPGM[i][j+1]
                          + dx*(1-dy)*hmemu.lnPGM[i+1][j]
                          + dx*dy*hmemu.lnPGM[i+1][j+1];

  // convert from (Mpc/h)^3 to (Mpc/h)^3/(c/H0=100)^3 (dimensioneless)
  return exp(out_lnP)/(cosmology.coverH0*cosmology.coverH0*cosmology.coverH0);
}

double external_p_gg(const double k, const double a)
{
  // convert from (x/Mpc/h - dimensioneless) to h/Mpc with x = c/H0 (Mpc)
  const double log10k = log10(k/cosmology.coverH0);
  const double z = 1.0/a-1.0;
  
  int i = 0;
  {
    size_t ilo = 0;
    size_t ihi = hmemu.lnPGG_nk-1;
    while (ihi>ilo+1) 
    {
      size_t ll = (ihi+ilo)/2;
      if(hmemu.lnPGG[ll][hmemu.lnPGG_nz] > log10k)
        ihi = ll;
      else
        ilo = ll;
    }
    i = ilo;
  }

  int j = 0;
  {
    size_t ilo = 0;
    size_t ihi = hmemu.lnPGG_nz-1;
    while (ihi>ilo+1) 
    {
      size_t ll = (ihi+ilo)/2;
      if(hmemu.lnPGG[hmemu.lnPGG_nk][ll] > z)
        ihi = ll;
      else
        ilo = ll;
    }
    j = ilo;
  }

  double dx = (log10k                           - hmemu.lnPGG[i][hmemu.lnPGG_nz])/
              (hmemu.lnPGG[i+1][hmemu.lnPGG_nz] - hmemu.lnPGG[i][hmemu.lnPGG_nz]);

  double dy = (z                                - hmemu.lnPGG[hmemu.lnPGG_nk][j])/
              (hmemu.lnPGG[hmemu.lnPGG_nk][j+1] - hmemu.lnPGG[hmemu.lnPGG_nk][j]);

  const double out_lnP =    (1-dx)*(1-dy)*hmemu.lnPGG[i][j]
                          + (1-dx)*dy*hmemu.lnPGG[i][j+1]
                          + dx*(1-dy)*hmemu.lnPGG[i+1][j]
                          + dx*dy*hmemu.lnPGG[i+1][j+1];

  // convert from (Mpc/h)^3 to (Mpc/h)^3/(c/H0=100)^3 (dimensioneless)
  return exp(out_lnP)/(cosmology.coverH0*cosmology.coverH0*cosmology.coverH0);
}