#include <assert.h>
#include <gsl/gsl_sf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "basics.h"
#include "halo.h"
#include "math.h"
#include "recompute.h"
#include "redmagic.h"
#include "structs.h"

double n_c_tab(double mh) {
  static double **tab = 0;
  static int Nm = 50;
  static double Mmin = 10., Mmax = 15., dm = 0.1;
  if (tab == 0) {
    FILE *f;
    int i;
    float f1, f2, f3;
    f = fopen(redm.HOD_FILE, "r");
    tab = create_double_matrix(0, 1, 0, Nm - 1);
    for (i = 0; i < Nm; i++) {
      fscanf(f, "%e %e %e\n", &f1, &f2, &f3);
      tab[0][i] = f2;
    }
    fclose(f);
    printf("n_c_tab\n");
  }
  if (log10(mh) <= Mmin || log10(mh) > Mmax) {
    return 0.;
  }
  return tab[0][(int)floor((log10(mh) - Mmin) / dm)];
}

double n_s_tab(double mh) {
  static double **tab = 0;
  static int Nm = 50;
  static double Mmin = 10., Mmax = 15., dm = 0.1;
  if (tab == 0) {
    FILE *f;
    int i;
    float f1, f2, f3;
    f = fopen(redm.HOD_FILE, "r");
    tab = create_double_matrix(0, 1, 0, Nm - 1);
    for (i = 0; i < Nm; i++) {
      fscanf(f, "%e %e %e\n", &f1, &f2, &f3);
      tab[0][i] = f3;
    }
    fclose(f);
    printf("n_s_tab\n");
  }
  if (log10(mh) <= Mmin || log10(mh) > Mmax) {
    return 0.;
  }
  return tab[0][(int)floor((log10(mh) - Mmin) / dm)];
}

double n_c_Reddick(double mh) {
  return redm.fc * 0.5 *
         (1.0 + gsl_sf_erf((log10(mh) - redm.hod[0]) / redm.hod[1]));
}

double n_s_Reddick(double mh) {
  return pow(mh / pow(10.0, redm.hod[2]), redm.hod[4]) *
         exp(-pow(10.0, redm.hod[3]) / mh);
}

double n_c_rm(double mh) {
  if (redm.parameterization == 3) {
    return (mh > pow(10.0, redm.hod[0])) ? 1. : 0.;
  }
  if (redm.parameterization == 2) {
    return n_c_tab(mh);
  }
  if (redm.parameterization == 1) {
    return n_c_Reddick(mh);
  }
  if (redm.hod[0] < 10 || redm.hod[0] > 16) {
    printf("n_c_rm: HOD parameters not set\n");
    exit(EXIT_FAILURE);
  }
  // default: Zehavi et al. + central fraction
  return redm.fc * 0.5 *
         (1.0 + gsl_sf_erf((log10(mh) - redm.hod[0]) / redm.hod[1]));
}

double n_s_rm(double mh) {
  if (redm.parameterization == 3) {
    return 0.0;
  }
  if (redm.parameterization == 2) {
    return n_s_tab(mh);
  }
  if (redm.parameterization == 1) {
    return n_s_Reddick(mh);
  }
  double ns =
      pow(fmax(mh - pow(10.0, redm.hod[3]), 0.) / pow(10.0, redm.hod[2]),
          redm.hod[4]) *
      n_c_rm(mh) / redm.fc;
  if (ns > 0) {
    return ns;
  } else {
    return 1.e-15;
  }
}

double int_bgal_rm(double m, void *params) {
  double *array = (double *)params;
  double a = array[0];
  return massfunc(exp(m), a) * B1(exp(m), a) *
         (n_c_rm(exp(m)) + n_s_rm(exp(m))) * exp(m);
}

double int_fsat_rm(double m, void *params) {
  double *array = (double *)params;
  double a = array[0];
  return massfunc(exp(m), a) * (n_s_rm(exp(m))) * exp(m);
}

double int_nc_rm(double m, void *params) {
  double *array = (double *)params;
  double a = array[0];
  return massfunc(exp(m), a) * n_c_rm(exp(m)) * exp(m);
}

double int_ngal_rm(double m, void *params) {
  double *array = (double *)params;
  double a = array[0];
  return massfunc(exp(m), a) * (n_c_rm(exp(m)) + n_s_rm(exp(m))) * exp(m);
}

double int_mc_rm(double m, void *params) {
  double *array = (double *)params;
  double a = array[0];
  return massfunc(exp(m), a) * n_c_rm(exp(m)) * exp(m) * exp(m);
}

double int_mmean_rm(double m, void *params) {
  double *array = (double *)params;
  double a = array[0];
  return massfunc(exp(m), a) * (n_c_rm(exp(m)) + n_s_rm(exp(m))) * exp(m) *
         exp(m);
}

double ngal_rm(double a) {

  static double **table;
  static cosmopara C;

  static double da = 0.0;
  static double M0 = 0.0;
  double aa;
  int i;
  double array[1];

  if (M0 != redm.hod[0] || recompute_cosmo3D(C)) {
    if (table == 0) {
      table = create_double_matrix(0, 1, 0, Ntable.N_a - 1);
      da = (1. / (tomo.clustering_zmin[0] + 1.) -
            1. / (tomo.clustering_zmax[tomo.clustering_Nbin - 1] + 1.)) /
           (Ntable.N_a - 1);
    }
    aa = 1. / (tomo.clustering_zmax[tomo.clustering_Nbin - 1] + 1.);
    for (i = 0; i < Ntable.N_a; i++, aa += da) {
      array[0] = aa;
      table[0][i] = int_gsl_integrate_medium_precision(
          int_ngal_rm, (void *)array, log(10.) * (redm.hod[0] - 2.),
          log(limits.M_max), NULL, 5000);
    }
    M0 = redm.hod[0];
    update_cosmopara(&C);
  }
  return interpol(table[0], Ntable.N_a,
                  1. / (tomo.clustering_zmax[tomo.clustering_Nbin - 1] + 1.),
                  1. / (tomo.clustering_zmin[0] + 1.), da, a, 1.0, 1.0);
}

double bgal_rm(double a) {
  static double **table;
  static cosmopara C;
  static double M0 = 0.0;
  static double da = 0.0;
  double aa;
  int i;
  double array[1];
  if (M0 != redm.hod[0] || recompute_cosmo3D(C)) {
    if (table == 0) {
      table = create_double_matrix(0, 1, 0, Ntable.N_a - 1);
      da = (1. / (tomo.clustering_zmin[0] + 1.) -
            1. / (tomo.clustering_zmax[tomo.clustering_Nbin - 1] + 1.)) /
           (Ntable.N_a - 1);
    }
    //j = 0;
    aa = 1. / (tomo.clustering_zmax[tomo.clustering_Nbin - 1] + 1.);
    for (i = 0; i < Ntable.N_a; i++, aa += da) {
      array[0] = aa;
      table[0][i] =
          int_gsl_integrate_medium_precision(int_bgal_rm, (void *)array,
                                             log(10.) * (redm.hod[0] - 2.),
                                             log(limits.M_max), NULL, 5000) /
          int_gsl_integrate_medium_precision(int_ngal_rm, (void *)array,
                                             log(10.) * (redm.hod[0] - 2.),
                                             log(limits.M_max), NULL, 5000);
    }
    M0 = redm.hod[0];
    update_cosmopara(&C);
  }
  return interpol(table[0], Ntable.N_a,
                  1. / (tomo.clustering_zmax[tomo.clustering_Nbin - 1] + 1.),
                  1. / (tomo.clustering_zmin[0] + 1.), da, a, 1.0, 1.0);
}

double mmean_rm(double a) {
  double array[2] = {a, 0.0};
  return int_gsl_integrate_medium_precision(int_mmean_rm, (void *)array,
                                            log(10.) * (redm.hod[0] - 2.),
                                            log(limits.M_max), NULL, 5000) /
         int_gsl_integrate_medium_precision(int_ngal_rm, (void *)array,
                                            log(10.) * (redm.hod[0] - 2.),
                                            log(limits.M_max), NULL, 5000);
}

double mc_rm(double a) {
  double array[2] = {a, 0.0};
  return int_gsl_integrate_medium_precision(int_mc_rm, (void *)array,
                                            log(10.) * (redm.hod[0] - 2.),
                                            log(limits.M_max), NULL, 5000) /
         int_gsl_integrate_medium_precision(int_nc_rm, (void *)array,
                                            log(10.) * (redm.hod[0] - 2.),
                                            log(limits.M_max), NULL, 5000);
}

double fsat_rm(double a) {
  double array[2] = {a, 0.0};
  return int_gsl_integrate_medium_precision(int_fsat_rm, (void *)array,
                                            log(10.) * (redm.hod[0] - 2.),
                                            log(limits.M_max), NULL, 5000) /
         int_gsl_integrate_medium_precision(int_ngal_rm, (void *)array,
                                            log(10.) * (redm.hod[0] - 2.),
                                            log(limits.M_max), NULL, 5000);
}

/********************* galaxy power spectrum from halo model + HOD
 * ***************************/
/* ******************** no halo exclusion or radial biasing so far
 * ***************************/

double u_g_rm(double k, double m,
              double a) // Fourier transformed normalized galaxy density
                        // profile, NFW with rescaled concentraion so far
{
  // analytic FT of NFW profile, from Cooray & Sheth 01
  if (redm.cg <= 0) {
    printf("galaxy concentration parameter ill-defined, exiting.\n");
    exit(EXIT_FAILURE);
  }
  return u_nfw_c(conc(m, a) * redm.cg, k, m, a);
}

/**** 1-halo galaxy-galaxy spectrum ******/
double int_for_G02_rm(double logm, void *para) {
  double *array = (double *)para;
  double m = exp(logm);
  double u, ns, k = array[0], a = array[1];
  u = u_g_rm(k, m, a);
  ns = n_s_rm(m);
  // printf("%e %e\n",k,m);
  return massfunc(m, a) * m *
         (u * u * ns * ns * 1.0 + 2.0 * u * ns * n_c_rm(m));
}

double G02_rm(double k, double a) { // needs to be devided by ngal(nz, a)^2
  double array[2] = {k, a};
  return int_gsl_integrate_medium_precision(int_for_G02_rm, (void *)array,
                                            log(10.) * (redm.hod[0] - 2.),
                                            log(limits.M_max), NULL, 5000);
}

double P_gg_rm(double k, double a) {
  static double logkmin = 0., logkmax = 0., dk = 0., da = 0., amin = 0.,
                amax = 0.;
  static int N_a = 10, N_k_nlin = 50;
  static double M0 = 0.0;
  static cosmopara C;

  static double **table_P_gg = 0;

  double klog, val, aa, kk;
  int i, j;

  if (M0 != redm.hod[0] ||
      recompute_cosmo3D(
          C)) { // extend this by halo model parameters if these are sampled
                // independently of cosmology parameters
    M0 = redm.hod[0];
    update_cosmopara(&C);
    if (table_P_gg == 0) {
      table_P_gg = create_double_matrix(0, N_a - 1, 0, N_k_nlin - 1);
    }
    amax = fmin(1. / (1 + redshift.clustering_zdistrpar_zmin - 0.01), 0.999);
    amin = 1. / (1 + redshift.clustering_zdistrpar_zmax + 0.01);
    da = (amax - amin) / (N_a - 1.);
    logkmin = log(limits.k_min_cH0);
    logkmax = log(limits.k_max_cH0);
    // printf("tabulating P_gg %e %e, %e %e\n", exp(logkmin),exp(logkmax),
    // amin,amax);
    dk = (logkmax - logkmin) / (N_k_nlin - 1.);
    aa = amin;
    for (i = 0; i < N_a; i++, aa += da) {
      klog = logkmin;
      for (j = 0; j < N_k_nlin; j++, klog += dk) {
        kk = exp(klog);
        //        table_P_gg[i][j] =
        //        log(Pdelta(kk,aa)*pow(bgal_rm(aa),2.0)+G02_rm(kk,aa)/pow(ngal_rm(aa),2.0));
        table_P_gg[i][j] = log(G02_rm(kk, aa) / pow(ngal_rm(aa), 2.0));
      }
    }
  }
  klog = log(k);
  if (klog <= logkmin || klog > logkmax - dk) {
    return 0.;
  }
  val = interpol2d(table_P_gg, N_a, amin, amax, da, a, N_k_nlin, logkmin,
                   logkmax, dk, klog, 0.0, 0.0);
  return exp(val);
}

/*
double int_for_C_cl_HOD_rm(double a, void *params)
{
  double res,ell, fK, k;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_p_2");

  ell       = ar[0]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;

  res=
W_HOD(a,ar[1])*W_HOD(a,ar[1])*dchi_da(a)/fK/fK;//pf_photoz(1./a-1.,(int)ar[1])*pf_photoz(1./a-1.,(int)ar[1])/(a*a)*hoverh0(a)/fK/fK;
  if (res !=0){res= res*P_gg_rm(k,a);}
  return res;
}

double C_cl_HOD_rm_tomo (double l, int nz)
{
  double da,ap,res, array[2];
  array[0] = l; array[1] = (double)nz;
  return
int_gsl_integrate_low_precision(int_for_C_cl_HOD_rm,(void*)array,1./(1+tomo.clustering_zmax[nz]),1./(1+tomo.clustering_zmin[nz]),NULL,1000);
}
*/

/**** 1-halo galaxy-matter spectrum ******/
double int_GM02_rm(double logm, void *para) {
  double *array = (double *)para;
  double m = exp(logm);
  double k = array[0], a = array[1];
  return massfunc(m, a) * m * m / (cosmology.rho_crit * cosmology.Omega_m) *
         u_nfw_c(conc(m, a), k, m, a) *
         (u_g_rm(k, m, a) * n_s_rm(m) + n_c_rm(m)); // u_nfw_c(conc(m,a),k,m,a)
}

double GM02_rm(double k, double a) { // needs to be devided by ngal_rm(nz, a)
  double array[2] = {k, a};
  return int_gsl_integrate_medium_precision(int_GM02_rm, (void *)array,
                                            log(10.) * (redm.hod[0] - 2.),
                                            log(limits.M_max), NULL, 5000);
}

double P_gm_rm(double k, double a) {
  static double logkmin = 0., logkmax = 0., dk = 0., da = 0., amin = 0.,
                amax = 0.;
  static int N_a = 10, N_k_nlin = 50;
  static double M0 = 0.0;
  static cosmopara C;

  static double **table_P_gm = 0;

  double klog, val, aa, kk;
  int i, j;

  if (M0 != redm.hod[0] ||
      recompute_cosmo3D(
          C)) { // extend this by halo model parameters if these are sampled
                // independently of cosmology parameters
    M0 = redm.hod[0];
    update_cosmopara(&C);
    if (table_P_gm == 0) {
      table_P_gm = create_double_matrix(0, N_a - 1, 0, N_k_nlin - 1);
    }
    amax = fmin(1. / (1 + tomo.clustering_zmin[0] - 0.01), 0.999);
    amin = 1. / (1 + tomo.clustering_zmax[tomo.clustering_Nbin - 1] + 0.01);
    da = (amax - amin) / (N_a - 1.);
    logkmin = log(limits.k_min_cH0);
    logkmax = log(limits.k_max_cH0);
    dk = (logkmax - logkmin) / (N_k_nlin - 1.);
    aa = amin;
    for (i = 0; i < N_a; i++, aa += da) {
      klog = logkmin;
      for (j = 0; j < N_k_nlin; j++, klog += dk) {
        kk = exp(klog);
        //        table_P_gm[i][j] =
        //        log(Pdelta(kk,aa)*bgal_rm(aa)+GM02_rm(kk,aa)/ngal_rm(aa));
        table_P_gm[i][j] = log(GM02_rm(kk, aa) / ngal_rm(aa));
      }
    }
  }
  // printf("P_gm\n");
  klog = log(k);
  val = interpol2d(table_P_gm, N_a, amin, amax, da, a, N_k_nlin, logkmin,
                   logkmax, dk, klog, 0.0, 0.0);
  return exp(val);
}

/*
double int_for_C_gl_HOD_rm(double a, void *params)
{
  double res,ell, fK, k;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_p_2");

  ell       = ar[0]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;


  res=
W_HOD(a,ar[1])*W_kappa(a,fK,ar[2])*dchi_da(a)/fK/fK;//pf_photoz(1./a-1.,(int)
(ar[1]))*g_tomo(a,(int) ar[2])/(a*a*a)/fK; if (res !=0){res= res*P_gm_rm(k,a);}
  return res;
}

double C_gl_HOD_rm_tomo (double l, int nzl, int nzs){
  double da,ap,res, array[3] = {l, (double)nzl,(double) nzs};
  return
int_gsl_integrate_low_precision(int_for_C_gl_HOD_rm,(void*)array,1./(1+tomo.clustering_zmax[nzl]),1./(1+tomo.clustering_zmin[nzl]),NULL,1000);
}
*/
