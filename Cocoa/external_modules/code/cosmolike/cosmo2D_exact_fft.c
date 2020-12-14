#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../cfftlog/cfftlog.h"

#include "basics.h"
#include "cosmo2D_exact_fft.h"
#include "cosmo2D_fourier.h"
#include "cosmo3D.h"
#include "radial_weights.h"
#include "recompute.h"
#include "redshift_spline.h"
#include "structs.h"

#include "../log.c/src/log.h"

double G_taper(double k) {
  double s_bao = 5.5 / cosmology.coverH0;
  return exp(-k * k * s_bao * s_bao);
}

double int_for_C_cl_lin(double a, void *params) {
  double res, ell, fK, k;
  double *ar = (double *)params;

  struct chis chidchi = chi_all(a);
  double hoverh0 = hoverh0v2(a, chidchi.dchida);

  ell = ar[2] + 0.5;
  fK = f_K(chidchi.chi);
  k = ell / fK;

  res = W_gal(a, ar[0], chidchi.chi, hoverh0) *
    W_gal(a, ar[1], chidchi.chi, hoverh0) * chidchi.dchida / (fK * fK);

  return res * p_lin(k, a);
}

// galaxy clustering power spectrum of galaxy bins ni, nj
double C_cl_lin_nointerp(double l, int ni, int nj) {
  double array[3] = {1.0 * ni, 1.0 * nj, l};
  return int_gsl_integrate_low_precision(int_for_C_cl_lin, (void *)array,
                fmax(amin_lens(ni), amin_lens(nj)), 0.99999, NULL, 1000);
}

// test for replacing Plin with Pdelta and rescaled
double int_for_C_cl_nl_rescale(double a, void *params) {
  double res, ell, fK, k;
  double *ar = (double *)params;

  double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  double hoverh0 = hoverh0v2(a, chidchi.dchida);

  ell = ar[2] + 0.5;
  fK = f_K(chidchi.chi);
  k = ell / fK;

  res = W_gal(a, ar[0], chidchi.chi, hoverh0) *
    W_gal(a, ar[1], chidchi.chi, hoverh0) * chidchi.dchida / (fK * fK);
  return res * Pdelta(k, 0.9999) * growfac_a * growfac_a;
}

// galaxy clustering power spectrum of galaxy bins ni, nj
double C_cl_nl_rescaled_nointerp(double l, int ni, int nj) {
  double array[3] = {1.0 * ni, 1.0 * nj, l};
  return int_gsl_integrate_low_precision(int_for_C_cl_nl_rescale,
    (void *)array, fmax(amin_lens(ni), amin_lens(nj)), 0.99999, NULL, 1000);
}

// Integrand for galaxy density
void f_chi_for_Psi_cl(double *chi_ar, int Nchi, double *f_chi_ar, int ni) {
  const double real_coverH0 = cosmology.coverH0/cosmology.h0; // unit Mpc
  for (int i = 0; i < Nchi; i++) {
    // first convert unit of chi from Mpc to c/H0
    const double a = a_chi(chi_ar[i]/real_coverH0);
    const double z = 1. / a - 1.;
    if ((z < tomo.clustering_zmin[ni]) || (z > tomo.clustering_zmax[ni])) {
      f_chi_ar[i] = 0.;
    } else {
      // get rid of unphysical negatives
      const double tmp1 = pf_photoz(z, ni);
      const double pf = (tmp1 < 0.) ? 0 : tmp1;
      f_chi_ar[i] = chi_ar[i] * pf * growfac(a) *
                    gbias.b1_function(z, ni) * hoverh0(a)/real_coverH0;
    }
  }
}

// Integrand for galaxy density RSD
void f_chi_for_Psi_cl_RSD(double *chi_ar, int Nchi, double *f_chi_RSD_ar,
int ni) {
  const double real_coverH0 = cosmology.coverH0 / cosmology.h0;
  for (int i = 0; i < Nchi; i++) {
    // first convert unit of chi from Mpc to c/H0
    const double a = a_chi(chi_ar[i] / real_coverH0);
    const double z = 1. / a - 1.;
    if ((z < tomo.clustering_zmin[ni]) || (z > tomo.clustering_zmax[ni])) {
      f_chi_RSD_ar[i] = 0.;
    } else {
      // get rid of unphysical negatives
      const double tmp1 = pf_photoz(z, ni);
      const double pf = (tmp1 < 0.) ? 0 : tmp1;
      struct growths tmp2 = growfac_all(a);
      f_chi_RSD_ar[i] = -chi_ar[i] * pf * tmp2.D * tmp2.f *
                        hoverh0(a) / real_coverH0;
    }
  }
}

// Integrand for lensing magnification of galaxy density
void f_chi_for_Psi_cl_Mag(double *chi_ar, int Nchi, double *f_chi_Mag_ar,
int ni) {
  const double real_coverH0 = cosmology.coverH0 / cosmology.h0;
  for (int i = 0; i < Nchi; i++) {
    // first convert unit of chi from Mpc to c/H0
    const double a = a_chi(chi_ar[i] / real_coverH0);
    const double z = 1. / a - 1.;
    const double fK = f_K(chi_ar[i] / real_coverH0);
    if ((z > tomo.clustering_zmax[ni])) {
      f_chi_Mag_ar[i] = 0.;
    } else {
      const double wmag = W_mag(a, fK, (double)ni);
      const double window_M = wmag / fK / (real_coverH0 * real_coverH0);
      f_chi_Mag_ar[i] = window_M * growfac(a); // unit [Mpc^-2]
    }
  }
}

// Mixture of non-Limber and Limber of C_cl (galaxy clustering)
void C_cl_mixed(int L, int LMAX, int ni, int nj, double *Cl, double dev,
double tolerance) {

  static double **k1_ar, **k2_ar, **Fk1_ar, **Fk2_ar;
  static double **Fk1_Mag_ar, **Fk2_Mag_ar;
  static double *chi_ar;

  const int Nell_block = 20; //COCOA: original value = 100
  const int Nchi = 500; //COCOA: original value = 1000
  const int LMAX_NOLIMBER = 400;
  int ell_ar[Nell_block];
  double f1_chi_ar[Nchi], f2_chi_ar[Nchi];
  double f1_chi_RSD_ar[Nchi], f2_chi_RSD_ar[Nchi];
  double f1_chi_Mag_ar[Nchi], f2_chi_Mag_ar[Nchi];

  const double real_coverH0 = cosmology.coverH0 / cosmology.h0;
  const double chi_min = chi(1./(1.+0.002))*real_coverH0;
  const double chi_max = chi(1./(1.+4.))*real_coverH0;
  const double dlnchi = log(chi_max / chi_min) / (Nchi - 1.);
  const double dlnk = dlnchi;

  // COCOA: no need to create/destroy arrays with same size at every call
  if(k1_ar == 0) {
    k1_ar = malloc(Nell_block * sizeof(double*));
    k2_ar = malloc(Nell_block * sizeof(double*));
    Fk1_ar = malloc(Nell_block * sizeof(double*));
    Fk2_ar = malloc(Nell_block * sizeof(double*));
    Fk1_Mag_ar = malloc(Nell_block * sizeof(double*));
    Fk2_Mag_ar = malloc(Nell_block * sizeof(double*));
    for (int i = 0; i < Nell_block; i++) {
      k1_ar[i] = malloc(Nchi * sizeof(double));
      k2_ar[i] = malloc(Nchi * sizeof(double));
      Fk1_ar[i] = malloc(Nchi * sizeof(double));
      Fk2_ar[i] = malloc(Nchi * sizeof(double));
      Fk1_Mag_ar[i] = malloc(Nchi * sizeof(double));
      Fk2_Mag_ar[i] = malloc(Nchi * sizeof(double));
    }
    chi_ar = malloc(Nchi * sizeof(double));
  }

  for (int i = 0; i < Nchi; i++) {
    // chi_min and chi_max are cosmology dependent
    chi_ar[i] = chi_min * exp(dlnchi * i);
  }

  #pragma omp parallel for
  for (int i = 0; i < Nell_block; i++) {
    for (int j = 0; j < Nchi; j++) {
      k1_ar[i][j] = 0.0;
      k2_ar[i][j] = 0.0;
      Fk1_ar[i][j] = 0.;
      Fk2_ar[i][j] = 0.;
      Fk1_Mag_ar[i][j] = 0.;
      Fk2_Mag_ar[i][j] = 0.;
    }
  }

  f_chi_for_Psi_cl(chi_ar, Nchi, f1_chi_ar, ni);
  if (ni != nj) {
    f_chi_for_Psi_cl(chi_ar, Nchi, f2_chi_ar, nj);
  }

  f_chi_for_Psi_cl_RSD(chi_ar, Nchi, f1_chi_RSD_ar, ni);
  if (ni != nj) {
    f_chi_for_Psi_cl_RSD(chi_ar, Nchi, f2_chi_RSD_ar, nj);
  }

  f_chi_for_Psi_cl_Mag(chi_ar, Nchi, f1_chi_Mag_ar, ni);
  if (ni != nj) {
    f_chi_for_Psi_cl_Mag(chi_ar, Nchi, f2_chi_Mag_ar, nj);
  }

  int i_block = 0;

  config my_config, my_config_RSD, my_config_Mag;
  my_config.nu = 1.;
  my_config.c_window_width = 0.25;
  my_config.derivative = 0;
  my_config.N_pad = 200;
  my_config.N_extrap_low = 0;
  my_config.N_extrap_high = 0;

  my_config_RSD.nu = 1.01;
  my_config_RSD.c_window_width = 0.25;
  my_config_RSD.derivative = 2;
  my_config_RSD.N_pad = 500;
  my_config_RSD.N_extrap_low = 0;
  my_config_RSD.N_extrap_high = 0;

  my_config_Mag.nu = 1.;
  my_config_Mag.c_window_width = 0.25;
  my_config_Mag.derivative = 0;
  my_config_Mag.N_pad = 500;
  my_config_Mag.N_extrap_low = 0;
  my_config_Mag.N_extrap_high = 0;

  // COCOA: LMAX_NOLIMBER might avoid infinite loop in weird models
  while ((fabs(dev) > tolerance) && (L < LMAX_NOLIMBER)) {
    for (int i = 0; i < Nell_block; i++) {
      ell_ar[i] = i + i_block * Nell_block;
    }

    cfftlog_ells(chi_ar, f1_chi_ar, Nchi, &my_config, ell_ar, Nell_block, k1_ar,
                 Fk1_ar);
    if (ni != nj) {
      cfftlog_ells(chi_ar, f2_chi_ar, Nchi, &my_config, ell_ar, Nell_block,
                   k2_ar, Fk2_ar);
    }

    cfftlog_ells_increment(chi_ar, f1_chi_RSD_ar, Nchi, &my_config_RSD, ell_ar,
                           Nell_block, k1_ar, Fk1_ar);
    if (ni != nj) {
      cfftlog_ells_increment(chi_ar, f2_chi_RSD_ar, Nchi, &my_config_RSD,
                             ell_ar, Nell_block, k2_ar, Fk2_ar);
    }

    // Add in lensing magnification contribution
    cfftlog_ells(chi_ar, f1_chi_Mag_ar, Nchi, &my_config_Mag, ell_ar,
                 Nell_block, k1_ar, Fk1_Mag_ar);
    if (ni != nj) {
      cfftlog_ells(chi_ar, f2_chi_Mag_ar, Nchi, &my_config_Mag, ell_ar,
                   Nell_block, k2_ar, Fk2_Mag_ar);
    }

    #pragma omp parallel for
    for (int i = 0; i < Nell_block; i++) {
      const double ell_prefactor = ell_ar[i] * (ell_ar[i] + 1.);
      for (int j = 0; j < Nchi; j++) {
        Fk1_ar[i][j] += (ell_prefactor / (k1_ar[i][j] * k1_ar[i][j]) *
                         (gbias.b_mag[ni]) * Fk1_Mag_ar[i][j]);
        if (ni != nj) {
          Fk2_ar[i][j] += (ell_prefactor / (k2_ar[i][j] * k2_ar[i][j]) *
                           (gbias.b_mag[nj]) * Fk2_Mag_ar[i][j]);
        }
      }
    }

    {
      const int i = 0;
      double tmp[Nchi];
      #pragma omp parallel for
      for (int j = 0; j < Nchi; j++) {
        const double k1_cH0 = k1_ar[i][j] * real_coverH0;
        if (ni == nj) {
          tmp[j] = (Fk1_ar[i][j]) * (Fk1_ar[i][j]) * k1_cH0 * k1_cH0 *
                     k1_cH0 * p_lin(k1_cH0, 1.0);
        } else {
          tmp[j] = (Fk1_ar[i][j]) * (Fk2_ar[i][j]) * k1_cH0 * k1_cH0 *
                     k1_cH0 * p_lin(k1_cH0, 1.0);
        }
      }
      double cl_temp = 0.;
      for (int j = 0; j < Nchi; j++) {
        cl_temp += tmp[j];
      }
      Cl[ell_ar[i]] = cl_temp * dlnk * 2. / M_PI +
                      C_cl_tomo_nointerp(1. * ell_ar[i], ni, nj) -
                      C_cl_lin_nointerp(1. * ell_ar[i], ni, nj);
    }
    #pragma omp parallel for
    for (int i = 1; i < Nell_block; i++) {
      double cl_temp = 0.;
      for (int j = 0; j < Nchi; j++) {
        const double k1_cH0 = k1_ar[i][j] * real_coverH0;
        if (ni == nj) {
          cl_temp += (Fk1_ar[i][j]) * (Fk1_ar[i][j]) * k1_cH0 * k1_cH0 *
                     k1_cH0 * p_lin(k1_cH0, 1.0);
        } else {
          cl_temp += (Fk1_ar[i][j]) * (Fk2_ar[i][j]) * k1_cH0 * k1_cH0 *
                     k1_cH0 * p_lin(k1_cH0, 1.0);
        }
      }
      Cl[ell_ar[i]] = cl_temp * dlnk * 2. / M_PI +
                      C_cl_tomo_nointerp(1. * ell_ar[i], ni, nj) -
                      C_cl_lin_nointerp(1. * ell_ar[i], ni, nj);
    }

    i_block++;
    L = i_block * Nell_block - 1;
    dev = Cl[L] / C_cl_tomo_nointerp((double)L, ni, nj) - 1.;
  }
  L++;

  Cl[L] = C_cl_tomo((double)L, ni, nj);
  #pragma omp parallel for
  for (int l = L+1; l < LMAX; l++) {
    Cl[l] = C_cl_tomo((double)l, ni, nj);
  }
}

double w_tomo_nonLimber(int nt, int ni, int nj) {
  static int LMAX = 100000;
  static int NTHETA = 0;
  static double **Pl = 0;
  static double *w_vec = 0;
  static cosmopara C;
  static nuisancepara N;
  static galpara G;

  if (like.Ntheta == 0) {
    log_fatal("like.Ntheta not initialized");
    exit(1);
  }
  if (ni != nj) {
    log_fatal("ni != nj tomography not supported");
    exit(1);
  }

  if (Pl == 0) {
    NTHETA = like.Ntheta;
    Pl = create_double_matrix(0, like.Ntheta - 1, 0, LMAX - 1);
    w_vec = create_double_vector(0, tomo.clustering_Nbin * like.Ntheta - 1);
    double xmin[like.Ntheta];
    double xmax[like.Ntheta];
    const double logdt = (log(like.vtmax) - log(like.vtmin)) / like.Ntheta;

    #pragma omp parallel for
    for (int i = 0; i < NTHETA; i++) {
      double *Pmin = create_double_vector(0, LMAX+1);
      double *Pmax = create_double_vector(0, LMAX+1);

      xmin[i] = cos(exp(log(like.vtmin) + (i + 0.0) * logdt));
      xmax[i] = cos(exp(log(like.vtmin) + (i + 1.0) * logdt));

      gsl_sf_legendre_Pl_array(LMAX, xmin[i], Pmin);
      gsl_sf_legendre_Pl_array(LMAX, xmax[i], Pmax);

      Pl[i][0] = 1.0;

      const double a = (1.0/(xmin[i] - xmax[i]))*(1. / (4. * M_PI));
      for (int l = 1; l < LMAX; l++) {
        Pl[i][l] = a*(Pmin[l + 1] - Pmax[l + 1] - Pmin[l - 1] + Pmax[l - 1]);
      }

      free_double_vector(Pmin,0,LMAX+1);
      free_double_vector(Pmax,0,LMAX+1);
    }
  }

  if (recompute_clustering(C, G, N, ni, nj)) {
    const int L = 0;
    // required fractional accuracy in C(l)
    const double tolerance = 0.01;
    // initialize to large value in order to start while loop
    const double dev = 10. * tolerance;
    // dev will be the actual difference between exact and Limber calcuation
    double** Cl = malloc(tomo.clustering_Nbin * sizeof(double*));
    for (int nz = 0; nz < tomo.clustering_Nbin; nz++) {
      Cl[nz] = calloc(LMAX, sizeof(double));
      C_cl_mixed(L, LMAX, nz, nz, Cl[nz], dev, tolerance); // no threading
    }
    #pragma omp parallel for
    for (int nz = 0; nz < tomo.clustering_Nbin; nz++) {
      for (int i = 0; i < NTHETA; i++) {
        w_vec[nz * like.Ntheta + i] = 0;
        for (int l = 1; l < LMAX; l++) {
          w_vec[nz * like.Ntheta + i] += Pl[i][l] * Cl[nz][l];
        }
      }
    }
    for (int nz = 0; nz < tomo.clustering_Nbin; nz++) {
      free(Cl[nz]);
    }
    free(Cl);
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }

  return w_vec[ni * like.Ntheta + nt];
}
