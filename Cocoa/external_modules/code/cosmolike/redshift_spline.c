#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_spline.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "basics.h"
#include "cosmo3D.h"
#include "recompute.h"
#include "redshift_spline.h"
#include "structs.h"

#include "log.c/src/log.h"

//#define Z_SPLINE_TYPE gsl_interp_akima
#define Z_SPLINE_TYPE gsl_interp_cspline
// ------------------------------------------------------------------------
// integration boundary routines
// ------------------------------------------------------------------------

double amin_source(int i) {
  if (i == -1 || redshift.shear_photoz == 1 || redshift.shear_photoz == 2 ||
      redshift.shear_photoz == 4) {
    return 1. / (redshift.shear_zdistrpar_zmax + 1.);
  }
  if (redshift.shear_photoz == 0) {
    return 1. / (1 + tomo.shear_zmax[i]);
  }
  return 1. /
         (1 + fmin(tomo.shear_zmax[i] + 5. * nuisance.sigma_zphot_shear[i] +
                       fabs(nuisance.bias_zphot_shear[i]),
                   redshift.shear_zdistrpar_zmax));
}

double amax_source(int i __attribute__((unused))) {
  return 1. / (1. + fmax(redshift.shear_zdistrpar_zmin, 0.001));
}

double amax_source_IA(int i) {
  if (i == -1 || redshift.shear_photoz == 1 || redshift.shear_photoz == 2 ||
      redshift.shear_photoz == 4) {
    return 1. / (1. + fmax(redshift.shear_zdistrpar_zmin, 0.001));
  }
  if (redshift.shear_photoz == 0) {
    return 1. / (1. + fmax(tomo.shear_zmin[i], 0.001));
  }
  return 1. /
         (1 + fmax(tomo.shear_zmin[i] - 5. * nuisance.sigma_zphot_shear[i] -
                       fabs(nuisance.bias_zphot_shear[i]),
                   0.001));
}

double amin_lens(int i) {
  if (i == -1 || redshift.clustering_photoz == 1 ||
      redshift.clustering_photoz == 2) {
    return 1. / (redshift.clustering_zdistrpar_zmax + 1.);
  }
  if (redshift.clustering_photoz == 0) {
    return 1. / (1 + tomo.clustering_zmax[i]);
  }
  if (redshift.clustering_photoz == 4) {
    return 1. / (1 + tomo.clustering_zmax[i] +
                 2. * fabs(nuisance.bias_zphot_clustering[i]));
  }
  return 1. / (1 + fmin(tomo.clustering_zmax[i] +
                            5. * nuisance.sigma_zphot_clustering[i] +
                            fabs(nuisance.bias_zphot_clustering[i]),
                        redshift.clustering_zdistrpar_zmax));
}

double amax_lens(int i) {
  if (gbias.b_mag[i] != 0) {
    return 1. / (1. + fmax(redshift.shear_zdistrpar_zmin, 0.001));
  }
  if (i == -1 || redshift.clustering_photoz == 1 ||
    redshift.clustering_photoz == 2
  ) {
    return 1. / (1. + fmax(redshift.clustering_zdistrpar_zmin, 0.001));
  }
  if (redshift.clustering_photoz == 0) {
    return 1. / (1. + fmax(tomo.clustering_zmin[i], 0.001));
  }
  if (redshift.clustering_photoz == 4) {
    return 1. / (1 + fmax(tomo.clustering_zmin[i] -
      2. * fabs(nuisance.bias_zphot_clustering[i]), 0.001));
  }
  return 1. / (1 + fmax(tomo.clustering_zmin[i] -
    5. * nuisance.sigma_zphot_clustering[i] -
    fabs(nuisance.bias_zphot_clustering[i]), 0.001));
}

// ------------------------------------------------------------------------
//redshift overlap tests, allowed tomography combinations
// ------------------------------------------------------------------------

// test whether the (l,zl) bin is in the linear clustering regime
// return 1 if true, 0 otherwise
int test_kmax(double l, int zl)
{
  static double chiref[10] = {-1.};
  if (chiref[0] < 0)
  {
    for (int i = 0; i < tomo.clustering_Nbin; i++)
    {
      chiref[i] = chi(1. /
              (1. + 0.5 * (tomo.clustering_zmin[i] + tomo.clustering_zmax[i])));
    }
  }
  const double R_min = like.Rmin_bias; // set minimum scale to which
                                       // we trust our bias model, in Mpc/h
  const double kmax = constants.twopi / R_min * cosmology.coverH0;
  if ((l + 0.5) / chiref[zl] < kmax)
  {
    return 1;
  }
  return 0;
}

// test whether source bin zs is behind lens bin zl
int test_zoverlap(int zl, int zs) {
  // if (tomo.clustering_zmin[zl] >= tomo.shear_zmax[zs]) {return 0;}
  if (ggl_efficiency(zl, zs) > survey.ggl_overlap_cut) {

    return 1;
  }
  if (redshift.shear_photoz < 4 &&
      tomo.clustering_zmax[zl] <= tomo.shear_zmin[zs]) {
    return 1;
  }
  if (redshift.shear_photoz == 4 && redshift.clustering_photoz != 4 &&
      tomo.clustering_zmax[zl] < zmean_source(zs)) {
    return 1;
  }
  if (redshift.shear_photoz == 4 && redshift.clustering_photoz == 4 &&
      zmean(zl) + 0.1 < zmean_source(zs)) {
    return 1;
  }
  return 0;
}

// test whether source bin zs is behind lens bin zl
int test_zoverlap_cov(int zl, int zs) {
  // if (tomo.clustering_zmin[zl] >= tomo.shear_zmax[zs]) {return 0;}
  if (ggl_efficiency(zl, zs) > survey.ggl_overlap_cut) {
    return 1;
  }
  if (redshift.shear_photoz < 4 &&
      tomo.clustering_zmax[zl] <= tomo.shear_zmin[zs]) {
    return 1;
  }
  if (redshift.shear_photoz == 4 && redshift.clustering_photoz != 4 &&
      tomo.clustering_zmax[zl] < zmean_source(zs)) {
    return 1;
  }
  if (redshift.shear_photoz == 4 && redshift.clustering_photoz == 4 &&
      zmean(zl) + 0.1 < zmean_source(zs)) {
    return 1;
  }
  return 0;
}

// test whether source bin zs is behind lens bin zl
int test_zoverlap_c(int zc, int zs) {
  if (redshift.shear_photoz < 4 &&
      tomo.cluster_zmax[zc] <= tomo.shear_zmin[zs]) {
    return 1;
  }
  if (redshift.shear_photoz == 4 && tomo.cluster_zmax[zc] < zmean_source(zs)) {
    return 1;
  }
  return 0;
}

int N_ggl(int zl, int zs) {
  static int N[10][10] = {{-42}};
  if (N[0][0] < 0) {
    int i, j, n;
    n = 0;
    for (i = 0; i < tomo.clustering_Nbin; i++) {
      for (j = 0; j < tomo.shear_Nbin; j++) {
        if (test_zoverlap(i, j)) {
          N[i][j] = n;
          n++;
        } else {
          N[i][j] = -1;
        }
      }
    }
  }
  return N[zl][zs];
}

void write_gglensing_zbins(char *surveyname) {
  FILE *F1;
  char filename[400];
  sprintf(filename, "%s.gglensing_zbins", surveyname);
  printf("filename gglensing %s\n", filename);
  F1 = fopen(filename, "w");
  fprintf(F1, "#zl zs  accept?  overlap metric\n");
  fprintf(F1, "#pipeline supplied ggl_overlap_cut = %.3f\n",
          survey.ggl_overlap_cut);
  for (int i = 0; i < tomo.clustering_Nbin; i++) {
    for (int j = 0; j < tomo.shear_Nbin; j++) {
      fprintf(F1, "%d   %d   %d    %.3f\n", i, j, test_zoverlap(i, j),
              ggl_efficiency(i, j));
    }
  }
  fclose(F1);
}

int ZL(int Nbin) {
  static int N[100] = {-42};
  if (N[0] < -1) {
    int i, j, n = 0;
    for (i = 0; i < tomo.clustering_Nbin; i++) {
      for (j = 0; j < tomo.shear_Nbin; j++) {
        if (test_zoverlap(i, j)) {
          N[n] = i;
          n++;
        }
      }
    }
  }
  return N[Nbin];
}

int ZS(int Nbin) {
  static int N[100] = {-42};
  if (N[0] < -1) {
    int i, j, n = 0;
    for (i = 0; i < tomo.clustering_Nbin; i++) {
      for (j = 0; j < tomo.shear_Nbin; j++) {
        if (test_zoverlap(i, j)) {
          N[n] = j;
          n++;
        }
      }
    }
  }
  return N[Nbin];
}

int N_cgl(int zc, int zs) {
  static int N[10][10] = {{-42}};
  if (N[0][0] < 0) {
    int i, j, n;
    n = 0;
    for (i = 0; i < tomo.cluster_Nbin; i++) {
      for (j = 0; j < tomo.shear_Nbin; j++) {
        if (test_zoverlap_c(i, j)) {
          N[i][j] = n;
          n++;
        } else {
          N[i][j] = -1;
        }
      }
    }
  }
  return N[zc][zs];
}

int ZC(int Nbin) {
  static int N[100] = {-42};
  if (N[0] < -1) {
    int i, j, n = 0;
    for (i = 0; i < tomo.cluster_Nbin; i++) {
      for (j = 0; j < tomo.shear_Nbin; j++) {
        if (test_zoverlap_c(i, j)) {
          N[n] = i;
          n++;
        }
      }
    }
  }
  return N[Nbin];
}

int ZSC(int Nbin) {
  static int N[100] = {-42};
  if (N[0] < -1) {
    int i, j, n = 0;
    for (i = 0; i < tomo.cluster_Nbin; i++) {
      for (j = 0; j < tomo.shear_Nbin; j++) {
        if (test_zoverlap_c(i, j)) {
          N[n] = j;
          n++;
        }
      }
    }
  }
  return N[Nbin];
}

// find shear tomography bin number N_shear of tomography combination (z1,z2)
int N_shear(int z1, int z2) {
  static int N[10][10] = {{-42}};
  if (N[0][0] < -1) {
    int i, j, n = 0;
    for (i = 0; i < tomo.shear_Nbin; i++) {
      for (j = i; j < tomo.shear_Nbin; j++) {
        N[i][j] = n;
        N[j][i] = n;
        n++;
      }
    }
  }
  return N[z1][z2];
}

// find z1 of tomography combination (z1,z2) constituting
// shear tomography bin Nbin
int Z1(int Nbin) {
  static int N[55] = {-42};
  if (N[0] < -1) {
    int i, j, n = 0;
    for (i = 0; i < tomo.shear_Nbin; i++) {
      for (j = i; j < tomo.shear_Nbin; j++) {
        N[n] = i;
        n++;
      }
    }
  }
  return N[Nbin];
}

// find z2 of tomography combination (z1,z2) constituting
// shear tomography bin Nbin
int Z2(int Nbin) {
  static int N[55] = {-42};
  if (N[0] < -1) {
    int i, j, n = 0;
    for (i = 0; i < tomo.shear_Nbin; i++) {
      for (j = i; j < tomo.shear_Nbin; j++) {
        N[n] = j;
        n++;
      }
    }
  }
  return N[Nbin];
}

// find zcl1 of tomography combination (zcl1,zcl2)
// constituting galaxy clustering tomography bin Nbin
int Zcl1(int Nbin) {
  static int N[55] = {-42};
  if (N[0] < -1) {
    int i, j, n = 0;
    for (i = 0; i < tomo.clustering_Nbin; i++) {
      for (j = i; j < tomo.clustering_Nbin; j++) {
        N[n] = i;
        n++;
      }
    }
  }
  return N[Nbin];
}

// find zcl2 of tomography combination (zcl1,zcl2)
// constituting galaxy clustering tomography bin Nbin
int Zcl2(int Nbin) {
  static int N[55] = {-42};
  if (N[0] < -1) {
    int i, j, n = 0;
    for (i = 0; i < tomo.clustering_Nbin; i++) {
      for (j = i; j < tomo.clustering_Nbin; j++) {
        N[n] = j;
        n++;
      }
    }
  }
  return N[Nbin];
}

// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// routines for redshift distributions, including photo-zs * (optional)
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------

// ------------------------------------------------------------------------
// Start shear routines
// ------------------------------------------------------------------------

double sigma_zphot_shear(double z, int nz) {
  return (1. + z) * nuisance.sigma_zphot_shear[nz];
}
double bias_zphot_shear(double z, int nz) {
  return (1. + z) * nuisance.bias_zphot_shear[nz];
}

// return nz(z,j) based on redshift file with structure z[i]
// nz[0][i] .. nz[tomo.shear_Nbin-1][i]
double zdistr_histo_n(double z, void *params)
{
  double *array = (double *)params;
  static double **tab;
  FILE *ein;
  static double zhisto_max, zhisto_min, dz;

  if (tab == 0) {
    double *z_v;
    int i, zbins;
    zbins = line_count(redshift.shear_REDSHIFT_FILE);
    tab = create_double_matrix(0, tomo.shear_Nbin - 1, 0, zbins - 1);
    z_v = create_double_vector(0, zbins - 1);
    ein = fopen(redshift.shear_REDSHIFT_FILE, "r");
    for (i = 0; i < zbins; i++) {
      fscanf(ein, "%le", &z_v[i]);
      if (i > 0 && z_v[i] < z_v[i - 1]) {
        break;
      }
      for (int k = 0; k < tomo.shear_Nbin; k++) {
        fscanf(ein, "%le", &tab[k][i]);
      }
    }
    fclose(ein);
    dz = (z_v[i - 1] - z_v[0]) / (1. * i - 1.);
    zhisto_max = z_v[i - 1] + dz;
    zhisto_min = z_v[0];
    zbins = i;
    // now, set tomography bin boundaries
    for (int k = 0; k < tomo.shear_Nbin; k++) {
      double max = tab[k][0];
      for (i = 1; i < zbins; i++) {
        if (tab[k][i] > max) {
          max = tab[k][i];
        }
      }
      i = 0;
      while (tab[k][i] < 1.e-8 * max && i < zbins - 2) {
        i++;
      }
      tomo.shear_zmin[k] = z_v[i];
      i = zbins - 1;
      while (tab[k][i] < 1.e-8 * max && i > 0) {
        i--;
      }
      tomo.shear_zmax[k] = z_v[i];
      log_info("tomo.shear_zmin[%d] = %.3f,tomo.shear_zmax[%d] = %.3f", k,
             tomo.shear_zmin[k], k, tomo.shear_zmax[k]);
    }
    free_double_vector(z_v, 0, zbins - 1);
    if (zhisto_max < tomo.shear_zmax[tomo.shear_Nbin - 1] || zhisto_min > tomo.shear_zmin[0]) {
      log_fatal("zhisto_min = %e,zhisto_max = %e", zhisto_min, zhisto_max);
      log_fatal("tomo.shear_zmin[0] = %e, tomo.shear_zmax[N-1] = %e", tomo.shear_zmin[0], tomo.shear_zmax[tomo.shear_Nbin - 1]);
      log_fatal("%s parameters incompatible with tomo.shear bin choice",redshift.shear_REDSHIFT_FILE);
      exit(1);
    }
  }

  if ((z >= zhisto_min) && (z < zhisto_max)) {
    return tab[(int)array[0]][(int)floor((z - zhisto_min) / dz)];
  }
  return 0.;
}

// return nz(z) based on redshift file with one redshift distribution
double zdistr_histo_1(double z, void *params __attribute__((unused))) {
  static double *tab = 0;
  FILE *ein;
  static double zhisto_max, zhisto_min, dz;
  if (tab == 0) {
    double *z_v, space1, space2;
    int i, zbins;
    zbins = line_count(redshift.shear_REDSHIFT_FILE);
    tab = create_double_vector(0, zbins - 1);
    z_v = create_double_vector(0, zbins - 1);
    ein = fopen(redshift.shear_REDSHIFT_FILE, "r");
    for (i = 0; i < zbins; i++) {
      fscanf(ein, "%le %le %le %le\n", &z_v[i], &space1, &space2, &tab[i]);
      if (i > 0 && z_v[i] < z_v[i - 1]) {
        break;
      }
    }
    fclose(ein);
    dz = (z_v[i - 1] - z_v[0]) / (1. * i - 1.);
    zhisto_max = z_v[i - 1] + dz;
    zhisto_min = z_v[0];
    redshift.shear_zdistrpar_zmin = zhisto_min;
    redshift.shear_zdistrpar_zmax = zhisto_max;
    free_double_vector(z_v, 0, zbins - 1);
  }
  if ((z >= zhisto_min) && (z < zhisto_max)) {
    return tab[(int)floor((z - zhisto_min) / dz)];
  }
  return 0.;
}

double n_of_z(double z, int nz) {
  double sigma = 0.15;
  double x = 1.0 + 0.5 * nz - z;
  return 1. / sqrt(2. * M_PI * sigma * sigma) *
         exp(-x * x / (2. * sigma * sigma));
}

// returns n(ztrue | j), works only with binned
// distributions; j>= 0 -> tomography bin j
double zdistr_photoz(double zz, int j) {
  static double **table = 0;
  static double *z_v = 0;
  static double da = 0.0;
  static double zhisto_max, zhisto_min;
  static nuisancepara N;
  static int zbins = -1;
  static gsl_spline *photoz_splines[11];
  //static gsl_interp_accel *photoz_accel[11];

  if (redshift.shear_photoz == -1) {
    return n_of_z(zz, j);
  }
  if ((redshift.shear_photoz != 4 && recompute_zphot_shear(N)) || table == 0) {
    update_nuisance(&N);
    if (table == 0) {
      int zbins1 = line_count(redshift.shear_REDSHIFT_FILE);
      if (redshift.shear_photoz != 4) {
        zbins = zbins1 * 20;
      } else {
        zbins = zbins1;
      }
      table = create_double_matrix(0, tomo.shear_Nbin, 0, zbins - 1);
      z_v = create_double_vector(0, zbins - 1);
      for (int i = 0; i < tomo.shear_Nbin + 1; i++) {
        photoz_splines[i] = gsl_spline_alloc(Z_SPLINE_TYPE, zbins);
        //photoz_accel[i] = gsl_interp_accel_alloc();
      }

      if (redshift.shear_photoz == 4) { // if multihisto, force zmin, zmax, tomo
                                        // bins to match supplied file
        FILE *ein;
        double space;
        int i, k;
        ein = fopen(redshift.shear_REDSHIFT_FILE, "r");
        for (i = 0; i < zbins1; i++) {
          fscanf(ein, "%le", &z_v[i]);
          if (i > 0 && z_v[i] < z_v[i - 1]) {
            break;
          }
          for (k = 0; k < tomo.shear_Nbin; k++) {
            fscanf(ein, "%le", &space);
          }
        }
        fclose(ein);
        redshift.shear_zdistrpar_zmin = fmax(z_v[0], 1.e-5);
        redshift.shear_zdistrpar_zmax =
            z_v[i - 1] + (z_v[i - 1] - z_v[0]) / (zbins1 - 1.);
      }
    }

    zhisto_max = redshift.shear_zdistrpar_zmax;
    zhisto_min = redshift.shear_zdistrpar_zmin;
    da = (zhisto_max - zhisto_min) / (1.0 * zbins);
    for (int i = 0; i < zbins; i++) {
      z_v[i] = zhisto_min + (i + 0.5) * da;
    }
    double array[4], NORM[11], norm, x1, x2, eta, outfrac;
    // the outlier fraction (outfrac) should be specified externally. This is a
    // temporary hack.
    for (int i = 0; i < tomo.shear_Nbin; i++) {
      array[0] = tomo.shear_zmin[i];
      array[1] = tomo.shear_zmax[i];
      norm = 0.0;
      outfrac = 0.05;
      // eta = 0.68; //this is the approximate mixing parameter for the pseudo
      // voigt for equal FWHM L and G distributions. Note than when using this
      // formulation, the non-equivalent profile re-scalings must be accounted
      // for below.
      eta = 0.736; // this is the approximate mixing parameter for the pseudo
                   // voigt for equally rescaled L and G distributions. Note
                   // that this leads to a slightly broader profile than when
                   // using equal FWHM L and G distributions (at a given
                   // gaussian sigma).
      // Currently, the input sigma parameter has not been rescaled. It still
      // gives the width for the Gaussian component of the psuedo Voigt profile.
      // In general, for a given FWHM, the corresponding sigma is smaller for
      // the Voigt profile. For equal FWHM in the component profiles, eta = 0.68
      // and the effective FWHM is 1.64 that of the Guassian case. For
      // equivalent rescalings, [which is the currently implemented case]
      // eta=0.736 and the effective FWHM is 1.79 that of the Gaussian case.
      switch (redshift.shear_photoz) {
      case 0: // no photo-zs, split 'true' n(z) histogram in bins
        norm = int_gsl_integrate_medium_precision(
            zdistr_histo_1, (void *) array, array[0], array[1], NULL, 1024);
        for (int k = 0; k < zbins; k++) {
          table[i + 1][k] = 0;
          if (z_v[k] >= array[0] && z_v[k] <= array[1]) {
            table[i + 1][k] = zdistr_histo_1(z_v[k], (void *)array) / norm;
          }
        }
        break;
      case 1: // Pseudo Voigt (Lorentzian + Gaussian) photo-zs
        norm = 0.;
        if (sigma_zphot_shear(0., i) == 0.) {
          log_fatal("Source galaxy photo-z model "
                 "underdetermined!\nredshift.shear_photoz =1, but "
                 "nuisance.sigma_zphot_shear[%d] not set",
                 i);
          exit(1);
        }
        for (int k = 0; k < zbins; k++) {
          x1 = (array[0] - z_v[k] + bias_zphot_shear(z_v[k], i)) /
               (sqrt(2.) * sigma_zphot_shear(z_v[k], i));
          x2 = (array[1] - z_v[k] + bias_zphot_shear(z_v[k], i)) /
               (sqrt(2.) * sigma_zphot_shear(z_v[k], i));
          table[i + 1][k] =
              zdistr_histo_1(z_v[k], (void *)array) *
              (eta * 0.3183 * (atan(x2) - atan(x1)) +
               ((1.0 - eta) * 0.5 * (gsl_sf_erf(x2) - gsl_sf_erf(x1))));
          // this creates a pseudo Voigt profile by adding a Gaussian and
          // Lorentian (convolved with a tophat) with the correct weights. See,
          // eg, the Wikipedia article on Voigt profiles for an explanation of
          // where the numbers come from.
          norm += table[i + 1][k] * da;
        }
        for (int k = 0; k < zbins; k++) {
          table[i + 1][k] /= norm;
        }
        break;
      case 2: // Pseudo Voigt (Lorentzian + Gaussian) + outlier (currently flat)
              // photo-zs
        norm = 0.;

        if (sigma_zphot_shear(0., i) == 0.) {
          log_fatal("Source galaxy photo-z model "
                 "underdetermined!\nredshift.shear_photoz =2, but "
                 "nuisance.sigma_zphot_shear[%d] not set",i);
          exit(1);
        }
        for (int k = 0; k < zbins; k++) {
          x1 = (array[0] - z_v[k] + bias_zphot_shear(z_v[k], i)) /
               (sqrt(2.) * sigma_zphot_shear(z_v[k], i));
          x2 = (array[1] - z_v[k] + bias_zphot_shear(z_v[k], i)) /
               (sqrt(2.) * sigma_zphot_shear(z_v[k], i));
          // table[i+1][k] =
          // zdistr_histo_1(zi,(void*)array)*(1.0-outfrac)*(eta*0.3183*(atan(x2)-atan(x1))
          // + ((1.0-eta)*0.5*(gsl_sf_erf(x2)-gsl_sf_erf(x1)))) + outfrac; the
          // outliers are currently assumed to have a flat distribution from
          // zmin to zmax. I think the above implementation is not correctly
          // normalized for the outliers. Instead do it at the normalization
          // step. This should probably be updated to allow for a leakage
          // matrix.
          table[i + 1][k] =
              zdistr_histo_1(z_v[k], (void *)array) *
              (eta * 0.3183 * (atan(x2) - atan(x1)) +
               ((1.0 - eta) * 0.5 * (gsl_sf_erf(x2) - gsl_sf_erf(x1))));
          norm += table[i + 1][k] * da;
        }
        for (int k = 0; k < zbins; k++) {
          table[i + 1][k] =
              outfrac / da + (1. - outfrac) * table[i + 1][k] / norm;
        }
        break;
      case 3: // Gaussian photo-zs
        norm = 0.;
        if (sigma_zphot_shear(0., i) == 0.) {
          log_fatal("Source galaxy photo-z model "
                 "underdetermined!\nredshift.shear_photoz =3, but "
                 "nuisance.sigma_zphot_shear[%d] not set",i);
          exit(1);
        }
        for (int k = 0; k < zbins; k++) {
          x1 = (array[0] - z_v[k] + bias_zphot_shear(z_v[k], i)) /
               (sqrt(2.) * sigma_zphot_shear(z_v[k], i));
          x2 = (array[1] - z_v[k] + bias_zphot_shear(z_v[k], i)) /
               (sqrt(2.) * sigma_zphot_shear(z_v[k], i));
          table[i + 1][k] = 0.5 * zdistr_histo_1(z_v[k], (void *)array) *
                            (gsl_sf_erf(x2) - gsl_sf_erf(x1));
          norm += table[i + 1][k] * da;
        }
        for (int k = 0; k < zbins; k++) {
          table[i + 1][k] /= norm;
        }
        break;
      case 4: // histogram file contains n(z) estimates for each bin
        array[0] = 1.0 * i;
        // call zdistr_histo_n once to update bin boundaries
        zdistr_histo_n(0., (void *)array);
        norm = 0.0;
        for (int k = 0; k < zbins; k++) {
          table[i + 1][k] = zdistr_histo_n(z_v[k], (void *)array); // norm;
          norm += table[i + 1][k] * da;
        }
        for (int k = 0; k < zbins; k++) {
          table[i + 1][k] /= norm;
        }
        break;
      default:
        log_fatal("redshift.shear_photoz = %d not supported in this cosmolike "
               "version",redshift.shear_photoz);
        exit(1);
      }
      NORM[i] = norm;
    }
    // calculate normalized overall redshift distribution (without bins), store
    // in table[0][:]
    norm = 0;
    for (int i = 0; i < tomo.shear_Nbin; i++) {
      norm += NORM[i];
    }

    for (int k = 0; k < zbins; k++) {
      table[0][k] = 0;
      for (int i = 0; i < tomo.shear_Nbin; i++) {
        table[0][k] += table[i + 1][k] * NORM[i] / norm;
      }
    }
    for (int i = -1; i < tomo.shear_Nbin; i++) {
      gsl_spline_init(photoz_splines[i + 1], z_v, table[i + 1], zbins);
    }
  }
  if (j >= tomo.shear_Nbin) {
    log_fatal("zdistr_photoz(z,%d) outside tomo.shear_Nbin range",j);
    exit(1);
  }
  if (redshift.shear_photoz == 4) {
    zz = zz - nuisance.bias_zphot_shear[j];
  }
  if (zz <= z_v[0] || zz >= z_v[zbins - 1])
    return 0.0;
  // VM - COBAYA_SAMPLER - SET ACCEL TO NULL (THREAD SAFE)
  return gsl_spline_eval(photoz_splines[j + 1], zz, NULL);
}

// ------------------------------------------------------------------------
// End shear routines
// ------------------------------------------------------------------------

double sigma_zphot_clustering(double z, int nz) {
  return (1. + z) * nuisance.sigma_zphot_clustering[nz];
}

double bias_zphot_clustering(double z, int nz) {
  return (1. + z) * nuisance.bias_zphot_clustering[nz];
}

// return pf(z) based on redshift file with one redshift distribution
double pf_histo(double z, void *params __attribute__((unused))) {
  static double *tab = 0;
  FILE *ein;

  static double zhisto_max, zhisto_min, dz;

  if (tab == 0) {
    double *z_v, space1, space2;
    int i, zbins;
    zbins = line_count(redshift.clustering_REDSHIFT_FILE);
    tab = create_double_vector(0, zbins - 1);
    z_v = create_double_vector(0, zbins - 1);
    ein = fopen(redshift.clustering_REDSHIFT_FILE, "r");
    for (i = 0; i < zbins; i++) {
      fscanf(ein, "%le %le %le %le\n", &z_v[i], &space1, &space2, &tab[i]);
      if (i > 0 && z_v[i] < z_v[i - 1]) {
        break;
      }
    }
    fclose(ein);
    dz = (z_v[i - 1] - z_v[0]) / (1. * i - 1.);
    zhisto_max = z_v[i - 1] + dz;
    zhisto_min = z_v[0];
    // redshift.clustering_zdistrpar_zmin = zhisto_min;
    // redshift.clustering_zdistrpar_zmax = zhisto_max;
    free_double_vector(z_v, 0, zbins - 1);
    if (zhisto_max < tomo.clustering_zmax[tomo.clustering_Nbin - 1] ||
        zhisto_min > tomo.clustering_zmin[0]) {
      log_fatal("Error in redshift_spline.c:pf_histo.c: %s parameters "
             "incompatible with tomo.clustering bin choice", redshift.clustering_REDSHIFT_FILE);
      exit(1);
    }
  }

  if ((z >= zhisto_min) && (z < zhisto_max)) {
    return tab[(int)floor((z - zhisto_min) / dz)];
  }
  return 0.;
}

// return pf(z,j) based on redshift file with structure
// z[i] nz[0][i] .. nz[tomo.clustering_Nbin-1][i]
double pf_histo_n(double z, void *params) {
  double *array = (double *)params;
  static double **tab;
  FILE *ein;
  static double zhisto_max, zhisto_min, dz;

  if (tab == 0) {
    double *z_v;
    int i, k, zbins;
    zbins = line_count(redshift.clustering_REDSHIFT_FILE);
    tab = create_double_matrix(0, tomo.clustering_Nbin - 1, 0, zbins - 1);
    z_v = create_double_vector(0, zbins - 1);
    ein = fopen(redshift.clustering_REDSHIFT_FILE, "r");
    for (i = 0; i < zbins; i++) {
      fscanf(ein, "%le", &z_v[i]);
      if (i > 0 && z_v[i] < z_v[i - 1]) {
        break;
      }
      for (k = 0; k < tomo.clustering_Nbin; k++) {
        fscanf(ein, " %le", &tab[k][i]);
      }
    }
    fclose(ein);
    dz = (z_v[i - 1] - z_v[0]) / (1. * i - 1.);
    zhisto_max = z_v[i - 1] + dz;
    zhisto_min = z_v[0];
    // now, set tomography bin boundaries
    for (k = 0; k < tomo.clustering_Nbin; k++) {
      double max = tab[k][0];
      for (i = 1; i < zbins; i++) {
        if (tab[k][i] > max) {
          max = tab[k][i];
        }
      }
      i = 0;
      while (tab[k][i] < 1.e-8 * max && i < zbins - 2) {
        i++;
      }
      tomo.clustering_zmin[k] = z_v[i];
      i = zbins - 1;
      while (tab[k][i] < 1.e-8 * max && i > 0) {
        i--;
      }
      tomo.clustering_zmax[k] = z_v[i];
      log_info("tomo.clustering_zmin[%d] = %.3f,tomo.clustering_zmax[%d] = %.3f",
          k, tomo.clustering_zmin[k], k, tomo.clustering_zmax[k]);
    }
    free_double_vector(z_v, 0, zbins - 1);
    if (zhisto_max < tomo.clustering_zmax[tomo.clustering_Nbin - 1] ||
        zhisto_min > tomo.clustering_zmin[0]) {
      log_fatal("%e %e   %e %e", zhisto_min, tomo.clustering_zmin[0], zhisto_max, tomo.clustering_zmax[tomo.clustering_Nbin - 1]);
      log_fatal("%s parameters incompatible with tomo.clustering bin choice", redshift.clustering_REDSHIFT_FILE);
      exit(1);
    }
  }

  if ((z >= zhisto_min) && (z < zhisto_max)) {
    return tab[(int)array[0]][(int)floor((z - zhisto_min) / dz)];
  }
  return 0.;
}

// returns n(ztrue, j), works only with binned distributions; j
// =-1 -> no tomography; j>= 0 -> tomography bin j
double pf_photoz(double zz, int j) {
  static double **table = 0;
  static double *z_v = 0;
  static double da = 0.0;
  static double zhisto_max, zhisto_min;
  static nuisancepara N;
  static int zbins = -1;
  static gsl_spline *photoz_splines[11];
  //static gsl_interp_accel *photoz_accel[11];

  if (redshift.clustering_photoz == -1) {
    return n_of_z(zz, j);
  }
  if ((redshift.clustering_photoz != 4 && recompute_zphot_clustering(N)) ||
      table == 0) {
    update_nuisance(&N);
    if (table == 0) {
      zbins = line_count(redshift.clustering_REDSHIFT_FILE);
      if (redshift.clustering_photoz != 4 && redshift.clustering_photoz != 0) {
        pf_histo(0.5, NULL);
        zbins *= 20;
      } // upsample if convolving with analytic photo-z model
      table = create_double_matrix(0, tomo.clustering_Nbin, 0, zbins - 1);
      z_v = create_double_vector(0, zbins - 1);
      for (int i = 0; i < tomo.clustering_Nbin + 1; i++) {
        photoz_splines[i] = gsl_spline_alloc(Z_SPLINE_TYPE, zbins);
        //photoz_accel[i] = gsl_interp_accel_alloc();
      }
      if (redshift.clustering_photoz == 4) { // if multihisto, force zmin, zmax,
                                             // tomo bins to match supplied file
        FILE *ein;
        double space;
        int i, k;
        ein = fopen(redshift.clustering_REDSHIFT_FILE, "r");
        for (i = 0; i < zbins; i++) {
          fscanf(ein, "%le", &z_v[i]);
          if (i > 0 && z_v[i] < z_v[i - 1]) {
            break;
          }
          for (k = 0; k < tomo.clustering_Nbin; k++) {
            fscanf(ein, "%le", &space);
          }
        }
        fclose(ein);
        redshift.clustering_zdistrpar_zmin = fmax(z_v[0], 1.e-5);
        redshift.clustering_zdistrpar_zmax =
            z_v[i - 1] + (z_v[i - 1] - z_v[0]) / (zbins - 1.);
      }
    }
    zhisto_max = redshift.clustering_zdistrpar_zmax;
    zhisto_min = redshift.clustering_zdistrpar_zmin;
    da = (zhisto_max - zhisto_min) / (1.0 * zbins);
    for (int i = 0; i < zbins; i++) {
      z_v[i] = zhisto_min + (i + 0.5) * da;
    }

    double array[4], NORM[11], norm, x1, x2, eta, outfrac, zi;
    // the outlier fraction (outfrac) should be specified externally. This is a
    // temporary hack.
    int i, k;

    for (i = 0; i < tomo.clustering_Nbin; i++) {
      array[0] = tomo.clustering_zmin[i];
      array[1] = tomo.clustering_zmax[i];
      norm = 0.0;
      outfrac = 0.05;
      // eta = 0.68; //this is the approximate mixing parameter for the pseudo
      // voigt for equal FWHM L and G distributions. Note than when using this
      // formulation, the non-equivalent profile re-scalings must be accounted
      // for below.
      eta = 0.736; // this is the approximate mixing parameter for the pseudo
                   // voigt for equally rescaled L and G distributions. Note
                   // that this leads to a slightly broader profile than when
                   // using equal FWHM L and G distributions (at a given
                   // gaussian sigma).
      // Currently, the input sigma parameter has not been rescaled. It still
      // gives the width for the Gaussian component of the psuedo Voigt profile.
      // In general, for a given FWHM, the corresponding sigma is smaller for
      // the Voigt profile. For equal FWHM in the component profiles, eta = 0.68
      // and the effective FWHM is 1.64 that of the Guassian case. For
      // equivalent rescalings, [which is the currently implemented case]
      // eta=0.736 and the effective FWHM is 1.79 that of the Gaussian case.
      switch (redshift.clustering_photoz) {
      case 0: // no photo-zs, split 'true' n(z) histogram in bins
        norm = int_gsl_integrate_medium_precision(
            pf_histo, (void *)array, array[0], array[1], NULL, 1024);
        for (k = 0; k < zbins; k++) {
          zi = z_v[k];
          table[i + 1][k] = 0.;
          if (zi >= array[0] && zi <= array[1]) {
            table[i + 1][k] = pf_histo(zi, (void *)array) / norm;
          }
        }
        break;
      case 1: // Pseudo Voigt (Lorentzian + Gaussian) photo-zs
        norm = 0.;
        if (sigma_zphot_clustering(0., i) == 0.) {
          log_fatal("Lens galaxy photo-z model "
                 "underdetermined!\nredshift.clustering_photoz =1, but "
                 "nuisance.sigma_zphot_clustering[%d] not set",i);
          exit(1);
        }
        for (k = 0; k < zbins; k++) {
          zi = z_v[k];
          x1 = (array[0] - zi + bias_zphot_clustering(zi, i)) /
               (sqrt(2.) * sigma_zphot_clustering(zi, i));
          x2 = (array[1] - zi + bias_zphot_clustering(zi, i)) /
               (sqrt(2.) * sigma_zphot_clustering(zi, i));
          table[i + 1][k] =
              pf_histo(zi, (void *)array) *
              (eta * 0.3183 * (atan(x2) - atan(x1)) +
               ((1.0 - eta) * 0.5 * (gsl_sf_erf(x2) - gsl_sf_erf(x1))));
          // this creates a pseudo Voigt profile by adding a Gaussian and
          // Lorentian (convolved with a tophat) with the correct weights. See,
          // eg, the Wikipedia article on Voigt profiles for an explanation of
          // where the numbers come from.
          norm += table[i + 1][k] * da;
        }
        for (k = 0; k < zbins; k++) {
          table[i + 1][k] /= norm;
        }
        break;
      case 2: // Pseudo Voigt (Lorentzian + Gaussian) + outlier (currently flat)
              // photo-zs
        norm = 0.;

        if (sigma_zphot_clustering(0., i) == 0.) {
          log_fatal("Lens galaxy photo-z model "
                 "underdetermined!\nredshift.clustering_photoz =2, but "
                 "nuisance.sigma_zphot_clustering[%d] not set",i);
          exit(1);
        }
        for (k = 0; k < zbins; k++) {
          zi = z_v[k];
          x1 = (array[0] - zi + bias_zphot_clustering(zi, i)) /
               (sqrt(2.) * sigma_zphot_clustering(zi, i));
          x2 = (array[1] - zi + bias_zphot_clustering(zi, i)) /
               (sqrt(2.) * sigma_zphot_clustering(zi, i));
          // table[i+1][k] =
          // pf_histo(zi,(void*)array)*(1.0-outfrac)*(eta*0.3183*(atan(x2)-atan(x1))
          // + ((1.0-eta)*0.5*(gsl_sf_erf(x2)-gsl_sf_erf(x1)))) + outfrac; the
          // outliers are currently assumed to have a flat distribution from
          // zmin to zmax. I think the above implementation is not correctly
          // normalized for the outliers. Instead do it at the normalization
          // step. This should probably be updated to allow for a leakage
          // matrix.
          table[i + 1][k] =
              pf_histo(zi, (void *)array) *
              (eta * 0.3183 * (atan(x2) - atan(x1)) +
               ((1.0 - eta) * 0.5 * (gsl_sf_erf(x2) - gsl_sf_erf(x1))));
          norm += table[i + 1][k] * da;
        }
        for (k = 0; k < zbins; k++) {
          table[i + 1][k] =
              outfrac / da + (1. - outfrac) * table[i + 1][k] / norm;
        }
        break;
      case 3: // Gaussian photo-zs
        if (sigma_zphot_clustering(0., i) == 0.) {
          log_fatal("Lens galaxy photo-z model "
                 "underdetermined!\nredshift.clustering_photoz =3, but "
                 "nuisance.sigma_zphot_clustering[%d] not set", i);
          exit(1);
        }
        norm = 0.;
        for (k = 0; k < zbins; k++) {
          zi = z_v[k];
          x1 = (array[0] - zi + bias_zphot_clustering(zi, i)) /
               (sqrt(2.) * sigma_zphot_clustering(zi, i));
          x2 = (array[1] - zi + bias_zphot_clustering(zi, i)) /
               (sqrt(2.) * sigma_zphot_clustering(zi, i));
          table[i + 1][k] = 0.5 * pf_histo(zi, (void *)array) *
                            (gsl_sf_erf(x2) - gsl_sf_erf(x1));
          norm += table[i + 1][k] * da;
        }
        for (k = 0; k < zbins; k++) {
          table[i + 1][k] /= norm;
        }
        break;
      case 4: // histogram file contains n(z) estimates for each bin
        array[0] = 1.0 * i;
        pf_histo_n(0., (void *)array);
        norm = int_gsl_integrate_medium_precision(
            pf_histo_n, (void *)array, tomo.clustering_zmin[i],
            tomo.clustering_zmax[i], NULL, 1024);
        if (norm == 0) {
          log_fatal("pf_photoz:norm(nz=%d)=0", i);
          exit(1);
        }
        for (k = 0; k < zbins; k++) {
          table[i + 1][k] = pf_histo_n(z_v[k], (void *)array) / norm;
        }
        break;
      default:
        log_fatal("%d not supported in this cosmolike version", redshift.clustering_photoz);
        exit(1);
      }
      NORM[i] = norm;
    }
    // calculate normalized overall redshift distribution (without bins), store
    // in table[0][:]
    norm = 0;
    for (i = 0; i < tomo.clustering_Nbin; i++) {
      norm += NORM[i];
    }
    for (k = 0, zi = zhisto_min; k < zbins; k++, zi += da) {
      table[0][k] = 0;
      for (i = 0; i < tomo.clustering_Nbin; i++) {
        table[0][k] += table[i + 1][k] * NORM[i] / norm;
      }
    }
    for (i = -1; i < tomo.clustering_Nbin; i++) {
      gsl_spline_init(photoz_splines[i + 1], z_v, table[i + 1], zbins);
    }
  }
  if (j >= tomo.clustering_Nbin) {
    log_fatal("pf_photoz(z,%d) outside tomo.clustering_Nbin range", j);
    exit(1);
  }
  if (redshift.clustering_photoz == 4) {
    zz = zz - nuisance.bias_zphot_clustering[j];
  }
  if (zz <= z_v[0] || zz >= z_v[zbins - 1]) {
    return 0.0;
  }
  //VM COBAYA_SAMPLER - CHANGE ACCEL TO NULL (THREAD SAFE)
  return gsl_spline_eval(photoz_splines[j + 1], zz, NULL);
}

// ------------------------------------------------------------------------
// routines calculating the number of source and lens galaxies per bin
// ------------------------------------------------------------------------

double int_nsource(double z, void *param __attribute__((unused))) {
  return zdistr_photoz(z, -1);
}

// returns n_gal for shear tomography bin j, works only with binned
// distributions; j =-1 -> no tomography; j>= 0 -> tomography bin j
double nsource(int j) {
  static double **table = 0;
  double array[3];
  int i;
  if (table == 0 || table[0][0] != survey.n_gal) {
    array[0] = zdistr_photoz(0., 0);
    table = create_double_matrix(0, tomo.shear_Nbin, 0, 1);
    if (redshift.shear_photoz == 4) {
      printf("redshift.shear_photoz = 4, using tabulated tomo.n_source =");
      for (i = 0; i < tomo.shear_Nbin; i++) {
        printf(" %e,", tomo.n_source[i]);
        table[i + 1][0] = tomo.n_source[i];
        if (tomo.n_source[i] < 0.01 || tomo.n_source[i] > 100.) {
          printf("\n\n!!!!!!!!!\n\nredshift.shear_photoz = 4, "
                 "tomo.n_source[%d] = %e, EXIT\n\n!!!!!!!!!\n",
                 i, tomo.n_source[i]);
          exit(0);
        }
      }
      printf("\n");
    } else {
      double res, norm;
      norm = int_gsl_integrate_medium_precision(
          int_nsource, NULL, redshift.shear_zdistrpar_zmin,
          redshift.shear_zdistrpar_zmax, NULL, 1024);
      for (i = 0; i < tomo.shear_Nbin; i++) {
        res = int_gsl_integrate_medium_precision(
            int_nsource, (void *)array, tomo.shear_zmin[i], tomo.shear_zmax[i],
            NULL, 1024);
        table[i + 1][0] = res / norm * survey.n_gal;
        if (table[i + 1][0] < 0.01 || table[i + 1][0] > 100) {
          log_fatal("nsource(%d) = %e\n outside expected range", i, table[i + 1][0]);
          exit(1);
        }
      }
    }
    table[0][0] = survey.n_gal;
  }
  return table[j + 1][0];
}

double int_nlens(double z, void *param __attribute__((unused))) {
  return pf_photoz(z, -1);
}

// returns n_gal for clustering tomography bin j, works only
// with binned distributions; j =-1 -> no tomography; j>= 0
// -> tomography bin j
double nlens(int j) {
  static double **table = 0;
  double array[3];
  int i;
  if (table == 0 || table[0][0] != survey.n_lens) {
    array[0] = pf_photoz(0., 0);
    table = create_double_matrix(0, tomo.clustering_Nbin, 0, 1);
    if (redshift.clustering_photoz == 4) {
      printf("redshift.clustering_photoz = 4, using tabulated tomo.n_lens =");
      for (i = 0; i < tomo.clustering_Nbin; i++) {
        printf(" %e,", tomo.n_lens[i]);
        table[i + 1][0] = tomo.n_lens[i];
      }
      printf("\n");
    } else {
      double res, norm;
      array[0] = tomo.clustering_zmin[0];
      array[1] = tomo.clustering_zmax[tomo.clustering_Nbin - 1];
      norm = int_gsl_integrate_medium_precision(int_nlens, NULL, array[0],
                                                array[1], NULL, 1024);
      for (i = 0; i < tomo.clustering_Nbin; i++) {
        array[0] = tomo.clustering_zmin[i];
        array[1] = tomo.clustering_zmax[i];
        res = int_gsl_integrate_medium_precision(
            int_nlens, (void *)array, array[0], array[1], NULL, 1024);
        table[i + 1][0] = res / norm * survey.n_lens;
        if (table[i + 1][0] < 0.001 || table[i + 1][0] > 100) {
          log_fatal("nlens(%d) = %e\n outside expected range", i, table[i + 1][0]);
          exit(1);
        }
      }
    }
    table[0][0] = survey.n_lens;
  }
  return table[j + 1][0];
}

double int_for_zmean(double z, void *params) {
  double *array = (double *)params;
  return z * pf_photoz(z, (int)array[0]);
}

double norm_for_zmean(double z, void *params) {
  double *array = (double *)params;
  return pf_photoz(z, (int)array[0]);
}

// mean true redshift of galaxies in tomography bin j
double zmean(int j) {
  static double **table = 0;
  if (table == 0) {
    {
      double array[1];
      array[0] = pf_photoz(0., 0);
    }
    table = create_double_matrix(0, tomo.clustering_Nbin, 0, 1);
    #pragma omp parallel for
    for (int i = 0; i < tomo.clustering_Nbin; i++) {
      double array[1];
      array[0] = 1.0 * i;
      table[i][0] =
          int_gsl_integrate_low_precision(
            int_for_zmean,
            (void *) array,
            tomo.clustering_zmin[i],
            tomo.clustering_zmax[i],
            NULL,
            1024
          ) /
          int_gsl_integrate_low_precision(
            norm_for_zmean,
            (void *) array,
            tomo.clustering_zmin[i],
            tomo.clustering_zmax[i],
            NULL,
            1024
          );
    }
  }
  return table[j][0];
}

double int_for_zmean_source(double z, void *params) {

  double *array = (double *)params;
  return z * zdistr_photoz(z, (int)array[0]);
}

// mean true redshift of source galaxies in tomography bin j
double zmean_source(int j) {
  static double **table = 0;
  if (table == 0) {
    {
      double array[1];
      array[0] = zdistr_photoz(0., 0);
    }
    table = create_double_matrix(0, tomo.shear_Nbin, 0, 1);
    {
      const int i = 0;
      double array[1];
      array[0] = 1.0 * i;
      table[i][0] = int_gsl_integrate_medium_precision(
          int_for_zmean_source, (void *)array, tomo.shear_zmin[i],
          tomo.shear_zmax[i], NULL, 1024);
    }
    #pragma omp parallel for
    for (int i = 1; i < tomo.shear_Nbin; i++) {
      double array[1];
      array[0] = 1.0 * i;
      table[i][0] = int_gsl_integrate_medium_precision(
          int_for_zmean_source, (void *)array, tomo.shear_zmin[i],
          tomo.shear_zmax[i], NULL, 1024);
    }
  }
  return table[j][0];
}

double int_for_ggl_efficiency(double z, void *params) {
  double *array = (double *)params;
  int zs = (int)array[1];
  return pf_photoz(z, (int)array[0]) * g_tomo(1. / (1. + z), zs) * (1. + z) *
         f_K(chi(1. / (1. + z)));
}

double max_g_tomo(int zs) {
  double g, max = 0;
  for (double z = 0.; z < tomo.shear_zmax[zs]; z += tomo.shear_zmax[zs] / 50.) {
    g = g_tomo(1. / (1. + z), zs) * (1. + z) * f_K(chi(1. / (1. + z)));
    if (g > max) {
      max = g;
    }
  }
  return max;
}

double ggl_efficiency(int zl, int zs) {
  static double **table = 0;
  if (table == 0) {
    {
      double init = pf_photoz(0, 0);
      init = zdistr_photoz(0, 0);
    }
    table = create_double_matrix(0, tomo.clustering_Nbin, 0, tomo.shear_Nbin);
    {
      const int i = 0;
      double array[2];
      array[0] = 1.0 * i;
      for (int j = 0; j < tomo.shear_Nbin; j++) {
        array[1] = 1.0 * j;
        table[i][j] =
            int_gsl_integrate_medium_precision(
              int_for_ggl_efficiency,
              (void *) array,
              tomo.clustering_zmin[i],
              tomo.clustering_zmax[i],
              NULL,
              1024
            ) / max_g_tomo(j);
      }
    }
    #pragma omp parallel for
    for (int i = 1; i < tomo.clustering_Nbin; i++) {
      double array[2];
      array[0] = 1.0 * i;
      for (int j = 0; j < tomo.shear_Nbin; j++) {
        array[1] = 1.0 * j;
        table[i][j] =
            int_gsl_integrate_medium_precision(
              int_for_ggl_efficiency,
              (void *) array,
              tomo.clustering_zmin[i],
              tomo.clustering_zmax[i],
              NULL,
              1024
            ) /max_g_tomo(j);
      }
    }
  }
  return table[zl][zs];
}

// ------------------------------------------------------------------------
// integrands for bin-averaged lens efficiencies
// ------------------------------------------------------------------------

double int_for_g_tomo(double aprime, void *params) {
  double *ar = (double *)params;
  const int zbin = (int) ar[0];
  const double chi1 = ar[1];
  const double chi_prime = chi(aprime);

  return zdistr_photoz(1. / aprime - 1., zbin) * f_K(chi_prime - chi1) /
    f_K(chi_prime) / (aprime * aprime);
}

// for tomography bin zbin
double g_tomo(double a, int zbin) {
  static nuisancepara N;
  static cosmopara C;
  static double **table = 0;
  static double da = 0.0;
  if (table == 0 || recompute_zphot_shear(N) || recompute_cosmo3D(C)) {
    if (table == 0) {
      table = create_double_matrix(0, tomo.shear_Nbin, 0, Ntable.N_a - 1);
    }

    da = (0.999999 - 1. / (redshift.shear_zdistrpar_zmax + 1.)) /(Ntable.N_a - 1);

    { // j=-1, no tomography is being done
      const int j = -1;
      double array[2];
      array[0] = (double)j; // if j=-1, no tomography is being done
      double aa = 1. / (redshift.shear_zdistrpar_zmax + 1.);
      for (int i = 0; i < Ntable.N_a; i++, aa += da) {
        array[1] = chi(aa);
        table[j + 1][i] = int_gsl_integrate_medium_precision(
          int_for_g_tomo,
          (void *) array,
          1. / (redshift.shear_zdistrpar_zmax + 1.),
          aa,
          NULL,
          4000
        );
      }
    }
    if (tomo.shear_Nbin > 0) {
      const int j = 0;
      double array[2];
      array[0] = (double) j;
      double aa = 1. / (redshift.shear_zdistrpar_zmax + 1.);
      for (int i = 0; i < Ntable.N_a; i++, aa += da) {
        array[1] = chi(aa);
        table[j + 1][i] = int_gsl_integrate_medium_precision(
          int_for_g_tomo,
          (void *) array,
          1. / (redshift.shear_zdistrpar_zmax + 1.),
          aa,
          NULL,
          4000
        );
      }
      #pragma omp parallel for
      for (int j = 1; j < tomo.shear_Nbin; j++) {
        double array[2];
        array[0] = (double) j;
        double aa = 1. / (redshift.shear_zdistrpar_zmax + 1.);
        for (int i = 0; i < Ntable.N_a; i++, aa += da) {
          array[1] = chi(aa);
          table[j + 1][i] = int_gsl_integrate_medium_precision(
            int_for_g_tomo,
            (void *) array,
            1. / (redshift.shear_zdistrpar_zmax + 1.),
            aa,
            NULL,
            4000
          );
        }
      }
	}
    update_nuisance(&N);
    update_cosmopara(&C);
  }
  if (a <= 1. / (redshift.shear_zdistrpar_zmax + 1.) || a > 1.0 - da) {
    return 0.0;
  }
  return interpol(
    table[zbin + 1],
    Ntable.N_a,
	1. / (redshift.shear_zdistrpar_zmax + 1.),
	0.999999,
	da,
	a,
	1.0,
	1.0
  ); // zbin =-1 is non-tomography
}

// \int n(z') W(z,z')^2 routines for source clustering
double int_for_g2_tomo(double aprime, void *params) {
  double *ar = (double *) params;
  const int zbin = (int) ar[0];
  const double chi1 = ar[1];
  const double chi_prime = chi(aprime);
  return zdistr_photoz(1./aprime-1.,zbin) * f_K(chi_prime - chi1)/
    f_K(chi_prime - chi1)/(f_K(chi_prime)*f_K(chi_prime))/(aprime*aprime);
}

// for tomography bin zbin
double g2_tomo(double a, int zbin) {
  static nuisancepara N;
  static cosmopara C;
  static double **table = 0;
  static double da = 0.0;
  if (table == 0 || recompute_zphot_shear(N) || recompute_cosmo3D(C)) {
    if (table == 0) {
      table = create_double_matrix(0, tomo.shear_Nbin, 0, Ntable.N_a-1);
    }

    da = (0.999999-1./(redshift.shear_zdistrpar_zmax+1.))/(Ntable.N_a-1);

	{ // j=-1, no tomography is being done
	  const int j = -1;
      double array[2];
      array[0] = (double) j;
      double aa = 1./(redshift.shear_zdistrpar_zmax + 1.);
      for (int i = 0; i < Ntable.N_a; i++, aa+=da) {
        array[1] = chi(aa);
        table[j+1][i] = int_gsl_integrate_medium_precision(
           int_for_g2_tomo,
           (void*) array,1./(redshift.shear_zdistrpar_zmax+1.),
           aa,
           NULL,
           4000
        );
      }
	}
	if (tomo.shear_Nbin > 0) {
      const int j = 0;
      double array[2];
      array[0] = (double) j;
      double aa = 1./(redshift.shear_zdistrpar_zmax + 1.);
      for (int i = 0; i < Ntable.N_a; i++, aa+=da) {
        array[1] = chi(aa);
        table[j+1][i] = int_gsl_integrate_medium_precision(
           int_for_g2_tomo,
           (void*) array,1./(redshift.shear_zdistrpar_zmax+1.),
           aa,
           NULL,
           4000
        );
      }
      #pragma omp parallel for
      for (int j = 1; j < tomo.shear_Nbin; j++) {
	    double array[2];
        array[0] = (double) j;
        double aa = 1./(redshift.shear_zdistrpar_zmax+1.);
        for (int i = 0; i < Ntable.N_a; i++,aa+=da) {
          array[1] = chi(aa);
          table[j+1][i] = int_gsl_integrate_medium_precision(
             int_for_g2_tomo,
             (void*) array,1./(redshift.shear_zdistrpar_zmax+1.),
             aa,
             NULL,
             4000
          );
        }
      }
    }
    update_nuisance(&N);
    update_cosmopara(&C);
  }
  if (a <= 1./(redshift.shear_zdistrpar_zmax+1.) || a > 1.0-da) {
    return 0.0;
  }
  return interpol(
  	table[zbin+1],
  	Ntable.N_a,
  	1./(redshift.shear_zdistrpar_zmax+1.),
  	0.999999,
  	da,
  	a,
  	1.0,
  	1.0
  ); //zbin =-1 is non-tomography
}

double int_for_g_lens(double aprime, void *params) {
  double *ar = (double *)params;
  const int zbin = (int) ar[0];
  const double chi1 = ar[1];
  const double chi_prime = chi(aprime);
  return pf_photoz(1. / aprime - 1., zbin) * f_K(chi_prime - chi1) /
        f_K(chi_prime) / (aprime * aprime);
}

// for *lens* tomography bin zbin
double g_lens(double a, int zbin) {
  static nuisancepara N;
  static cosmopara C;
  static double **table = 0;
  static double da = 0.0;
  if (table == 0 || recompute_zphot_clustering(N) || recompute_cosmo3D(C)) {
    if (table == 0) {
      table = create_double_matrix(0, tomo.clustering_Nbin, 0, Ntable.N_a - 1);
    }

    da = (0.999999 - 1. / (redshift.clustering_zdistrpar_zmax + 1.)) /
         (Ntable.N_a - 1);
    {
      const int j = -1;
      double array[2];
      array[0] = (double) j; // j=-1, no tomography is being done
      double aa = 1./(redshift.clustering_zdistrpar_zmax + 1.);
      for (int i = 0; i < Ntable.N_a; i++, aa += da) {
        array[1] = chi(aa);
        table[j + 1][i] = int_gsl_integrate_medium_precision(
          int_for_g_lens,
          (void *) array,
          1. / (redshift.shear_zdistrpar_zmax + 1.),
          aa,
          NULL,
          4000
        );
      }
    }
    if (tomo.clustering_Nbin > 0) {
      {
        const int j = 0;
        double array[2];
        array[0] = (double) j;
        double aa = 1./(redshift.clustering_zdistrpar_zmax + 1.);
        for (int i = 0; i < Ntable.N_a; i++, aa += da) {
          array[1] = chi(aa);
          table[j + 1][i] = int_gsl_integrate_medium_precision(
            int_for_g_lens,
            (void *) array,
            1. / (redshift.shear_zdistrpar_zmax + 1.),
            aa,
            NULL,
            4000
          );
        }
      }
      #pragma omp parallel for
      for (int j = 1; j < tomo.clustering_Nbin; j++) {
        double array[2];
        array[0] = (double) j; // if j=-1, no tomography is being done
        double aa = 1./(redshift.clustering_zdistrpar_zmax + 1.);
        for (int i = 0; i < Ntable.N_a; i++, aa += da) {
          array[1] = chi(aa);
          table[j + 1][i] = int_gsl_integrate_medium_precision(
            int_for_g_lens,
            (void *)array,
            1. / (redshift.shear_zdistrpar_zmax + 1.),
            aa,
            NULL,
            4000
          );
        }
      }
    }
    update_nuisance(&N);
    update_cosmopara(&C);
  }
  if (a <= 1. / (redshift.clustering_zdistrpar_zmax + 1.) || a > 1.0 - da) {
    return 0.0;
  }
  return interpol(
    table[zbin + 1],
    Ntable.N_a,
    1. / (redshift.clustering_zdistrpar_zmax + 1.),
    0.999999,
    da,
    a,
    1.0,
    1.0
  ); // zbin =-1 is non-tomography
}

double g_cmb(double a) {
  static cosmopara C;
  static double chi_cmb = 0.;
  static double fchi_cmb = 0.;
  if (recompute_cosmo3D(C)) {
    chi_cmb = chi(1./1091.);
    fchi_cmb = f_K(chi_cmb);
    update_cosmopara(&C);
  }
  return f_K(chi_cmb-chi(a))/fchi_cmb;
}

double g_bg(double a __attribute__((unused)), int nzlens __attribute__((unused))) {
  printf("g_bg option no longer supported\n'");
  return 0.;
}
