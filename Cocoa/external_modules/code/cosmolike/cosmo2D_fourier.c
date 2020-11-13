#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../log.c/src/log.h"

#include "basics.h"
#include "cosmo3D.h"
#include "cosmo2D_fourier.h"
#include "bias.h"
#include "pt_cfastpt.h"
#include "radial_weights.h"
#include "recompute.h"
#include "redshift_spline.h"
#include "structs.h"

// ----------------------------------------------------------------------------
// Limber approximation integrands for angular power spectra
// use W_i W_j/fK^2 dchi_da(a) as Limber weight for power spectra
// ----------------------------------------------------------------------------

double int_for_C_cl_tomo(double a, void *params) {
#ifdef DEBUG
  if (a >= 1.0) {
    log_fatal("a>=1 in int_for_C_cl_tomo");
    exit(1);
  }
#endif

  double *ar = (double *)params;
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double ell = ar[2] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell / fK;

  const double res = W_gal(a, ar[0], chidchi.chi, hoverh0) *
    W_gal(a, ar[1], chidchi.chi, hoverh0) * chidchi.dchida/(fK * fK);
  return res * Pdelta(k, a);
}

double int_for_C_cl_tomo_b2(double a, void *params) {
#ifdef DEBUG
  if (a >= 1.0) {
    log_fatal("a>=1 in int_for_C_cl_tomo_b2");
    exit(1);
  }
#endif

  double *ar = (double *) params;
  const double b1 = gbias.b1_function(1. / a - 1., (int)ar[0]);
  const double b2 = gbias.b2[(int)ar[0]];
  const double bs2 = gbias.bs2[(int)ar[0]];
  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double g4 = growfac_a*growfac_a*growfac_a*growfac_a;
  const double ell = ar[2] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell / fK;
  const double PK = Pdelta(k, a);
  double s4 = 0.; // PT_sigma4(k);

  double res =
    W_HOD(a, ar[0], hoverh0)*W_HOD(a, ar[1], hoverh0)*chidchi.dchida/(fK*fK);

  if (res) {
    res = res * (b1 * b1 * PK +
                 g4 * (b1 * b2 * PT_d1d2(k) +
                       0.25 * b2 * b2 * (PT_d2d2(k) - 2. * s4) +
                       b1 * bs2 * PT_d1s2(k) +
                       0.5 * b2 * bs2 * (PT_d2s2(k) - 4. / 3. * s4) +
                       .25 * bs2 * bs2 * (PT_s2s2(k) - 8. / 9. * s4) +
                       b1 * b3nl_from_b1(b1) * PT_d1d3(k)));
  }
  res += (W_gal(a, ar[0], chidchi.chi, hoverh0) * W_mag(a, fK, ar[1]) +
          W_gal(a, ar[1], chidchi.chi, hoverh0) * W_mag(a, fK, ar[0])) *
         chidchi.dchida / fK / fK * PK;
  return res;
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// angular power spectra - without look-up tables
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// galaxy clustering power spectrum of galaxy bins ni, nj
double C_cl_tomo_nointerp(double l, int ni, int nj) {
  static int init = -1;
  double array[3] = {1.0 * ni, 1.0 * nj, l};
  if (gbias.b2[ni] || gbias.b2[nj]) {
    if (ni != nj) {
      if (init == -1) {
        log_info(
          "Called C_cl(l,z1=%d,z2=%d) with non-linear bias parameters set.",
          ni,
          nj
        );
        log_info(
          "Cross-clustering beyond linear bias for cross-tomography bins "
          "not yet supported."
        );
        log_info("Use linear bias only for z1 != z2 clustering.");
        init = 1;
      }
      // COCOA: WE CHANGE INTEGRATION FROM MEDIUM TO LOW PRECISION!!
      return int_gsl_integrate_low_precision(
        int_for_C_cl_tomo,
        (void *) array,
        fmax(amin_lens(ni),
        amin_lens(nj)),
        fmin(amax_lens(ni),
        amax_lens(nj)),
        NULL,
        1000
      );
    }
    // COCOA: WE CHANGE INTEGRATION FROM MEDIUM TO LOW PRECISION!!
    return int_gsl_integrate_low_precision(
      int_for_C_cl_tomo_b2,
      (void *) array,
      amin_lens(ni),
      amax_lens(ni),
      NULL,
      1000
    );
  } else {
    if (ni == nj) {
      return int_gsl_integrate_medium_precision(
        int_for_C_cl_tomo,
        (void *) array,
        amin_lens(ni),
        0.999999,
        NULL,
        1000
      );
    }
    // COCOA: WE CHANGE INTEGRATION FROM MEDIUM TO LOW PRECISION!!
    return int_gsl_integrate_low_precision(
      int_for_C_cl_tomo,
      (void *)array,
      amin_lens(nj),
      0.99999,
      NULL,
      1000
    ); // zi<=zj
  }
}

// ---------------------------------------------------------------------------
// angular power spectra - with look-up tables
// ---------------------------------------------------------------------------

// galaxy clustering power spectrum of galaxies in bins ni, nj
double C_cl_tomo(double l, int ni, int nj) {
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
#ifdef DEBUG
  if (ni < 0 ||
    ni >= tomo.clustering_Nbin ||
    nj < 0 ||
    nj >= tomo.clustering_Nbin
  ) {
    log_fatal("C_cl_tomo(l,%d,%d) outside tomo.clustering_Nbin range", ni, nj);
    exit(1);
  }
#endif

  if (recompute_clustering(C, G, N, ni, nj)) {
    if (table == 0) {
      table = create_double_matrix(
          0, tomo.clustering_Nbin * tomo.clustering_Nbin - 1, 0,
          Ntable.N_ell - 1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin) / (Ntable.N_ell);
    }
    for (int i = 0; i < tomo.clustering_Nbin * tomo.clustering_Nbin; i++) {
      table[i][0] = 123456789.0;
    }
    update_cosmopara(&C);
    update_nuisance(&N);
    update_galpara(&G);
  }

  int j = ni * tomo.clustering_Nbin + nj;

  if (table[j][0] > 123456780.0) {
    // still need to recompute this tomography bin combination
    {
      const int i = 0;
      const double llog = logsmin + i*ds;
      const double result = C_cl_tomo_nointerp(exp(llog), ni, nj);
      if (result <= 0) {
        table[j][i] = -100;
      }
      else {
        table[j][i] = log(result);
      }
      table[nj * tomo.clustering_Nbin + ni][i] = table[j][i];
    }
    #pragma omp parallel for
    for (int i = 1; i < Ntable.N_ell; i++) {
      const double llog = logsmin + i*ds;
      const double result = C_cl_tomo_nointerp(exp(llog), ni, nj);
      if (result <= 0) {
        table[j][i] = -100;
      }
      else {
        table[j][i] = log(result);
      }
      table[nj * tomo.clustering_Nbin + ni][i] = table[j][i];
    }
  }

  double f1 = exp(interpol_fitslope(table[j], Ntable.N_ell, logsmin, logsmax,
                                    ds, log(l), 1.));
  if (isnan(f1)) {
    f1 = 0.;
  }

  return f1;
}
