#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bias.h"
#include "basics.h"
#include "cosmo2D_fullsky_TATT.h"
#include "cosmo3D.h"
#include "pt_cfastpt.h"
#include "radial_weights.h"
#include "recompute.h"
#include "redshift_spline.h"
#include "structs.h"

#include "log.c/src/log.h"

// ell_max for transform to angular correlation functions
static int LMAX = 100000;
// ell_min for switching from exact evalution of C(ell) to interpolated look-up
// table
static int LMIN_tab = 20;
// number of grid point for C(ell) look-up tables
static int NTAB_TATT = 60;

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Limber Integrands for Shear and GGL
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// NLA/TA amplitude C1, nz argument only need if per-bin amplitude
double C1_TA(double a, double nz, double growfac_a) {
  // per-bin IA parameters
  if (like.IA == 3 || like.IA == 5) {
    return -nuisance.A_z[(int)nz]*cosmology.Omega_m*nuisance.c1rhocrit_ia/ growfac_a;
  }
  // power law evolution
  return -cosmology.Omega_m * nuisance.c1rhocrit_ia /
         growfac_a * nuisance.A_ia *
         pow(1. / (a * nuisance.oneplusz0_ia), nuisance.eta_ia);
}

// TA source bias parameter, nz argument only need if per-bin amplitude
double b_TA(double a __attribute__((unused)), double nz) {
  // per-bin IA parameters
  if (like.IA == 5) {
    return nuisance.b_ta_z[(int)nz];
  }
  // power law evolution
  return nuisance.b_ta_z[0];
}

// TT amplitude C2, nz argument only need if per-bin amplitude
double C2_TT(double a, double nz, double growfac_a) {
  // per-bin IA parameters
  if (like.IA == 5) {
    return 5. * nuisance.A2_z[(int)nz] * cosmology.Omega_m *
           nuisance.c1rhocrit_ia * pow(1.0/growfac_a, 2.0);
  }
  // power law evolution
  return 5. * nuisance.A2_ia * cosmology.Omega_m * nuisance.c1rhocrit_ia *
         (1.0 /(growfac_a*growfac_a)) *
         pow(1. / (a * nuisance.oneplusz0_ia), nuisance.eta_ia_tt);
}

double int_for_C_shear_shear_IA_EE(double a, void *params) {
  if (a >= 1.0) {
    log_fatal("a >= 1");
    exit(1);
  }

  double *ar = (double *) params;

  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  double hoverh0 = hoverh0v2(a, chidchi.dchida);

  const double ell = ar[2] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell / fK;

  // radial n_z weight for first source bin (for use with IA term)
  const double ws1 = W_source(a, ar[0], hoverh0);

  // radial n_z weight for second source bin (for use with IA term)
  const double ws2 = W_source(a, ar[1], hoverh0);
  const double wk1 = W_kappa(a, fK, ar[0]); // radial lens efficiency for first source bin
  const double wk2 = W_kappa(a, fK, ar[1]); // radial lens efficiency for second source bin

  // IA parameters for first source bin
  const double C1 = C1_TA(a, ar[0], growfac_a);
  const double b_ta = b_TA(a, ar[0]);
  const double C2 = C2_TT(a, ar[0], growfac_a);

  // IA parameters for second source bin
  const double C1_2 = C1_TA(a, ar[1], growfac_a);
  const double b_ta_2 = b_TA(a, ar[1]);
  const double C2_2 = C2_TT(a, ar[1], growfac_a);

  // GG cosmic shear
  const double pdelta_ak = Pdelta(k, a);
  double res = wk1 * wk2 * pdelta_ak;
  //COCOA: Took these evaluations of the parenthesis - to force them to update
  //COCOA: the static variables in the first call that is done outside OpenMP loop
  const double tmp1 = TATT_II_EE(k, a, C1, C2, b_ta, C1_2, C2_2, b_ta_2, growfac_a, pdelta_ak);
  const double tmp2 = TATT_GI_E(k, a, C1, C2, b_ta, growfac_a, pdelta_ak);
  const double tmp3 = TATT_GI_E(k, a, C1_2, C2_2, b_ta_2, growfac_a, pdelta_ak);
  if (C1 || C1_2 || C2 || C2_2) {
    // II contribution
    res += ws1 * ws2 * tmp1;
    // GI contribution
    res += ws1 * wk2 * tmp2 + ws2 * wk1 * tmp3;
  }
  return res * chidchi.dchida / (fK * fK);
}

double int_for_C_shear_shear_IA_BB(double a, void *params) {
  if (a >= 1.0) {
    log_fatal("a >= 1");
    exit(1);
  }

  double *ar = (double *) params;

  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);

  const double ell = ar[2] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell / fK;

  // radial n_z weight for first source bin (for use with IA term)
  const double ws1 = W_source(a, ar[0], hoverh0);

  // radial n_z weight for second source bin (for use with IA term)
  const double ws2 = W_source(a, ar[1], hoverh0);

  // IA parameters for first source bin
  const double C1 = C1_TA(a, ar[0], growfac_a);
  const double b_ta = b_TA(a, ar[0]);
  const double C2 = C2_TT(a, ar[0], growfac_a);

  // IA parameters for second source bin
  const double C1_2 = C1_TA(a, ar[1], growfac_a);
  const double b_ta_2 = b_TA(a, ar[1]);
  const double C2_2 = C2_TT(a, ar[1], growfac_a);

  //COCOA: Took these evaluations of the parenthesis - to force them to update
  //COCOA: the static variables in the first call that is done outside OpenMP loop
  double res = 0.;
  const double tmp1 = TATT_II_BB(k, a, C1, C2, b_ta, C1_2, C2_2, b_ta_2, growfac_a);
  if ((b_ta || C2) && (b_ta_2 || C2_2)) {
    res = ws1 * ws2 * tmp1;
  }
  return res * chidchi.dchida / (fK * fK);
}

double int_for_C_ggl_IA_TATT(double a, void *params) {
  if (a >= 1.0) {
    log_fatal("a >= 1");
    exit(1);
  }

  double *ar = (double *) params;

  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);

  const double ell = ar[2] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell / fK;

  // radial n_z weight for source bin (for use with IA term)
  const double ws = W_source(a, ar[1], hoverh0);
  const double wk = W_kappa(a, fK, ar[1]); // radial lens efficiency for source bin

  // IA parameters for first source bin
  const double C1 = C1_TA(a, ar[1], growfac_a);
  const double b_ta = b_TA(a, ar[1]);
  const double C2 = C2_TT(a, ar[1], growfac_a);

  // radial n_z weight for lens bin (for use with clustering term)
  const double w_density = W_HOD(a, ar[0], hoverh0);

  // lens efficiency *b_mag for lens bin (for lens magnification)
  const double w_mag = W_mag(a, fK, ar[0]) * gbias.b_mag[(int)ar[0]];

  // galaxy bias parameters for lens bin
  const double b1 = gbias.b1_function(1. / a - 1., (int) ar[0]);
  const double b2 = gbias.b2[(int)ar[0]];
  const double bs2 = gbias.bs2[(int)ar[0]];

  const double g4 = growfac_a*growfac_a*growfac_a*growfac_a;
  const double Pnl = Pdelta(k, a);
  double P_1loop = b1 * Pnl;
  if (w_density * b2 != 0) {
    P_1loop += g4 * (0.5 * b2 * PT_d1d2(k) + 0.5 * bs2 * PT_d1s2(k) +
                     0.5 * b3nl_from_b1(b1) * PT_d1d3(k));
  }

  // 1-loop P_gm ggl terms
  double res = w_density * wk * P_1loop;
  // lens magnification x G term
  res += w_mag * wk * Pnl;
  // (linear bias lens density + lens magnification) with TATT_GI terms

  //COCOA: Took these evaluations of the parenthesis - to force them to update
  //COCOA: the static variables in the first call that is done outside OpenMP loop
  const double tmp1 = TATT_GI_E(k, a, C1, C2, b_ta, growfac_a, Pnl);
  if (C1 || C2) {
    res += (b1 * w_density + w_mag) * ws * tmp1;
  }
  return res * chidchi.dchida / fK / fK;
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Angular Power Spectra
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// Cocoa: CHANGED THE INTEGRATION FROM MEDIUM TO LOW PRECISION!!
double C_EE_TATT(double l, int ni, int nj) {
  double array[3] = {(double)ni, (double)nj, l};
  return int_gsl_integrate_low_precision(
    int_for_C_shear_shear_IA_EE,
    (void *) array,
    fmax(amin_source(ni), amin_source(nj)),
    amax_source(ni),
    NULL,
    1000
  );
}

// Cocoa: CHANGED THE INTEGRATION FROM MEDIUM TO LOW PRECISION!!
double C_BB_TATT(double l, int ni, int nj) {
  double array[3] = {(double)ni, (double)nj, l};
  return int_gsl_integrate_low_precision(
    int_for_C_shear_shear_IA_BB,
    (void *) array,
    fmax(amin_source(ni), amin_source(nj)),
    fmin(amax_source_IA(ni), amax_source_IA(nj)),
    NULL,
    1000
  );
}

// Cocoa: CHANGED THE INTEGRATION FROM MEDIUM TO LOW PRECISION!!
double C_ggl_TATT(double l, int nl, int ns) {
  double array[3] = {(double)nl, (double)ns, l};
  return int_gsl_integrate_low_precision(
    int_for_C_ggl_IA_TATT,
    (void *) array,
    amin_lens(nl),
    amax_lens(nl),
    NULL,
    1000
  );
}


// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Tabulated Angular Power Spectra
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// shear power spectrum of source galaxies in bins ni, nj
double C_EE_tab(double l, int ni, int nj) {
  static cosmopara C;
  static nuisancepara N;
  static double **table, *sig;
  static int osc[100];
  static double ds = .0, logsmin = .0, logsmax = .0;
  if (ni < 0 || ni >= tomo.shear_Nbin || nj < 0 || nj >= tomo.shear_Nbin) {
    log_fatal("C_shear_shear_EE(l,%d,%d) outside tomo.shear_Nbin range", ni, nj);
    exit(1);
  }
  if (recompute_shear(C, N)) {
    if (table == 0) {
      table = create_double_matrix(0, tomo.shear_Npowerspectra - 1, 0,
                                   NTAB_TATT - 1);
      sig = create_double_vector(0, tomo.ggl_Npowerspectra - 1);
      logsmin = log(fmax(LMIN_tab - 1., 1.0));
      logsmax = log(LMAX + 1);
      ds = (logsmax - logsmin) / (NTAB_TATT - 1.);
    }

    {
      const int k = 0;
      sig[k] = 1.;
      osc[k] = 0;
      if (C_EE_TATT(500., Z1(k), Z2(k)) < 0) {
        sig[k] = -1.;
      }
      #pragma omp parallel for
      for (int i = 0; i < NTAB_TATT; i++) {
        const double llog = logsmin + i*ds;
        table[k][i] = C_EE_TATT(exp(llog), Z1(k), Z2(k));
        if (table[k][i] * sig[k] < 0.) {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0) {
        for (int i = 0; i < NTAB_TATT; i++) {
          table[k][i] = log(sig[k] * table[k][i]);
        }
      }
    }
    #pragma omp parallel for
    for (int k = 1; k < tomo.shear_Npowerspectra; k++) {
      double llog = logsmin;
      sig[k] = 1.;
      osc[k] = 0;
      if (C_EE_TATT(500., Z1(k), Z2(k)) < 0) {
        sig[k] = -1.;
      }
      for (int i = 0; i < NTAB_TATT; i++, llog += ds) {
        table[k][i] = C_EE_TATT(exp(llog), Z1(k), Z2(k));
        if (table[k][i] * sig[k] < 0.) {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0) {
        for (int i = 0; i < NTAB_TATT; i++) {
          table[k][i] = log(sig[k] * table[k][i]);
        }
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  if (log(l) < logsmin || log(l) > logsmax) {
    log_fatal("C_EE_tab: l = %e outside look-up table range [%e,%e]", l, exp(logsmin), exp(logsmax));
    exit(1);
  }
  int k = N_shear(ni, nj);
  double f1;
  if (osc[k] == 0) {
    f1 = sig[k] * exp(interpol(table[k], NTAB_TATT, logsmin, logsmax, ds,
                               log(l), 0., 0.));
  } else {
    f1 = interpol(table[k], NTAB_TATT, logsmin, logsmax, ds, log(l), 0., 0.);
  }
  if (isnan(f1)) {
    f1 = 0.;
  }
  return f1;
}

// shear power spectrum of source galaxies in bins ni, nj
double C_BB_tab(double l, int ni, int nj) {
  static cosmopara C;
  static nuisancepara N;
  static double **table, *sig;
  static int osc[100];
  static double ds = .0, logsmin = .0, logsmax = .0;
  if (ni < 0 || ni >= tomo.shear_Nbin || nj < 0 || nj >= tomo.shear_Nbin) {
    log_fatal("C_shear_shear_BB(l,%d,%d) outside tomo.shear_Nbin range",
           ni, nj);
    exit(1);
  }
  if (recompute_shear(C, N)) {
    if (table == 0) {
      table = create_double_matrix(
        0,
        tomo.shear_Npowerspectra - 1,
        0,
        NTAB_TATT - 1
      );
      sig = create_double_vector(0, tomo.ggl_Npowerspectra - 1);
      logsmin = log(fmax(LMIN_tab - 1., 1.0));
      logsmax = log(LMAX + 1);
      ds = (logsmax - logsmin) / (NTAB_TATT - 1.);
    }

    {
      const int k = 0;
      sig[k] = 1.;
      osc[k] = 0;
      if (C_BB_TATT(500., Z1(k), Z2(k)) < 0) {
        sig[k] = -1.;
      }
      #pragma omp parallel for
      for (int i = 0; i < NTAB_TATT; i++) {
        const double llog = logsmin + i*ds;
        table[k][i] = C_BB_TATT(exp(llog), Z1(k), Z2(k));
        if (table[k][i] * sig[k] < 0.) {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0) {
        for (int i = 0; i < NTAB_TATT; i++) {
          table[k][i] = log(sig[k] * table[k][i]);
        }
      }
    }
    #pragma omp parallel for
    for (int k = 1; k < tomo.shear_Npowerspectra; k++) {
      sig[k] = 1.;
      osc[k] = 0;
      if (C_BB_TATT(500., Z1(k), Z2(k)) < 0) {
        sig[k] = -1.;
      }
      for (int i = 0; i < NTAB_TATT; i++) {
        const double llog = logsmin + i*ds;
        table[k][i] = C_BB_TATT(exp(llog), Z1(k), Z2(k));
        if (table[k][i] * sig[k] < 0.) {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0) {
        for (int i = 0; i < NTAB_TATT; i++) {
          table[k][i] = log(sig[k] * table[k][i]);
        }
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  const double lgl = log(l);
  if (log(l) < logsmin || log(l) > logsmax) {
    log_fatal("C_BB_tab: l = %e outside look-up table range [%e,%e]", l,
           exp(logsmin), exp(logsmax));
    exit(1);
  }
  int k = N_shear(ni, nj);
  double f1;
  if (osc[k] == 0) {
    f1 = sig[k] * exp(
      interpol(
        table[k],
        NTAB_TATT,
        logsmin,
        logsmax,
        ds,
        lgl,
        0.,
        0.
      )
    );
  } else {
    f1 = interpol(table[k], NTAB_TATT, logsmin, logsmax, ds, lgl, 0., 0.);
  }
  if (isnan(f1)) {
    f1 = 0.;
  }
  return f1;
}

// G-G lensing power spectrum, lens bin ni, source bin nj
double C_ggl_TATT_tab(double l, int ni, int nj) {
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double **table, *sig;
  static int osc[100];
  static double ds = .0, logsmin = .0, logsmax = .0;
  if (ni < 0 || ni >= tomo.clustering_Nbin || nj < 0 || nj >= tomo.shear_Nbin) {
    log_fatal("C_ggl_TATT_tab(l,%d,%d) outside tomo.X_Nbin range", ni, nj);
    exit(1);
  }
  if (recompute_ggl(C, G, N, ni)) {
    if (table == 0) {
      table = create_double_matrix(0, tomo.ggl_Npowerspectra-1, 0,
        NTAB_TATT-1);
      sig = create_double_vector(0, tomo.ggl_Npowerspectra-1);
      logsmin = log(fmax(LMIN_tab - 1., 1.0));
      logsmax = log(LMAX + 1);
      ds = (logsmax - logsmin) / (NTAB_TATT - 1.);
    }
    {
      const int k = 0;
      sig[k] = 1.;
      osc[k] = 0;
      if (C_ggl_TATT(500., ZL(k), ZS(k)) < 0) {
        sig[k] = -1.;
      }
      #pragma omp parallel for
      for (int i = 0; i < NTAB_TATT; i++) {
        const double llog = logsmin + i*ds;
        table[k][i] = C_ggl_TATT(exp(llog), ZL(k), ZS(k));
        if (table[k][i] * sig[k] < 0.) {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0) {
        for (int i = 0; i < NTAB_TATT; i++) {
          table[k][i] = log(sig[k] * table[k][i]);
        }
      }
    }
    #pragma omp parallel for
    for (int k = 1; k < tomo.ggl_Npowerspectra; k++) {
      sig[k] = 1.;
      osc[k] = 0;
      if (C_ggl_TATT(500., ZL(k), ZS(k)) < 0) {
        sig[k] = -1.;
      }
      for (int i = 0; i < NTAB_TATT; i++) {
        const double llog = logsmin + i*ds;
        table[k][i] = C_ggl_TATT(exp(llog), ZL(k), ZS(k));
        if (table[k][i] * sig[k] < 0.) {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0) {
        for (int i = 0; i < NTAB_TATT; i++) {
          table[k][i] = log(sig[k] * table[k][i]);
        }
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
    update_galpara(&G);
  }
  const double lgl = log(l);
  if (lgl < logsmin || lgl > logsmax) {
    log_fatal(
      "C_ggl_TATT_tab: l = %e outside look-up table range [%e,%e]",
      l,
      exp(logsmin),
      exp(logsmax)
    );
    exit(1);
  }
  int k = N_ggl(ni, nj);
  double f1 = 0.;
  if (test_zoverlap(ni, nj) && osc[k] == 0) {
    f1 = sig[k] * exp(
      interpol(
        table[k],
        NTAB_TATT,
        logsmin,
        logsmax,
        ds,
        lgl,
        0.,
        0.
      )
    );
  }
  if (test_zoverlap(ni, nj) && osc[k] == 1) {
    f1 = interpol(table[k], NTAB_TATT, logsmin, logsmax, ds, lgl, 0., 0.);
  }
  if (isnan(f1)) {
    f1 = 0;
  }
  return f1;
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Correlation Functions (real space)
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// shear tomography correlation functions
double xi_pm_TATT(int pm, int nt, int ni, int nj) {
  static double **Glplus = 0;
  static double **Glminus = 0;
  static double *xi_vec_plus = 0;
  static double *xi_vec_minus = 0;
  static cosmopara C;
  static nuisancepara N;
  if (like.Ntheta == 0) {
    log_fatal("like.Ntheta not initialized");
    exit(1);
  }
  if (Glplus == 0) {
    Glplus = create_double_matrix(0, like.Ntheta - 1, 0, LMAX - 1);
    Glminus = create_double_matrix(0, like.Ntheta - 1, 0, LMAX - 1);
    xi_vec_plus = create_double_vector(0,
      tomo.shear_Npowerspectra * like.Ntheta - 1);
    xi_vec_minus = create_double_vector(0,
      tomo.shear_Npowerspectra * like.Ntheta - 1);

    double xmin[like.Ntheta];
    double xmax[like.Ntheta];
    double logdt = (log(like.vtmax) - log(like.vtmin)) / like.Ntheta;
    for (int i = 0; i < like.Ntheta; i++) {
      xmin[i] = cos(exp(log(like.vtmin) + (i + 0.0) * logdt));
      xmax[i] = cos(exp(log(like.vtmin) + (i + 1.0) * logdt));
    }

    #pragma omp parallel for
    for (int i = 0; i < like.Ntheta; i++) {
      double *Pmin = create_double_vector(0, LMAX + 1);
      double *Pmax = create_double_vector(0, LMAX + 1);
      double *dPmin = create_double_vector(0, LMAX + 1);
      double *dPmax = create_double_vector(0, LMAX + 1);

      gsl_sf_legendre_Pl_deriv_array(LMAX, xmin[i], Pmin, dPmin);
      gsl_sf_legendre_Pl_deriv_array(LMAX, xmax[i], Pmax, dPmax);

      for (int l = 3; l < LMAX; l++) {
        Glplus[i][l] =
            (2. * l + 1) / (2. * M_PI * l * l * (l + 1) * (l + 1)) *
            (-l * (l - 1.) / 2 * (l + 2. / (2 * l + 1)) *
                    (Pmin[l - 1] - Pmax[l - 1]) -
                l * (l - 1.) * (2. - l) / 2 *
                    (xmin[i] * Pmin[l] - xmax[i] * Pmax[l]) +
                l * (l - 1) / (2 * l + 1) * (Pmin[l + 1] - Pmax[l + 1])
                + (4 - l) * (dPmin[l] - dPmax[l]) +
                (l + 2) * (xmin[i] * dPmin[l - 1] - xmax[i] * dPmax[l - 1] -
                           Pmin[l - 1] + Pmax[l - 1])
                + 2 * (l - 1) *
                      (xmin[i] * dPmin[l] - xmax[i] * dPmax[l] - Pmin[l] +
                       Pmax[l]) -
                2 * (l + 2) * (dPmin[l - 1] - dPmax[l - 1])
                    ) / (xmin[i] - xmax[i]);

        Glminus[i][l] =
            (2. * l + 1) / (2. * M_PI * l * l * (l + 1) * (l + 1)) *
            ( -l * (l - 1.) / 2 * (l + 2. / (2 * l + 1)) *
                    (Pmin[l - 1] - Pmax[l - 1]) -
                l * (l - 1.) * (2. - l) / 2 *
                    (xmin[i] * Pmin[l] - xmax[i] * Pmax[l]) +
                l * (l - 1) / (2 * l + 1) * (Pmin[l + 1] - Pmax[l + 1])
                + (4 - l) * (dPmin[l] - dPmax[l]) +
                (l + 2) * (xmin[i] * dPmin[l - 1] - xmax[i] * dPmax[l - 1] -
                           Pmin[l - 1] + Pmax[l - 1])
                - 2 * (l - 1) *
                      (xmin[i] * dPmin[l] - xmax[i] * dPmax[l] - Pmin[l] +
                       Pmax[l]) +
                2 * (l + 2) * (dPmin[l - 1] - dPmax[l - 1])) /
            (xmin[i] - xmax[i]);
      }

      free_double_vector(Pmin, 0, LMAX + 1);
      free_double_vector(Pmax, 0, LMAX + 1);
      free_double_vector(dPmin, 0, LMAX + 1);
      free_double_vector(dPmax, 0, LMAX + 1);
    }
  }
  if (recompute_shear(C, N))
  {
    {
      int nz = 0;
      double Cl_EE[LMAX];
      double Cl_BB[LMAX];
      const double Z1NZ = Z1(nz);
      const double Z2NZ = Z2(nz);
      // only compute BB if the TATT parameters allow for B-mode terms
      if ( nuisance.b_ta_z[0] ||
           nuisance.b_ta_z[Z1(nz)] ||
           nuisance.b_ta_z[Z2(nz)] ||
           nuisance.A2_ia ||
           nuisance.A2_z[Z1(nz)] ||
           nuisance.A2_z[Z2(nz)])
      {
        {
          const int l = 2;
          Cl_EE[l] = C_EE_TATT(1.0 * l, Z1NZ, Z2NZ);
          Cl_BB[l] = C_BB_TATT(1.0 * l, Z1NZ, Z2NZ);
        }
        #pragma omp parallel for
        for (int l = 3; l < LMIN_tab; l++) {
          Cl_EE[l] = C_EE_TATT(1.0 * l, Z1NZ, Z2NZ);
          Cl_BB[l] = C_BB_TATT(1.0 * l, Z1NZ, Z2NZ);
        }
        {
          const int l = LMIN_tab;
          Cl_EE[l] = C_EE_tab(1.0 * l, Z1NZ, Z2NZ);
          Cl_BB[l] = C_BB_tab(1.0 * l, Z1NZ, Z2NZ);
        }
        #pragma omp parallel for
        for (int l = LMIN_tab + 1; l < LMAX; l++) {
          Cl_EE[l] = C_EE_tab(1.0 * l, Z1NZ, Z2NZ);
          Cl_BB[l] = C_BB_tab(1.0 * l, Z1NZ, Z2NZ);
        }
      } else
      {
        {
          const int l = 2;
          Cl_EE[l] = C_EE_TATT(1.0 * l, Z1NZ, Z2NZ);
          Cl_BB[l] = 0.0;
        }
        #pragma omp parallel for
        for (int l = 3; l < LMIN_tab; l++) {
          Cl_EE[l] = C_EE_TATT(1.0 * l, Z1NZ, Z2NZ);
          Cl_BB[l] = 0.0;
        }
        {
          const int l = LMIN_tab;
          Cl_EE[l] = C_EE_tab(1.0 * l, Z1NZ, Z2NZ);
          Cl_BB[l] = 0.0;
        }
        #pragma omp parallel for
        for (int l = LMIN_tab+1; l < LMAX; l++) {
          Cl_EE[l] = C_EE_tab(1.0 * l, Z1NZ, Z2NZ);
          Cl_BB[l] = 0.0;
        }
      }
      #pragma omp parallel for
      for (int i = 0; i < like.Ntheta; i++) {
        xi_vec_plus[nz * like.Ntheta + i] = 0;
        xi_vec_minus[nz * like.Ntheta + i] = 0;
        for (int l = 2; l < LMAX; l++) {
          xi_vec_plus[nz * like.Ntheta + i] +=
              Glplus[i][l] * (Cl_EE[l] + Cl_BB[l]);
          xi_vec_minus[nz * like.Ntheta + i] +=
              Glminus[i][l] * (Cl_EE[l] - Cl_BB[l]);
        }
      }
    }
    #pragma omp parallel for
    for (int nz = 1; nz < tomo.shear_Npowerspectra; nz++) {
      double *Cl_EE = calloc(LMAX, sizeof(double));
      double *Cl_BB = calloc(LMAX, sizeof(double));
      const double Z1NZ = Z1(nz);
      const double Z2NZ = Z2(nz);
      // only compute BB if the TATT parameters allow for B-mode terms
      if (nuisance.b_ta_z[0] || nuisance.b_ta_z[Z1(nz)] ||
          nuisance.b_ta_z[Z2(nz)] || nuisance.A2_ia || nuisance.A2_z[Z1(nz)] ||
          nuisance.A2_z[Z2(nz)]) {
        for (int l = 2; l < LMIN_tab; l++) {
          Cl_EE[l] = C_EE_TATT(1.0 * l, Z1NZ, Z2NZ);
          Cl_BB[l] = C_BB_TATT(1.0 * l, Z1NZ, Z2NZ);
        }
        for (int l = LMIN_tab; l < LMAX; l++) {
          Cl_EE[l] = C_EE_tab(1.0 * l, Z1NZ, Z2NZ);
          Cl_BB[l] = C_BB_tab(1.0 * l, Z1NZ, Z2NZ);
        }
      }
      else {
        for (int l = 2; l < LMIN_tab; l++) {
          Cl_EE[l] = C_EE_TATT(1.0 * l, Z1NZ, Z2NZ);
          Cl_BB[l] = 0.0;
        }
        for (int l = LMIN_tab; l < LMAX; l++) {
          Cl_EE[l] = C_EE_tab(1.0 * l, Z1NZ, Z2NZ);
          Cl_BB[l] = 0.0;
        }
      }
      for (int i = 0; i < like.Ntheta; i++) {
        xi_vec_plus[nz * like.Ntheta + i] = 0;
        xi_vec_minus[nz * like.Ntheta + i] = 0;
        for (int l = 2; l < LMAX; l++) {
          xi_vec_plus[nz * like.Ntheta + i] +=
              Glplus[i][l] * (Cl_EE[l] + Cl_BB[l]);
          xi_vec_minus[nz * like.Ntheta + i] +=
              Glminus[i][l] * (Cl_EE[l] - Cl_BB[l]);
        }
      }
      free(Cl_EE);
      free(Cl_BB);
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }

  if (pm > 0) {
    return xi_vec_plus[N_shear(ni, nj) * like.Ntheta + nt];
  } else {
    return xi_vec_minus[N_shear(ni, nj) * like.Ntheta + nt];
  }
}

double w_gamma_t_TATT(int nt, int ni, int nj) {
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
  if (Pl == 0) {
    Pl = create_double_matrix(0, like.Ntheta - 1, 0, LMAX - 1);
    w_vec = create_double_vector(0, tomo.ggl_Npowerspectra * like.Ntheta - 1);
    NTHETA = like.Ntheta;
    double xmin[like.Ntheta];
    double xmax[like.Ntheta];
    const double logdt = (log(like.vtmax) - log(like.vtmin)) / like.Ntheta;
    for (int i = 0; i < like.Ntheta; i++) {
      xmin[i] = cos(exp(log(like.vtmin) + (i + 0.0) * logdt));
      xmax[i] = cos(exp(log(like.vtmin) + (i + 1.0) * logdt));
      log_info("bin %d: theta_min = %e [rad], theta_max = %e [rad]", i,
          exp(log(like.vtmin) + (i + 0.0) * logdt),
          exp(log(like.vtmin) + (i + 1.0) * logdt));
    }
    #pragma omp parallel for
    for (int i = 0; i < NTHETA; i++) {
      double Pmin[LMAX + 1];
      double Pmax[LMAX + 1];

      gsl_sf_legendre_Pl_array(LMAX, xmin[i], Pmin);
      gsl_sf_legendre_Pl_array(LMAX, xmax[i], Pmax);

      for (int l = 1; l < LMAX; l++) {
        Pl[i][l] = (2. * l + 1) /
                   (4. * M_PI * l * (l + 1) * (xmin[i] - xmax[i])) *
                   ((l + 2. / (2 * l + 1.)) * (Pmin[l - 1] - Pmax[l - 1]) +
                    (2 - l) * (xmin[i] * Pmin[l] - xmax[i] * Pmax[l]) -
                    2. / (2 * l + 1.) * (Pmin[l + 1] - Pmax[l + 1]));
      }
    }
  }
  if (recompute_ggl(C, G, N, ni)) {
    // COCOA: CANT OPENMP HERE - VARYNG ONLY BIAS CAUSES RACE CONDITION
    for (int nz = 0; nz < tomo.ggl_Npowerspectra; nz++)
    {
      double *Cl = calloc(LMAX, sizeof(double));
      {
        const int l = 2;
        Cl[l] = C_ggl_TATT(1.0 * l, ZL(nz), ZS(nz));
      }
      #pragma omp parallel for
      for (int l = 3; l < LMIN_tab; l++) {
        Cl[l] = C_ggl_TATT(1.0 * l, ZL(nz), ZS(nz));
      }
      {
        const int l = LMIN_tab;
        Cl[l] = C_ggl_TATT_tab(1.0 * l, ZL(nz), ZS(nz));
      }
      #pragma omp parallel for
      for (int l = LMIN_tab+1; l < LMAX; l++) {
        Cl[l] = C_ggl_TATT_tab(1.0 * l, ZL(nz), ZS(nz));
      }
      #pragma omp parallel for
      for (int i = 0; i < NTHETA; i++) {
        w_vec[nz * like.Ntheta + i] = 0;
        for (int l = 2; l < LMAX; l++) {
          w_vec[nz * like.Ntheta + i] += Pl[i][l] * Cl[l];
        }
      }
      free(Cl);
    }
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  return w_vec[N_ggl(ni, nj) * like.Ntheta + nt];
}
