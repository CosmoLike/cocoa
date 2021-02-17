#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../cfftlog/cfftlog.h"

#include "bias.h"
#include "basics.h"
#include "cosmo3D.h"
#include "cosmo2D.h"
#include "IA.h"
#include "pt_cfastpt.h"
#include "radial_weights.h"
#include "recompute.h"
#include "redshift_spline.h"
#include "structs.h"

#include "log.c/src/log.h"

static int LMAX = 100000;
static int LMIN_tab = 20;

typedef double (*C_tomo_pointer)(double l, int n1, int n2);

// ----------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Limber Approximation
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// TODO: MAKE THIS MORE GENERAL
double beam_planck(double l) {
  double fwhm_arcmin = 5.4;
  double sigma = fwhm_arcmin/sqrt(8.*log(2.0))*constants.arcmin;
  return exp(-0.5*l*l*sigma*sigma);
}

// ----------------------------------------------------------------------------
// Integrands for Angular Power Spectra
// ----------------------------------------------------------------------------

double int_for_C_shear_tomo(double a, void *params) {
  if (a >= 1.0) {
    log_fatal("a>=1");
    exit(1);
  }
  double *ar = (double *) params;

  const double ell_prefactor = (ar[2] - 1.)*(ar[2])*(ar[2] + 1.)*(ar[2] + 2.);
  struct chis chidchi = chi_all(a);
  const double ell = ar[2]+0.5;
  const double fK = f_K(chi(a));
  const double k = ell/fK;

  const double res =
    W_kappa(a, fK, ar[0])*W_kappa(a, fK, ar[1])*chidchi.dchida/(fK*fK);
  return res*Pdelta(k,a);
}

double int_for_C_gl_tomo(double a, void *params) { // Add RSD
  if (a >= 1.0) {
    log_fatal("a>=1");
    exit(1);
  }
  double *ar = (double *) params;

  const double ell_prefactor1 = (ar[2])*(ar[2] + 1.);
  double ell_prefactor2 = (ar[2] - 1.)*ell_prefactor1*(ar[2] + 2.);
  if(ell_prefactor2 <= 0.) {
    ell_prefactor2 = 0.;
  }
  else {
    ell_prefactor2 = sqrt(ell_prefactor2);
  }

  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double ell = ar[2]+0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double wgal = W_gal(a, ar[0], chidchi.chi, hoverh0) +
    W_mag(a, fK, ar[0])*(ell_prefactor1/ell/ell - 1.);

  const double res = wgal*W_kappa(a, fK, ar[1])*chidchi.dchida/(fK * fK);
  return res*Pdelta(k, a)*ell_prefactor2/(ell * ell);
}

double int_for_C_gl_tomo_b2(double a, void *params) {
  if (a >= 1.0) {
    log_fatal("a>=1");
    exit(1);
  }
  double *ar = (double *) params;

  const double b1 = gbias.b1_function(1./a - 1., (int)ar[0]);
  const double b2 = gbias.b2[(int)ar[0]];
  const double bs2 = gbias.bs2[(int)ar[0]];
  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double g4 = growfac_a*growfac_a*growfac_a*growfac_a;

  const double ell = ar[2]+0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;

  double res =
    W_HOD(a, ar[0], hoverh0)*W_kappa(a, fK, ar[1])*chidchi.dchida/(fK*fK);

  res *= b1*Pdelta(k,a) + g4*(0.5*b2*PT_d1d2(k) + 0.5*bs2*PT_d1s2(k) +
    0.5*b3nl_from_b1(b1)*PT_d1d3(k));

  res += (W_mag(a, fK, ar[0])*W_kappa(a, fK, ar[1])*chidchi.dchida/(fK*fK))*
    b1*Pdelta(k,a);

  return res;
}

double int_for_C_cl_tomo(double a, void *params) {
  if (a >= 1.0) {
    log_fatal("a>=1");
    exit(1);
  }
  double *ar = (double *) params;
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
  if (a >= 1.0) {
    log_fatal("a>=1 in int_for_C_cl_tomo_b2");
    exit(1);
  }

  double *ar = (double *) params;
  const double b1 = gbias.b1_function(1. / a - 1., (int)ar[0]);
  const double b2 = gbias.b2[(int) ar[0]];
  const double bs2 = gbias.bs2[(int) ar[0]];
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

// galaxy position x kappa CMB
double int_for_C_gk(double a, void *params) {
  if (a >= 1.0) {
    log_fatal("a>=1");
    exit(1);
  }
  double *ar = (double *) params;

  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double ell = ar[1]+0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;

  const double res = W_gal(a, ar[0], chidchi.chi, hoverh0)/(fK*fK);
  return res*W_k(a,fK)*chidchi.dchida*Pdelta(k,a);
}

double int_for_C_gk_b2(double a, void *params) {
  if (a >= 1.0) {
    log_fatal("a>=1");
    exit(1);
  }
  double *ar = (double *) params;

  const double b1 = gbias.b1_function(1./a-1.,(int)ar[0]);
  const double b2 = gbias.b2[(int)ar[0]];
  const double bs2 = gbias.bs2[(int)ar[0]];
  const double growfac_a = growfac(a);
  const double g4 = growfac_a*growfac_a*growfac_a*growfac_a;

  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double ell = ar[1]+0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;

  const double res = W_HOD(a,ar[0],hoverh0)*W_k(a,fK)*chidchi.dchida/(fK*fK);
  return res*(b1*Pdelta(k,a)+g4*(0.5*b2*PT_d1d2(k)+0.5*bs2*PT_d1s2(k)));
}

// kappa CMB x kappaCMB
double int_for_C_kk(double a, void *params) {
  if (a >= 1.0) {
    log_fatal("a>=1");
    exit(1);
  }
  double *ar = (double *) params;

  struct chis chidchi = chi_all(a);
  const double ell = ar[0] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;

  const double res = pow(W_k(a,fK), 2)*chidchi.dchida/(fK*fK);
  return res*Pdelta(k,a);
}

// shear x kappa CMB
double int_for_C_ks(double a, void *params) {
  if (a >= 1.0) {
    log_fatal("a>=1");
    exit(1);
  }
  double *ar = (double *) params;

  struct chis chidchi = chi_all(a);
  const double ell = ar[1] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;

  const double res = W_kappa(a, fK, ar[0])/(fK*fK);
  return res*W_k(a, fK)*chidchi.dchida*Pdelta(k,a);
}

double int_for_C_ks_IA_mpp(double a, void *params) { // for like.IA == 4
  double *ar = (double *) params;
  struct chis chidchi = chi_all(a);

  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double ell = ar[1]+0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double ws1 = W_source(a, ar[0], hoverh0);
  const double wk1 = W_kappa(a, fK, ar[0]);
  const double wk2 = W_k(a, fK);
  const double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac(a)*
  nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia), nuisance.eta_ia);

  const double res= -ws1*wk2*norm + wk1*wk2;
  return res*Pdelta(k,a)*chidchi.dchida/(fK*fK);
}

// ----------------------------------------------------------------------------
// Angular Power Spectra - no interpolation
// ----------------------------------------------------------------------------

//shear tomography power spectra of source galaxy bins ni, nj
double C_shear_tomo_nointerp(double l, int ni, int nj) {
  double array[3] = {(double) ni, (double) nj, l};
  int j,k;
  if (ni <= nj) {
    j = nj;
    k = ni;
  }
  else{
    j = ni;
    k = nj;
  }
  return int_gsl_integrate_medium_precision(
    int_for_C_shear_tomo,
    (void*) array,
    amin_source(j),
    0.99999,
    NULL,
    1000
  );
}

//G-G lensing power spectrum, lens bin ni, source bin nj
double C_gl_tomo_nointerp(double l, int ni, int nj)  {
  if(l == 1.) {
    return 0.;
  }
  double array[3] = {(double) ni, (double) nj, l};
  if (gbias.b2[ni] || gbias.b2[nj]) {
    return int_gsl_integrate_low_precision(
      int_for_C_gl_tomo_b2,
      (void*) array,
      amin_lens(ni),
      amax_lens(ni),
      NULL,
      1000
    );
  }
  return int_gsl_integrate_medium_precision(
    int_for_C_gl_tomo,
    (void*) array,
    amin_lens(ni),
    0.99999,
    NULL,
    1000
  );
}

// galaxy clustering power spectrum bins ni, nj
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
    // COCOA: CHANGED INTEGRATION FROM MEDIUM TO LOW PRECISION
    return int_gsl_integrate_low_precision(
      int_for_C_cl_tomo,
      (void *) array,
      amin_lens(nj),
      0.99999,
      NULL,
      1000
    ); // zi<=zj
  }
}

// galaxy position x kappa CMB, lens z-bin nl
double C_gk_tomo_nointerp(double l, int nl) {
  double array[2] = {(double)nl,l};
  if (gbias.b2[nl] || gbias.b2[nl]) {
    return int_gsl_integrate_medium_precision(
      int_for_C_gk_b2,
      (void*) array,
      amin_lens(nl),
      amax_lens(nl),
      NULL,
      1000
    );
  }
  return int_gsl_integrate_medium_precision(
    int_for_C_gk,
    (void*) array,
    amin_lens(nl),
    0.99999,
    NULL,
    1000
  );
}

double C_ks_IA_nointerp(double s, int ni) {
  double array[2] = {(double) ni, s};
  if (like.IA == 4) {
    return int_gsl_integrate_medium_precision(
      int_for_C_ks_IA_mpp,
      (void*) array,
      amin_source(ni),
      amax_source(ni),
      NULL,
      1000
    );
  }
  log_fatal("C_ks_IA does not support like.IA = %d", like.IA);
  exit(1);
}

// shear x kappa CMB, for source z-bin ns
double C_ks_tomo_nointerp(double l, int ns) {
  if (like.IA) {
    return C_ks_IA_nointerp(l,ns);
  }
  double array[2] = {(double) ns, l};
  return int_gsl_integrate_medium_precision(
    int_for_C_ks,
    (void*) array,
    amin_source(ns),
    amax_source(ns),
    NULL,
    1000
  );
}

// kappa CMB x kappa CMB
double C_kk_nointerp(double l) {
  double array[1] = {l};
  return int_gsl_integrate_medium_precision(
    int_for_C_kk,
    (void*) array,
    limits.a_min*(1.+1.e-5),
    1. - 1.e-5,
    NULL,
    1000
  );
}

// ----------------------------------------------------------------------------
// Angular Power Spectra - with interpolation
// ----------------------------------------------------------------------------

double C_shear_tomo(double l, int ni, int nj) {
  static cosmopara C;
  static nuisancepara N;
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  if (ni < 0 || ni >= tomo.shear_Nbin || nj < 0 || nj >= tomo.shear_Nbin) {
    log_fatal("C_shear_tomo(l,%d,%d) outside tomo.shear_Nbin range", ni, nj);
    exit(1);
  }
  if (recompute_shear(C,N)) {
    if (table == 0) {
      table = create_double_matrix(0, tomo.shear_Npowerspectra-1, 0,
        Ntable.N_ell - 1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell-1.);
    }
    {
      const int k = 0;
      {
        const int i = 0;
        const double llog = logsmin + i*ds;
        table[k][i]= log(C_shear_tomo_nointerp(exp(llog), Z1(k), Z2(k)));
      }
      #pragma omp parallel for
      for (int i=1; i<Ntable.N_ell; i++) {
        const double llog = logsmin + i*ds;
        table[k][i]= log(C_shear_tomo_nointerp(exp(llog), Z1(k), Z2(k)));
      }
    }
    #pragma omp parallel for
    for (int k=1; k<tomo.shear_Npowerspectra; k++) {
      for (int i=0; i<Ntable.N_ell; i++) {
        const double llog = logsmin + i*ds;
        table[k][i]= log(C_shear_tomo_nointerp(exp(llog), Z1(k), Z2(k)));
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  double f1 = exp(
    interpol_fitslope(
      table[N_shear(ni,nj)],
      Ntable.N_ell,
      logsmin,
      logsmax,
      ds,
      log(l),
      1.
    )
  );
  if (isnan(f1)) {
    f1 = 0.;
  }
  return f1;
}

double C_gl_tomo(double l, int ni, int nj) {
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  if (ni < 0 || ni >= tomo.clustering_Nbin || nj < 0 || nj >= tomo.shear_Nbin) {
    log_fatal("C_gl_tomo(l,%d,%d) outside tomo.X_Nbin range", ni, nj);
    exit(1);
  }
  if (recompute_ggl(C,G,N,ni)) {
    if(limits.P_2_s_min == 1) { // COCOA: added this logical test
      log_fatal("dev (logical) error: interpolation starting w/ l=1"
        " creates a race condition");
      exit(1);
    }
    if (table == 0) {
      table =
        create_double_matrix(0, tomo.ggl_Npowerspectra-1, 0, Ntable.N_ell - 1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell - 1.);
    }
    {
      const int k = 0;
      {
        const int i = 0;
        const double llog = logsmin + i*ds;
        table[k][i]= log(C_gl_tomo_nointerp(exp(llog), ZL(k), ZS(k)));
      }
      #pragma omp parallel for
      for (int i=1; i<Ntable.N_ell; i++) {
          const double llog = logsmin + i*ds;
          table[k][i]= log(C_gl_tomo_nointerp(exp(llog), ZL(k), ZS(k)));
      }
    }
    #pragma omp parallel for
    for (int k=1; k<tomo.ggl_Npowerspectra; k++) {
      for (int i=0; i<Ntable.N_ell; i++) {
          const double llog = logsmin + i*ds;
          table[k][i]= log(C_gl_tomo_nointerp(exp(llog), ZL(k), ZS(k)));
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
    update_galpara(&G);
  }
  if(test_zoverlap(ni, nj)) {
    double f1 = exp(
      interpol_fitslope(
        table[N_ggl(ni,nj)],
        Ntable.N_ell,
        logsmin,
        logsmax,
        ds,
        log(l),
        1.
      )
    );
    if (isnan(f1)) {
      f1 = 0;
    }
    return f1;
  } else {
    return 0.0;
  }
}

double C_cl_tomo(double l, int ni, int nj) {
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  if (ni < 0 || ni >= tomo.clustering_Nbin || nj < 0 ||
    nj >= tomo.clustering_Nbin) {
    log_fatal("C_cl_tomo(l,%d,%d) outside tomo.clustering_Nbin range", ni, nj);
    exit(1);
  }
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

  const int j = ni * tomo.clustering_Nbin + nj;

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
  double f1 = exp(
    interpol_fitslope(
      table[j],
      Ntable.N_ell,
      logsmin,
      logsmax,
      ds,
      log(l),
      1.
    )
  );
  if (isnan(f1)) {
    f1 = 0.;
  }
  return f1;
}

double C_gk_tomo(double l, int ni) {
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  if (ni < 0 || ni >= tomo.clustering_Nbin) {
    log_fatal("Bin %d outside tomo.clustering_Nbin range", ni);
    exit(1);
  }
  if (recompute_clustering(C,G,N,ni,ni)) {
    if (table == 0) {
      table = create_double_matrix(0, tomo.clustering_Nbin-1, 0,
        Ntable.N_ell-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell);
    }
    {
      const int k = 0;
      table[k][0] = log(C_gk_tomo_nointerp(exp(logsmin),k));
      #pragma omp parallel for
      for (int i = 1; i < Ntable.N_ell; i++) {
        const double llog = logsmin + i*ds;
        table[k][i] = log(C_gk_tomo_nointerp(exp(llog),k));
      }
    }
    #pragma omp parallel for
    for (int k=1; k<tomo.clustering_Nbin; k++) {
      for (int i=0; i<Ntable.N_ell; i++) {
        const double llog = logsmin + i*ds;
        table[k][i]= log(C_gk_tomo_nointerp(exp(llog),k));
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
    update_galpara(&G);
  }
  double f1 = exp(interpol(table[ni], Ntable.N_ell, logsmin,
    logsmax, ds, log(l), 1., 1.));
  if (isnan(f1)) {
    f1 = 0.;
  }
  return f1;
}

double C_ks_tomo(double l, int ni) {
  static cosmopara C;
  static nuisancepara N;
  static double **table, *sig;
  static int osc[100];
  static double ds = .0, logsmin = .0, logsmax = .0;
  if (ni < 0 || ni >= tomo.clustering_Nbin){
    log_fatal("Bin %d outside tomo.clustering_Nbin range", ni);
    exit(1);
  }
  if (recompute_shear(C, N)) {
    if (table == 0) {
      table = create_double_matrix(0, tomo.shear_Nbin-1, 0, Ntable.N_ell-1);
      sig = create_double_vector(0,tomo.shear_Nbin-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell);
    }
    {
      const int k = 0;
      sig[k] = 1.;
      osc[k] = 0;
      const double res = C_ks_tomo_nointerp(500.,k);
      if (res < 0){
        sig[k] = -1.;
      }
      #pragma omp parallel for
      for (int i=0; i<Ntable.N_ell; i++) {
        const double llog = logsmin + i*ds;
        table[k][i]= C_ks_tomo_nointerp(exp(llog), k);
        if (res*sig[k] < 0.) {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0){
        for(int i = 0; i < Ntable.N_ell; i++){
          table[k][i] = log(sig[k]*table[k][i]);
        }
      }
    }
    #pragma omp parallel for
    for (int k=1; k<tomo.shear_Nbin; k++) {
      sig[k] = 1.;
      osc[k] = 0;
      const double res = C_ks_tomo_nointerp(500.,k);
      if (res < 0){
        sig[k] = -1.;
      }
      for (int i=0; i<Ntable.N_ell; i++) {
        const double llog = logsmin + i*ds;
        table[k][i]= C_ks_tomo_nointerp(exp(llog), k);
        if (res*sig[k] < 0.) {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0){
        for(int i = 0; i < Ntable.N_ell; i++){
          table[k][i] = log(sig[k]*table[k][i]);
        }
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  double f1 = 0.;
  if (osc[ni] == 0) {
    f1 = sig[ni]*exp(interpol_fitslope(table[ni], Ntable.N_ell, logsmin,
      logsmax, ds, log(l), 1.));
  }
  if (osc[ni] == 1) {
    f1 = interpol_fitslope(table[ni], Ntable.N_ell, logsmin, logsmax, ds,
      log(l), 1.);
  }
  if (isnan(f1)) {
    f1 = 0.;
  }
  return f1;
}

double C_kk(double l) {
  static cosmopara C;
  static double *table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  if (recompute_cosmo3D(C)) {
    if (table == 0) {
       table   = create_double_vector(0, Ntable.N_ell - 1);
       logsmin = log(limits.P_2_s_min);
       logsmax = log(limits.P_2_s_max);
       ds = (logsmax - logsmin)/(Ntable.N_ell);
    }
    {
      const int i = 0;
      const double llog = logsmin + i*ds;
       table[i]= log(C_kk_nointerp(exp(llog)));
    }
    #pragma omp parallel for
    for (int i=1; i<Ntable.N_ell; i++) {
      const double llog = logsmin + i*ds;
       table[i] = log(C_kk_nointerp(exp(llog)));
    }
    update_cosmopara(&C);
  }
  double f1 = exp(interpol(table, Ntable.N_ell, logsmin, logsmax,
    ds, log(l), 1., 1.));
  if (isnan(f1)) {
    f1 = 0.;
  }
  return f1;
}

double C_gk_wrapper(double l, int ni) {
  return C_gk_tomo(l, ni)*beam_planck(l);
}

double C_ks_wrapper(double l, int ni) {
  return C_ks_tomo(l, ni)*beam_planck(l);
}

// ----------------------------------------------------------------------------
// Correlation Functions (real space) - flat sky
// ----------------------------------------------------------------------------

void hankel_kernel_FT(double x, fftw_complex *res, double *arg, int argc __attribute__((unused))) {
  fftw_complex a1, a2, g1, g2;
  int mu;
  double mod, xln2, si, co, d1, d2, pref, q;
  q = arg[0];
  mu = (int)(arg[1] + 0.1);

  /* arguments for complex gamma */
  a1[0] = 0.5 * (1.0 + mu + q);
  a2[0] = 0.5 * (1.0 + mu - q);
  a1[1] = 0.5 * x;
  a2[1] = -a1[1];
  cdgamma(a1, &g1);
  cdgamma(a2, &g2);
  xln2 = x * constants.ln2;
  si = sin(xln2);
  co = cos(xln2);
  d1 = g1[0] * g2[0] + g1[1] * g2[1]; /* Re */
  d2 = g1[1] * g2[0] - g1[0] * g2[1]; /* Im */
  mod = g2[0] * g2[0] + g2[1] * g2[1];
  pref = exp(constants.ln2 * q) / mod;

  (*res)[0] = pref * (co * d1 - si * d2);
  (*res)[1] = pref * (si * d1 + co * d2);
}

void xipm_via_hankel(double **xi, double *logthetamin, double *logthetamax,
C_tomo_pointer C_tomo,int ni, int nj) {
  const double l_min = 0.0001;
  const double l_max = 5.0e6;
  static double loglmax = -123.0, loglmin, dlnl,  lnrc, arg[2];
  static int nc;

  double        l, kk, *lP, t;
  fftw_plan     plan1,plan;
  fftw_complex *f_lP,*conv;
  fftw_complex  kernel;
  int           i, count;
  lP   = fftw_malloc(Ntable.N_thetaH*sizeof(double));
  f_lP = fftw_malloc((Ntable.N_thetaH/2+1)*sizeof(fftw_complex));
  conv = fftw_malloc((Ntable.N_thetaH/2+1)*sizeof(fftw_complex));
  plan  = fftw_plan_dft_r2c_1d(Ntable.N_thetaH, lP, f_lP, FFTW_ESTIMATE);
  plan1 = fftw_plan_dft_c2r_1d(Ntable.N_thetaH, conv, lP, FFTW_ESTIMATE);
  if (loglmax==-123.0) {
    loglmax  = log(l_max);
    loglmin  = log(l_min);
    dlnl     = (loglmax-loglmin)/(1.0*Ntable.N_thetaH);
    lnrc     = 0.5*(loglmax+loglmin);
    nc       = Ntable.N_thetaH/2+1;
  }
  /* Power spectrum on logarithmic bins */
  for(i=0; i<Ntable.N_thetaH; i++) {
    l     = exp(lnrc+(i-nc)*dlnl);
    lP[i] = l*C_tomo(l,ni,nj);
  }
  /* go to log-Fourier-space */
  fftw_execute(plan);
  arg[0] = 0;   /* bias */
  for (count=0; count<=1; count++) {
    arg[1] = (count==0 ? 0 : 4);   /* order of Bessel function */
    /* perform the convolution, negative sign for kernel (complex conj.!) */
    for(i=0; i<Ntable.N_thetaH/2+1; i++) {
      kk = 2*constants.pi*i/(dlnl*Ntable.N_thetaH);
      hankel_kernel_FT(kk, &kernel, arg, 2);
      conv[i][0] = f_lP[i][0]*kernel[0]-f_lP[i][1]*kernel[1];
      conv[i][1] = f_lP[i][1]*kernel[0]+f_lP[i][0]*kernel[1];
    }
    /* force Nyquist- and 0-frequency-components to be double */
    conv[0][1] = 0;
    conv[Ntable.N_thetaH/2][1] = 0;
    /* go back to double space, i labels log-bins in theta */
    fftw_execute(plan1);
    for(i=0; i<Ntable.N_thetaH; i++) {
      t = exp((nc-i)*dlnl-lnrc);             /* theta=1/l */
      xi[count][Ntable.N_thetaH-i-1] = lP[i]/(t*2*constants.pi*Ntable.N_thetaH);
    }
  }

  *logthetamin = (nc-Ntable.N_thetaH+1)*dlnl-lnrc;
  *logthetamax = nc*dlnl-lnrc;
  /* clean up */
  fftw_free(conv);
  fftw_free(lP);
  fftw_free(f_lP);
  fftw_destroy_plan(plan);
  fftw_destroy_plan(plan1);
}

void twopoint_via_hankel(double **xi, double *logthetamin, double *logthetamax,
C_tomo_pointer C_tomo, int ni, int nj, int N_Bessel) {
  const double l_min = 0.0001;
  const double l_max = 5.0e6;
  double loglmax, loglmin, dlnl, lnrc, arg[2];
  static int nc;

  double        l, kk, *lP, t;
  fftw_plan     plan1,plan;
  fftw_complex *f_lP,*conv;
  fftw_complex  kernel;
  int           i;
  lP   = fftw_malloc(Ntable.N_thetaH*sizeof(double));
  f_lP = fftw_malloc((Ntable.N_thetaH/2+1)*sizeof(fftw_complex));
  conv = fftw_malloc((Ntable.N_thetaH/2+1)*sizeof(fftw_complex));
  plan  = fftw_plan_dft_r2c_1d(Ntable.N_thetaH, lP, f_lP, FFTW_ESTIMATE);
  plan1 = fftw_plan_dft_c2r_1d(Ntable.N_thetaH, conv, lP, FFTW_ESTIMATE);
  loglmax  = log(l_max);
  loglmin  = log(l_min);
  dlnl     = (loglmax-loglmin)/(1.0*Ntable.N_thetaH-1.);
  lnrc     = 0.5*(loglmax+loglmin);
  nc       = Ntable.N_thetaH/2+1;
  /* Power spectrum on logarithmic bins */
  for(i=0; i<Ntable.N_thetaH; i++) {
    l     = exp(lnrc+(i-nc)*dlnl);
    lP[i] = l*C_tomo(l,ni,nj);

  }

  /* go to log-Fourier-space */
  fftw_execute(plan);
  arg[0] = 0;   /* bias */
  arg[1] = N_Bessel;   /* order of Bessel function */
  /* perform the convolution, negative sign for kernel (complex conj.!) */
  for(i=0; i<Ntable.N_thetaH/2+1; i++) {
    kk = 2*constants.pi*i/(dlnl*Ntable.N_thetaH);
    hankel_kernel_FT(kk, &kernel, arg, 2);
    conv[i][0] = f_lP[i][0]*kernel[0]-f_lP[i][1]*kernel[1];
    conv[i][1] = f_lP[i][1]*kernel[0]+f_lP[i][0]*kernel[1];
  }
  /* force Nyquist- and 0-frequency-components to be double */
  conv[0][1] = 0;
  conv[Ntable.N_thetaH/2][1] = 0;
  /* go back to double space, i labels log-bins in theta */
  fftw_execute(plan1);
  for(i=0; i<Ntable.N_thetaH; i++) {
    t = exp((nc-i)*dlnl-lnrc);             /* theta=1/l */
    xi[0][Ntable.N_thetaH-i-1] = lP[i]/(t*2*constants.pi*Ntable.N_thetaH);
  }


  *logthetamin = (nc-Ntable.N_thetaH+1)*dlnl-lnrc;
  *logthetamax = nc*dlnl-lnrc;

  /* clean up */
  fftw_free(conv);
  fftw_free(lP);
  fftw_free(f_lP);
  fftw_destroy_plan(plan);
  fftw_destroy_plan(plan1);
}

double xi_pm_tomo(int pm, double theta, int ni, int nj) {
  static cosmopara C;
  static nuisancepara N;
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;

  if (recompute_shear(C,N)) {
    if (like.IA != 0 && like.IA != 3 && like.IA != 4) {
      log_fatal("xi_pm_tomo does not support like.IA = %d yet", like.IA);
      exit(1);
    }
    C_tomo_pointer C_pointer = &C_shear_tomo;
    if (like.IA == 3 || like.IA == 4) {
      C_pointer = &C_shear_shear_IA_tab;
    }
    if (table == 0) {
      table = create_double_matrix(0, 2*tomo.shear_Npowerspectra - 1, 0,
        Ntable.N_thetaH - 1);
    }
    // COCOA: DONT OPENMP A LOOP THAT CALLS FFTLOG
    for (int i = 0; i < tomo.shear_Npowerspectra; i++){
      double **tab;
      tab   = create_double_matrix(0, 1, 0, Ntable.N_thetaH-1);

      xipm_via_hankel(tab, &logthetamin, &logthetamax, C_pointer, Z1(i), Z2(i));

      for (int k = 0; k < Ntable.N_thetaH; k++) {
        table[2*i][k] = tab[0][k];
        table[2*i + 1][k] = tab[1][k];
      }

      free_double_matrix(tab,0, 1, 0, Ntable.N_thetaH - 1);
    }
    dlogtheta = (logthetamax - logthetamin)/((double) Ntable.N_thetaH);
    update_cosmopara(&C);
    update_nuisance(&N);
  }

  return interpol(
    table[2*N_shear(ni, nj) + (1 - pm)/2],
    Ntable.N_thetaH,
    logthetamin,
    logthetamax,
    dlogtheta,
    log(theta),
    0.0,
    0.0
  );
}

double w_gamma_t_tomo(double theta,int ni, int nj) {
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;

  if (recompute_ggl(C,G,N,ni)) {
    if (like.IA != 0 && like.IA != 3 && like.IA != 4) {
      printf("cosmo2D_real.c: w_gamma_t_tomo does not support like.IA = %d yet\nEXIT!\n",like.IA);
      exit(1);
    }
    C_tomo_pointer C_gl_pointer = &C_gl_tomo;
    if (like.IA ==3 || like.IA ==4) C_gl_pointer = &C_ggl_IA_tab;

    if (table == 0) {
      table  = create_double_matrix(0, tomo.ggl_Npowerspectra, 0,
        Ntable.N_thetaH - 1);
    }
    // COCOA: DONT OPENMP A LOOP THAT CALLS FFTLOG
    for (int i = 0; i <tomo.ggl_Npowerspectra; i++) {
      double **tab;
      tab   = create_double_matrix(0, 1, 0, Ntable.N_thetaH-1);

      twopoint_via_hankel(tab, &logthetamin, &logthetamax,C_gl_pointer, ZL(i),ZS(i),2);

      for (int k = 0; k < Ntable.N_thetaH; k++) {
        table[i][k] = tab[0][k];
      }
      free_double_matrix(tab,0, 1, 0, Ntable.N_thetaH);
    }
    dlogtheta = (logthetamax-logthetamin)/((double)Ntable.N_thetaH);
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  if(test_zoverlap(ni,nj)) {
    return interpol(
      table[N_ggl(ni,nj)],
      Ntable.N_thetaH,
      logthetamin,
      logthetamax,
      dlogtheta,
      log(theta),
      0.0,
      0.0
    );
  } else {
    return 0.0;
  }
}

double w_tomo(double theta, int ni, int nj) {
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  if (recompute_clustering(C,G,N,ni,nj)) {
    double **tab;
    int i, j,k;
    tab   = create_double_matrix(0, 1, 0, Ntable.N_thetaH-1);
    if (table==0) table = create_double_matrix(0, tomo.clustering_Nbin*tomo.clustering_Nbin-1, 0, Ntable.N_thetaH-1);
    for (i = 0; i < tomo.clustering_Nbin; i++){
      for (j= i; j < tomo.clustering_Nbin; j++){
        twopoint_via_hankel(tab, &logthetamin, &logthetamax,&C_cl_tomo, i,j,0);
        dlogtheta = (logthetamax-logthetamin)/((double)Ntable.N_thetaH);
        for (k = 0; k < Ntable.N_thetaH; k++){
          table[i*tomo.clustering_Nbin+j][k] = tab[0][k];
          table[j*tomo.clustering_Nbin+i][k] = tab[0][k];
        }
      }
    }
    free_double_matrix(tab,0, 1, 0, Ntable.N_thetaH-1);
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  return interpol(table[ni*tomo.clustering_Nbin+nj], Ntable.N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 1.0, 1.0);
}

// ----------------------------------------------------------------------------
// Correlation Functions (real space) - full sky
// ----------------------------------------------------------------------------

double xi_pm_tomo_fullsky(int pm, int nt, int ni, int nj) {
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
  if (like.IA != 0 && like.IA != 3 && like.IA != 4) {
    log_fatal("xi_pm_tomo does not support like.IA = %d yet", like.IA);
    exit(1);
  }
  if (Glplus == 0) {
    Glplus = create_double_matrix(0, like.Ntheta - 1, 0, LMAX - 1);
    Glminus = create_double_matrix(0, like.Ntheta - 1, 0, LMAX - 1);
    xi_vec_plus = create_double_vector(0,
      tomo.shear_Npowerspectra*like.Ntheta - 1);
    xi_vec_minus = create_double_vector(0,
      tomo.shear_Npowerspectra*like.Ntheta - 1);

    double xmin[like.Ntheta];
    double xmax[like.Ntheta];
    const double logdt = (log(like.vtmax) - log(like.vtmin))/like.Ntheta;
    for(int i = 0; i < like.Ntheta ; i++) {
      xmin[i] = cos(exp(log(like.vtmin)+(i + 0.0)*logdt));
      xmax[i] = cos(exp(log(like.vtmin)+(i + 1.0)*logdt));
    }

    #pragma omp parallel for
    for (int i = 0; i<like.Ntheta; i ++){
      double *Pmin = create_double_vector(0, LMAX + 1);
      double *Pmax = create_double_vector(0, LMAX + 1);
      double *dPmin = create_double_vector(0, LMAX + 1);
      double *dPmax = create_double_vector(0, LMAX + 1);

      const double x = cos(like.theta[i]);

      gsl_sf_legendre_Pl_deriv_array(LMAX, xmin[i], Pmin, dPmin);
      gsl_sf_legendre_Pl_deriv_array(LMAX, xmax[i], Pmax, dPmax);

      for (int l = 3; l < LMAX; l ++){
        Glplus[i][l] =(2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(
          -l*(l - 1.)/2*(l+2./(2*l+1)) * (Pmin[l-1] - Pmax[l-1])
          -l*(l - 1.)*(2.-l)/2 * (xmin[i]*Pmin[l] - xmax[i]*Pmax[l])
          +l*(l - 1.)/(2.*l+1) * (Pmin[l+1] - Pmax[l+1])
          +(4 - l) * (dPmin[l] - dPmax[l])
          +(l + 2) * (xmin[i]*dPmin[l-1] - xmax[i]*dPmax[l-1] - Pmin[l-1] + Pmax[l-1])
          +2*(l - 1) * (xmin[i]*dPmin[l] - xmax[i]*dPmax[l] - Pmin[l] + Pmax[l])
          -2*(l + 2) * (dPmin[l-1] - dPmax[l-1])
          )/(xmin[i] - xmax[i]);

        Glminus[i][l] =(2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(
          -l*(l - 1.)/2*(l+2./(2*l+1)) * (Pmin[l-1] - Pmax[l-1])
          -l*(l - 1.)*(2.-l)/2 * (xmin[i]*Pmin[l] - xmax[i]*Pmax[l])
          +l*(l - 1.)/(2.*l+1) * (Pmin[l+1] - Pmax[l+1])
          + (4 - l) * (dPmin[l] - dPmax[l])
          + (l + 2) * (xmin[i]*dPmin[l-1] - xmax[i]*dPmax[l-1] - Pmin[l-1] + Pmax[l-1])
          -2*(l - 1) * (xmin[i]*dPmin[l] - xmax[i]*dPmax[l] - Pmin[l] + Pmax[l])
          +2*(l + 2) * (dPmin[l-1] - dPmax[l-1])
          )/(xmin[i]-xmax[i]);
      }

      free_double_vector(Pmin, 0, LMAX + 1);
      free_double_vector(Pmax, 0, LMAX + 1);
      free_double_vector(dPmin, 0, LMAX + 1);
      free_double_vector(dPmax, 0, LMAX + 1);
    }
  }
  if (recompute_shear(C,N)) {
    C_tomo_pointer C_pointer = &C_shear_tomo;
    if (like.IA == 3 || like.IA == 4) {
      C_pointer = &C_shear_shear_IA_tab;
    }
    // COCOA: FROM TATT EXPERIENCE, WE CANT OPENMP HERE (POSSIBLE RACE)
    for (int nz = 0; nz < tomo.shear_Npowerspectra; nz++) {
      double *Cl = calloc(LMAX, sizeof(double));
      {
        const int l = 2;
        Cl[l] = C_shear_shear_IA_tab(1.0*l, Z1(nz), Z2(nz));
      }
      #pragma omp parallel for
      for (int l = 3; l < LMAX; l++) {
        Cl[l] = C_shear_shear_IA_tab(1.0*l, Z1(nz), Z2(nz));
      }
      #pragma omp parallel for
      for (int i = 0; i < like.Ntheta; i++) {
        xi_vec_plus[nz*like.Ntheta +i] = 0;
        xi_vec_minus[nz*like.Ntheta + i] = 0;
        for (int l = 2; l < LMAX; l++) {
          xi_vec_plus[nz*like.Ntheta + i] += Glplus[i][l]*Cl[l];
          xi_vec_minus[nz*like.Ntheta + i] += Glminus[i][l]*Cl[l];
        }
      }
      free(Cl);
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  if (pm > 0) {
    return xi_vec_plus[N_shear(ni,nj)*like.Ntheta + nt];
  } else {
    return xi_vec_minus[N_shear(ni,nj)*like.Ntheta + nt];
  }
}

double w_gamma_t_tomo_fullsky(int nt, int ni, int nj) {
  static double **Pl =0;
  static double *w_vec =0;
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  if (like.Ntheta == 0) {
    log_fatal("like.Ntheta not initialized");
    exit(1);
  }
  if (like.IA != 0 && like.IA != 3 && like.IA != 4) {
    log_fatal("w_gamma_t_tomo_fullsky does not support like.IA = %d yet", like.IA);
    exit(1);
  }
  if (Pl == 0) {
    Pl = create_double_matrix(0, like.Ntheta - 1, 0, LMAX - 1);
    w_vec = create_double_vector(0,tomo.ggl_Npowerspectra*like.Ntheta - 1);

    double xmin[like.Ntheta];
    double xmax[like.Ntheta];
    const double logdt = (log(like.vtmax) - log(like.vtmin))/like.Ntheta;
    for(int i = 0; i < like.Ntheta ; i++) {
      xmin[i] = cos(exp(log(like.vtmin)+(i + 0.0)*logdt));
      xmax[i] = cos(exp(log(like.vtmin)+(i + 1.0)*logdt));
    }

    #pragma omp parallel for
    for (int i = 0; i < like.Ntheta; i++) {
      double *Pmin = create_double_vector(0, LMAX + 1);
      double *Pmax = create_double_vector(0, LMAX + 1);

      gsl_sf_legendre_Pl_array(LMAX, xmin[i], Pmin);
      gsl_sf_legendre_Pl_array(LMAX, xmax[i], Pmax);

      for (int l = 1; l < LMAX; l++) {
        Pl[i][l] = (2.*l+1)/(4.*M_PI*l*(l+1)*(xmin[i]-xmax[i]))
          *((l+2./(2*l+1.))*(Pmin[l-1]-Pmax[l-1])
          + (2-l)*(xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
          - 2./(2*l+1.)*(Pmin[l+1]-Pmax[l+1]));
      }

      free_double_vector(Pmin, 0, LMAX + 1);
      free_double_vector(Pmax, 0, LMAX + 1);
    }
  }
  if (recompute_ggl(C,G,N,ni)) {
    C_tomo_pointer C_gl_pointer = &C_gl_tomo;
    if (like.IA == 3 || like.IA == 4) {
      C_gl_pointer = &C_ggl_IA_tab;
    }
    // COCOA: FROM TATT EXPERIENCE, WE CANT OPENMP HERE (POSSIBLE RACE)
    for (int nz = 0; nz <tomo.ggl_Npowerspectra; nz++) {
      double *Cl = calloc(LMAX, sizeof(double));
      {
        const int l = 1;
        Cl[l] = C_ggl_IA_tab(1.0*l, ZL(nz), ZS(nz));
      }
      #pragma omp parallel for
      for (int l = 2; l < LMAX; l++) {
        Cl[l] = C_ggl_IA_tab(1.0*l, ZL(nz), ZS(nz));
      }
      #pragma omp parallel for
      for (int i = 0; i < like.Ntheta; i++) {
        w_vec[nz*like.Ntheta + i] = 0;
        for (int l = 2; l < LMAX; l++) {
          w_vec[nz*like.Ntheta + i] += Pl[i][l]*Cl[l];
        }
      }
      free(Cl);
    }
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  return w_vec[N_ggl(ni,nj)*like.Ntheta + nt];
}

double w_tomo_fullsky(int nt, int ni, int nj) {
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
		Pl = create_double_matrix(0, like.Ntheta - 1, 0, LMAX - 1);
		w_vec = create_double_vector(0, tomo.clustering_Nbin*like.Ntheta - 1);

		double xmin[like.Ntheta];
    double xmax[like.Ntheta];
		const double logdt = (log(like.vtmax) - log(like.vtmin))/like.Ntheta;
		for(int i = 0; i < like.Ntheta ; i++) {
			xmin[i] = cos(exp(log(like.vtmin) + (i + 0.0)*logdt));
			xmax[i] = cos(exp(log(like.vtmin) + (i + 1.0)*logdt));
		}

		#pragma omp parallel for
		for (int i = 0; i < like.Ntheta; i++) {
			double *Pmin = create_double_vector(0, LMAX + 1);
			double *Pmax = create_double_vector(0, LMAX + 1);

			gsl_sf_legendre_Pl_array(LMAX, xmin[i], Pmin);
			gsl_sf_legendre_Pl_array(LMAX, xmax[i], Pmax);

			for (int l = 1; l < LMAX; l++){
				Pl[i][l] = 1./(4.*M_PI)*(Pmin[l+1] - Pmax[l+1] - Pmin[l-1] + Pmax[l-1])/
					(xmin[i]-xmax[i]);
			}

			free_double_vector(Pmin, 0, LMAX + 1);
			free_double_vector(Pmax, 0, LMAX + 1);
		}
	}
	if (recompute_clustering(C, G, N, ni, nj)) {
		// COCOA: FROM TATT EXPERIENCE, WE CANT OPENMP HERE (POSSIBLE RACE)
		for (int nz = 0; nz < tomo.clustering_Nbin; nz++) {
			double *Cl = calloc(LMAX, sizeof(double));
			{
        int l = 1;
        Cl[l] = C_cl_tomo_nointerp(l, nz, nz);
        l = 2; // just to be on the safe side (race condition) (l=1 just zero?)
        Cl[l] = C_cl_tomo_nointerp(l, nz, nz);
      }
      #pragma omp parallel for
      for (int l = 3; l < LMIN_tab; l++) {
        Cl[l] = C_cl_tomo_nointerp(l, nz, nz);
      }
 			{
        const int l = LMIN_tab;
        C_cl_tomo(1.0*l, nz, nz);
      }
      #pragma omp parallel for
      for (int l = LMIN_tab + 1; l < LMAX; l++) {
        C_cl_tomo(1.0*l, nz, nz);
      }
			#pragma omp parallel for
			for (int i = 0; i < like.Ntheta; i++) {
				w_vec[nz*like.Ntheta + i] = 0;
				for (int l = 1; l < LMAX; l++) {
					w_vec[nz*like.Ntheta + i] += Pl[i][l]*Cl[l];
				}
			}
			free(Cl);
		}
		update_cosmopara(&C);
		update_galpara(&G);
		update_nuisance(&N);
	}
	return w_vec[ni*like.Ntheta + nt];
}

double w_gk_tomo_fullsky(int nt, int ni) {
  static double **Pl =0;
  static double *w_vec =0;
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  if (like.Ntheta == 0) {
    log_fatal("like.Ntheta not initialized");
    exit(1);
  }
  if (Pl ==0){
    Pl = create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
    w_vec = create_double_vector(0,tomo.clustering_Nbin*like.Ntheta-1);

    double xmin[like.Ntheta];
    double xmax[like.Ntheta];
    const double logdt=(log(like.vtmax)-log(like.vtmin))/like.Ntheta;
    for(int i=0; i<like.Ntheta ; i++){
      xmin[i] = cos(exp(log(like.vtmin) + (i + 0.0)*logdt));
      xmax[i] = cos(exp(log(like.vtmin) + (i + 1.0)*logdt));
    }

    #pragma omp parallel for
    for (int i = 0; i < like.Ntheta; i ++) {
      double *Pmin = create_double_vector(0, LMAX+1);
      double *Pmax = create_double_vector(0, LMAX+1);

      gsl_sf_legendre_Pl_array(LMAX, xmin[i],Pmin);
      gsl_sf_legendre_Pl_array(LMAX, xmax[i],Pmax);

      for (int l = 1; l < LMAX; l ++){
        Pl[i][l] = 1./(4.*M_PI)*(Pmin[l+1]-Pmax[l+1]-Pmin[l-1]+Pmax[l-1])/
          (xmin[i]-xmax[i]);
      }

      free_double_vector(Pmin,0,LMAX+1);
      free_double_vector(Pmax,0,LMAX+1);
    }
  }
  if (recompute_gk(C,G,N,ni)) {
    // COCOA: FROM TATT EXPERIENCE, WE CANT OPENMP HERE (POSSIBLE RACE)
    for (int nz = 0; nz <tomo.clustering_Nbin; nz ++) {
      double *Cl = calloc(LMAX, sizeof(double));
      {
        int l = 1;
        Cl[l] = C_gk_wrapper(1.0*l, nz);
      }
      #pragma omp parallel for
      for (int l = 2; l < LMAX; l++){
        Cl[l] = C_gk_wrapper(1.0*l, nz);
      }
      #pragma omp parallel for
      for (int i = 0; i < like.Ntheta; i++){
        w_vec[nz*like.Ntheta+i] =0;
        for (int l = 1; l < LMAX; l++){
          w_vec[nz*like.Ntheta+i]+=Pl[i][l]*Cl[l];
        }
      }
    }
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  return w_vec[ni*like.Ntheta+nt];
}

double w_ks_tomo_fullsky(int nt, int ni) {
  static double ** Pl =0;
  static double *w_vec =0;
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  if (like.Ntheta ==0) {
    log_fatal("like.Ntheta not initialized");
    exit(1);
  }
  if (Pl == 0) {
    Pl = create_double_matrix(0, like.Ntheta - 1, 0, LMAX - 1);
    w_vec = create_double_vector(0, tomo.shear_Nbin*like.Ntheta - 1);

    double xmin[like.Ntheta];
    double xmax[like.Ntheta];
    const double logdt = (log(like.vtmax) - log(like.vtmin))/like.Ntheta;
    for(int i = 0; i < like.Ntheta; i++) {
      xmin[i] = cos(exp(log(like.vtmin) + (i + 0.0)*logdt));
      xmax[i] = cos(exp(log(like.vtmin) + (i + 1.0)*logdt));
    }
    #pragma omp parallel for
    for (int i = 0; i < like.Ntheta; i++) {
      double *Pmin = create_double_vector(0, LMAX + 1);
      double *Pmax = create_double_vector(0, LMAX + 1);

      gsl_sf_legendre_Pl_array(LMAX, xmin[i], Pmin);
      gsl_sf_legendre_Pl_array(LMAX, xmax[i], Pmax);

      for (int l = 1; l < LMAX; l ++){
        Pl[i][l] = (2.*l+1)/(4.*M_PI*l*(l+1)*(xmin[i]-xmax[i]))
          *((l+2./(2*l+1.))*(Pmin[l-1]-Pmax[l-1])
          +(2-l)*(xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
          -2./(2*l+1.)*(Pmin[l+1]-Pmax[l+1]));
      }

      free_double_vector(Pmin, 0, LMAX+1);
      free_double_vector(Pmax, 0, LMAX+1);
    }
  }
  if (recompute_ks(C, G, N, ni)) {
    // COCOA: FROM TATT EXPERIENCE, WE CANT OPENMP HERE (POSSIBLE RACE)
    for (int nz = 0; nz < tomo.shear_Nbin; nz++) {
      double *Cl = calloc(LMAX, sizeof(double));
      {
        int l = 1;
        Cl[l] = C_ks_wrapper(1.0*l, nz);
      }
      #pragma omp parallel for
      for (int l = 2; l < LMAX; l++) {
        Cl[l] = C_ks_wrapper(1.0*l, nz);
      }
      #pragma omp parallel for
      for (int i = 0; i < like.Ntheta; i++) {
        w_vec[nz*like.Ntheta+i] = 0;
        for (int l = 2; l < LMAX; l++){
          w_vec[nz*like.Ntheta+i] += Pl[i][l]*Cl[l];
        }
      }
      free (Cl);
    }
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  return w_vec[ni*like.Ntheta + nt];
}


// ----------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Non-Limber
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Integrands for Angular Power Spectra
// ---------------------------------------------------------------------------

double int_for_C_gl_lin(double a, void *params) {
  double *ar = (double *) params;

  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);

  const double ell = ar[2] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;

  const double ell_prefactor1 = (ar[2])*(ar[2] + 1.);
  double ell_prefactor2 = (ar[2] - 1.)*ell_prefactor1*(ar[2] + 2.);
  if(ell_prefactor2 <= 0.) {
    ell_prefactor2=0.;
  }
  else {
    ell_prefactor2=sqrt(ell_prefactor2);
  }

  const double chi_0 = f_K(ell/k);
  const double chi_1 = f_K((ell+1.)/k);
  if (chi_1 > chi(limits.a_min)){
    return 0;
  }
  const double a_0 = a_chi(chi_0);
  const double a_1 = a_chi(chi_1);

  double wgal = W_gal(a, ar[0], chidchi.chi, hoverh0);
  wgal += W_mag(a, fK, ar[0])*(ell_prefactor1/ell/ell -1.) ;
  const double res = (wgal + W_RSD(ell, a_0, a_1, ar[0]))*
    W_kappa(a, fK, ar[1])*chidchi.dchida/(fK*fK) * ell_prefactor2/(ell*ell);
  return res*p_lin(k,a);
}

double int_for_C_gl_IA_lin(double a, void *params) {
  double *ar = (double *) params;

  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);

  const double ell = ar[2]+0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;

  const double ell_prefactor1 = (ar[2])*(ar[2]+1.);
  double ell_prefactor2 = (ar[2]-1.)*ell_prefactor1*(ar[2]+2.);
  if(ell_prefactor2 <= 0.) {
    ell_prefactor2=0.;
  }
  else {
    ell_prefactor2=sqrt(ell_prefactor2);
  }

  const double chi_0 = f_K(ell/k);
  const double chi_1 = f_K((ell+1.)/k);
  if (chi_1 > chi(limits.a_min)){
    return 0;
  }
  const double a_0 = a_chi(chi_0);
  const double a_1 = a_chi(chi_1);

  const double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia*
    growfac(0.9999)/growfac(a)*nuisance.A_ia*
    pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);

  const double res = (W_gal(a, ar[0], hoverh0, chidchi.chi) +
    W_RSD(ell, a_0, a_1, ar[0]) + W_mag(a, fK, ar[0])*
    (ell_prefactor1/ell/ell -1.))*(W_kappa(a, fK, ar[1]) -
    W_source(a, ar[1], hoverh0)*norm, hoverh0)*chidchi.dchida/(fK*fK)*ell_prefactor2/(ell*ell);
  return res*p_lin(k,a);
}

void f_chi_for_Psi_sh(double* chi_ar, int Nchi, double* f_chi_ar, int ns) {
  double g0 =1./growfac(1.);
  double a, z, fK;
  int i;
  double real_coverH0 = cosmology.coverH0 / cosmology.h0;
  double window_L, wkappa;
  for(i=0;i<Nchi;i++) {
    a = a_chi(chi_ar[i] / real_coverH0) ; // first convert unit of chi from Mpc to c/H0
    z = 1./a - 1.;
    struct chis chidchi = chi_all(a);
    const double hoverh0 = hoverh0v2(a, chidchi.dchida);

    fK = f_K(chi_ar[i]/real_coverH0);
    if(z>tomo.shear_zmax[ns])
    {
      f_chi_ar[i] = 0.;
    }
    else
    {
      wkappa = W_kappa(a, fK, (double)ns);
      window_L = wkappa/ fK / (real_coverH0*real_coverH0);
      // printf("win_L, wkappa, f_K, real_coverH0, %lg %lg %lg %lg\n", window_L, wkappa, fK,real_coverH0);
      f_chi_ar[i] = window_L*growfac(a)*g0; // unit [Mpc^-2]
    }
  }
}

void f_chi_for_Psi_sh_IA(double* chi_ar, int Nchi, double* f_chi_IA_ar, int ns) {
  double g0 =1./growfac(1.);
  double a, z, fK;
  int i;
  double real_coverH0 = cosmology.coverH0 / cosmology.h0;
  double window_ia, wsource;
  double norm;
  for(i=0;i<Nchi;i++) {
    a = a_chi(chi_ar[i] / real_coverH0) ; // first convert unit of chi from Mpc to c/H0
    z = 1./a - 1.;
    struct chis chidchi = chi_all(a);
    const double hoverh0 = hoverh0v2(a, chidchi.dchida);

    fK = f_K(chi_ar[i]/real_coverH0);
    norm = cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(0.9999)/growfac(a)*nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);
    if( (z<tomo.shear_zmin[ns]) || (z>tomo.shear_zmax[ns]) )
    {
      f_chi_IA_ar[i] = 0.;
    }
    else
    {
      wsource = W_source(a, (double) ns, hoverh0);
      wsource = (wsource>0.)? wsource:0.;
      window_ia = -wsource * norm / fK / (real_coverH0*real_coverH0);
      f_chi_IA_ar[i] = window_ia*growfac(a)*g0; // unit [Mpc^-2]
    }
  }
}

double int_for_C_cl_lin(double a, void *params) {
  double *ar = (double*) params;

  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);

  const double ell = ar[2] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell / fK;

  const double res = W_gal(a, ar[0], chidchi.chi, hoverh0) *
    W_gal(a, ar[1], chidchi.chi, hoverh0) * chidchi.dchida / (fK * fK);

  return res * p_lin(k, a);
}

// test for replacing Plin with Pdelta and rescaled
double int_for_C_cl_nl_rescale(double a, void *params) {
  double *ar = (double*) params;

  double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  double hoverh0 = hoverh0v2(a, chidchi.dchida);

  const double ell = ar[2] + 0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell / fK;

  const double res = W_gal(a, ar[0], chidchi.chi, hoverh0) *
    W_gal(a, ar[1], chidchi.chi, hoverh0) * chidchi.dchida / (fK * fK);

  return res * Pdelta(k, 0.9999) * growfac_a * growfac_a;
}

// Integrand for galaxy density
void f_chi_for_Psi_cl(double *chi_ar, int Nchi, double *f_chi_ar, int ni) {
  const double real_coverH0 = cosmology.coverH0/cosmology.h0; // unit Mpc
  for (int i = 0; i < Nchi; i++) {
    const double a = a_chi(chi_ar[i]/real_coverH0 /* convert unit to c/H0 */);
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
    const double a = a_chi(chi_ar[i]/real_coverH0 /* convert unit to c/H0 */);
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
    const double a = a_chi(chi_ar[i]/real_coverH0 /* convert unit to c/H0 */);
    const double z = 1. / a - 1.;
    const double fK = f_K(chi_ar[i]/real_coverH0 /* convert unit to c/H0 */);
    if ((z > tomo.clustering_zmax[ni])) {
      f_chi_Mag_ar[i] = 0.;
    } else {
      const double wmag = W_mag(a, fK, (double) ni);
      const double window_M = wmag / fK / (real_coverH0 * real_coverH0);
      f_chi_Mag_ar[i] = window_M * growfac(a); // unit [Mpc^-2]
    }
  }
}

// ----------------------------------------------------------------------------
// Angular Power Spectra
// ----------------------------------------------------------------------------

//galaxy clustering power spectrum of galaxy bins ni, nj
double C_gl_lin_nointerp(double l, int ni, int nj)
{
  double array[3] = {1.0*ni,1.0*nj,l};
  return int_gsl_integrate_medium_precision(int_for_C_gl_lin,(void*)array,amin_lens(ni),amax_lens(ni),NULL,1000);
}

//galaxy clustering power spectrum of galaxy bins ni, nj
double C_gl_lin_IA_nointerp(double l, int ni, int nj)
{
  double array[3] = {1.0*ni,1.0*nj,l};
  return int_gsl_integrate_medium_precision(int_for_C_gl_IA_lin,(void*)array,amin_lens(ni),0.9999,NULL,1000);
}

// galaxy clustering power spectrum of galaxy bins ni, nj
double C_cl_lin_nointerp(double l, int ni, int nj) {
  double array[3] = {1.0 * ni, 1.0 * nj, l};
  return int_gsl_integrate_low_precision(int_for_C_cl_lin, (void *)array,
                fmax(amin_lens(ni), amin_lens(nj)), 0.99999, NULL, 1000);
}

// galaxy clustering power spectrum of galaxy bins ni, nj
double C_cl_nl_rescaled_nointerp(double l, int ni, int nj) {
  double array[3] = {1.0 * ni, 1.0 * nj, l};
  return int_gsl_integrate_low_precision(int_for_C_cl_nl_rescale,
    (void *)array, fmax(amin_lens(ni), amin_lens(nj)), 0.99999, NULL, 1000);
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
  const double chi_min = chi(1./(1.+0.002))*real_coverH0; // DIMENSIONELESS
  const double chi_max = chi(1./(1.+4.))*real_coverH0; // DIMENSIONELESS
  const double dlnchi = log(chi_max/chi_min) / (Nchi - 1.);
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
      Fk1_ar[i][j] = 0.0;
      Fk2_ar[i][j] = 0.0;
      Fk1_Mag_ar[i][j] = 0.0;
      Fk2_Mag_ar[i][j] = 0.0;
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

// Mixture of non-Limber and Limber of C_cl (G-G lensing)
void C_gl_mixed(int L, int LMAX, int nl, int ns, double *Cl, double dev, double tolerance) {
  // nl = 4;
  // ns = 0;
  int i,j,i_block;
  long l;
  // run 100 ells at a time, and see if switching to Limber is needed.
  // Save runtime for Limber, and save re-creation time of fftw_plan.
  int Nell_block = 100, Nchi = 1000;
  int ell_ar[Nell_block];
  double **k1_ar, **k2_ar, **Fk1_ar, **Fk2_ar;
  double **Fk1_Mag_ar;

  k1_ar = malloc(Nell_block * sizeof(double *));
  k2_ar = malloc(Nell_block * sizeof(double *));
  Fk1_ar = malloc(Nell_block * sizeof(double *));
  Fk2_ar = malloc(Nell_block * sizeof(double *));

  Fk1_Mag_ar = malloc(Nell_block * sizeof(double *));
  for(i=0;i<Nell_block;i++) {
    k1_ar[i] = malloc(Nchi * sizeof(double));
    k2_ar[i] = malloc(Nchi * sizeof(double));
    Fk1_ar[i] = malloc(Nchi * sizeof(double));
    Fk2_ar[i] = malloc(Nchi * sizeof(double));
    Fk1_Mag_ar[i] = malloc(Nchi * sizeof(double));
    for(j=0;j<Nchi;j++) {
      Fk1_ar[i][j] = 0.;
      Fk2_ar[i][j] = 0.;
      Fk1_Mag_ar[i][j] = 0.;
    }
  }

  double chi_ar[Nchi];
  double f1_chi_ar[Nchi], f1_chi_RSD_ar[Nchi], f1_chi_Mag_ar[Nchi];
  double f2_chi_ar[Nchi], f2_chi_IA_ar[Nchi];

  // double f2_chi_temp[Nchi];

  double chi_min = 10., chi_max = 7000.;
  // double chi_min = 6., chi_max = 6000.;

  double dlnchi = log(chi_max/chi_min) / (Nchi - 1.);
  double dlnk = dlnchi;

  for(i=0; i<Nchi; i++) {
    chi_ar[i] = chi_min * exp(dlnchi*i);
  }
  f_chi_for_Psi_cl(chi_ar, Nchi, f1_chi_ar, nl);
  f_chi_for_Psi_cl_RSD(chi_ar, Nchi, f1_chi_RSD_ar, nl);
  f_chi_for_Psi_cl_Mag(chi_ar, Nchi, f1_chi_Mag_ar, nl);
  f_chi_for_Psi_sh(chi_ar, Nchi, f2_chi_ar, ns);
  f_chi_for_Psi_sh_IA(chi_ar, Nchi, f2_chi_IA_ar, ns);
  for(j=0;j<Nchi;j++) {
    f2_chi_ar[j] += f2_chi_IA_ar[j];
  }

  i_block = 0;
  double cl_temp;

  config my_config, my_config_RSD, my_config_Mag, my_config_L;
  my_config.nu = 1.;
  my_config.c_window_width = 0.25;
  my_config.derivative = 0;
  my_config.N_pad = 200;
  my_config.N_extrap_low = 0;
  my_config.N_extrap_high = 0;

  my_config_RSD.nu = 1.01;
  my_config_RSD.c_window_width = 0.25;
  my_config_RSD.derivative = 2;
  my_config_RSD.N_pad = 200;
  my_config_RSD.N_extrap_low = 0;
  my_config_RSD.N_extrap_high = 0;

  my_config_Mag.nu = 1.;
  my_config_Mag.c_window_width = 0.25;
  my_config_Mag.derivative = 0;
  my_config_Mag.N_pad = 1000;
  my_config_Mag.N_extrap_low = 0;
  my_config_Mag.N_extrap_high = 0;

  my_config_L.nu = 1.;
  my_config_L.c_window_width = 0.25;
  my_config_L.derivative = 0;
  my_config_L.N_pad = 1000.;
  my_config_L.N_extrap_low = 0;
  my_config_L.N_extrap_high = 0;

  double ell_prefactor, ell_prefactor2;

  double real_coverH0 = cosmology.coverH0 / cosmology.h0;
  double k1_cH0;

  while (fabs(dev) > tolerance){
    for(i=0;i<Nell_block;i++) {ell_ar[i]=i+i_block*Nell_block;}
    cfftlog_ells(chi_ar, f1_chi_ar, Nchi, &my_config, ell_ar, Nell_block, k1_ar, Fk1_ar);
    cfftlog_ells_increment(chi_ar, f1_chi_RSD_ar, Nchi, &my_config_RSD, ell_ar, Nell_block, k1_ar, Fk1_ar);

    // Add in lensing magnification contribution
    cfftlog_ells(chi_ar, f1_chi_Mag_ar, Nchi, &my_config_Mag, ell_ar, Nell_block, k1_ar, Fk1_Mag_ar);
    for(i=0;i<Nell_block;i++) {
      ell_prefactor = ell_ar[i]*(ell_ar[i]+1);
      for(j=0;j<Nchi;j++) {
        Fk1_ar[i][j]+= (ell_prefactor / (k1_ar[i][j]*k1_ar[i][j])* (gbias.b_mag[nl]) * Fk1_Mag_ar[i][j]);
      }
    }

    // shear part
    cfftlog_ells(chi_ar, f2_chi_ar, Nchi, &my_config_L, ell_ar, Nell_block, k2_ar, Fk2_ar);

    for(i=0;i<Nell_block;i++) {
      ell_prefactor2 =(ell_ar[i]-1.)*ell_ar[i]*(ell_ar[i]+1.)*(ell_ar[i]+2.);
      if(ell_prefactor2<=0.) {ell_prefactor2=0.;}
      else {ell_prefactor2 = sqrt(ell_prefactor2);}
      for(j=0;j<Nchi;j++) {
        Fk2_ar[i][j]*= (ell_prefactor2 / (k1_ar[i][j]*k1_ar[i][j]));
      }
    }

    for(i=0;i<Nell_block;i++) {
      cl_temp = 0.;
      for(j=0;j<Nchi;j++) {
        k1_cH0 = k1_ar[i][j] * real_coverH0;
        cl_temp += (Fk1_ar[i][j])*(Fk2_ar[i][j]) *k1_cH0*k1_cH0*k1_cH0 *p_lin(k1_cH0,1.0);
      }
      Cl[ell_ar[i]] = cl_temp * dlnk * 2./M_PI + C_ggl_IA(1.*ell_ar[i],nl,ns) - C_gl_lin_IA_nointerp(1.*ell_ar[i],nl,ns);
      dev = Cl[ell_ar[i]]/C_ggl_IA(1.0*ell_ar[i],nl,ns)-1.;
    }

    i_block++;
    L = i_block*Nell_block -1 ;
    dev = Cl[L]/C_ggl_IA(1.0*L,nl,ns)-1.;
  }
  L++;


  for (l = L; l < LMAX; l++){
    // Cl[l]=C_gl_tomo((double)l,nl,ns);
    Cl[l]=C_ggl_IA_tab((double)l,nl,ns);
  }
  // printf("finished bin %d %d\n", nl,ns);
  for(i=0;i<Nell_block;i++) {
    free(k1_ar[i]);free(k2_ar[i]);
    free(Fk1_ar[i]);free(Fk2_ar[i]);
    free(Fk1_Mag_ar[i]);
  }
  free(k1_ar);free(k2_ar);
  free(Fk1_ar);free(Fk2_ar);
  free(Fk1_Mag_ar);
  // fclose(OUT);
  // exit(0);
}

// ---------------------------------------------------------------------------
// Angular Functions (real Space)
// ---------------------------------------------------------------------------

double w_gamma_t_fullsky_nonLimber(int nt, int ni, int nj) {
  // if(1) return 0.;
  static int LMAX = 100000;
  static int NTHETA = 0;
  static double ** Pl =0;
  static double *Cl =0;
  static double *w_vec =0;
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  int i,l,nz;
  if (like.Ntheta ==0){
    printf("cosmo2D_fullsky.c:w_gamma_t_tomo: like.Ntheta not initialized\nEXIT\n"); exit(1);
  }
  if (Pl ==0){
    Pl =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
    Cl = create_double_vector(0,LMAX-1);
    w_vec = create_double_vector(0,tomo.ggl_Npowerspectra*like.Ntheta-1);
    NTHETA = like.Ntheta;
    double *xmin, *xmax, *Pmin, *Pmax, *dP;
    xmin= create_double_vector(0, like.Ntheta-1);
    xmax= create_double_vector(0, like.Ntheta-1);
    double logdt=(log(like.vtmax)-log(like.vtmin))/like.Ntheta;

    // double *xmid, *Pmid;
    double mythetamin, mythetamax;
    // xmid= create_double_vector(0, like.Ntheta-1);
    // Pmid= create_double_vector(0, LMAX+1);
    for(i=0; i<like.Ntheta ; i++){
      mythetamin = exp(log(like.vtmin)+(i+0.0)*logdt);
      mythetamax = exp(log(like.vtmin)+(i+1.0)*logdt);
      xmin[i]=cos(mythetamin);
      xmax[i]=cos(mythetamax);
      // xmid[i]= (2./3.) * (pow(thetamax,3) - pow(thetamin,3)) / (thetamax*thetamax - thetamin*thetamin);
    }
    Pmin= create_double_vector(0, LMAX+1);
      Pmax= create_double_vector(0, LMAX+1);

    for (i = 0; i<NTHETA; i ++){
      // printf("Tabulating Legendre coefficients %d/%d\n",i+1, NTHETA);
      gsl_sf_legendre_Pl_array(LMAX, xmin[i],Pmin);
      gsl_sf_legendre_Pl_array(LMAX, xmax[i],Pmax);
        // gsl_sf_legendre_Pl_array(LMAX, xmid[i],Pmid);
      for (int l = 2; l < LMAX; l ++){
        // Pl[i][l] = (2.*l+1)/(4.*M_PI*l*(l+1))*gsl_sf_legendre_Plm(l,2,cos(like.theta[i]));
        Pl[i][l] = (2.*l+1)/(4.*M_PI*l*(l+1)*(xmin[i]-xmax[i]))
        *((l+2./(2*l+1.))*(Pmin[l-1]-Pmax[l-1])
        +(2-l)*(xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        -2./(2*l+1.)*(Pmin[l+1]-Pmax[l]));
        // Pl[i][l] = (2.*l+1)/(4.*M_PI*l*(l+1))*Pmid[l]*Pmid[l];
      }
    }
    free_double_vector(xmin,0,like.Ntheta-1);
    free_double_vector(xmax,0,like.Ntheta-1);
    free_double_vector(Pmin,0,LMAX+1);
    free_double_vector(Pmax,0,LMAX+1);
  }
  if (recompute_ggl(C,G,N,ni)){
    //required fractional accuracy in C(l)
    double tolerance= 0.01;
    //dev will be the actual difference between exact and Limber calcuation
    double dev;

    // if (like.IA != 0 && like.IA != 3 && like.IA != 4){printf("cosmo2D_real.c: w_gamma_t_tomo does not support like.IA = %d yet\nEXIT!\n",like.IA); exit(1);}
    // C_tomo_pointer C_gl_pointer = &C_gl_tomo;
    // if (like.IA ==3 || like.IA ==4) C_gl_pointer = &C_ggl_IA_tab;

    for (nz = 0; nz <tomo.ggl_Npowerspectra; nz ++){
      int L = 1;
      // initialize to large value in order to start while loop
      dev=10.*tolerance;
      C_gl_mixed(L, LMAX, ZL(nz),ZS(nz), Cl, dev, tolerance);
      for (i = 0; i < NTHETA; i++){
        w_vec[nz*like.Ntheta+i] =0;
        for (l = 1; l < LMAX; l++){
          w_vec[nz*like.Ntheta+i]+=Pl[i][l]*Cl[l];
        }
      }
    }
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  return w_vec[N_ggl(ni,nj)*like.Ntheta+nt];
}

double w_tomo_fullsky_nonLimber(int nt, int ni, int nj) {
  static int LMAX = 100000;
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
    Pl = create_double_matrix(0, like.Ntheta - 1, 0, LMAX - 1);
    w_vec = create_double_vector(0, tomo.clustering_Nbin * like.Ntheta - 1);
    double xmin[like.Ntheta];
    double xmax[like.Ntheta];
    const double logdt = (log(like.vtmax) - log(like.vtmin)) / like.Ntheta;

    #pragma omp parallel for
    for (int i = 0; i < like.Ntheta; i++) {
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

      free_double_vector(Pmin, 0, LMAX + 1);
      free_double_vector(Pmax, 0, LMAX + 1);
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
      for (int i = 0; i < like.Ntheta; i++) {
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
