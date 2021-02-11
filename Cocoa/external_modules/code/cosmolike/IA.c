#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "basics.h"
#include "cosmo3D.h"
#include "cosmo2D_fourier.h"
#include "pt_cfastpt.h"
#include "radial_weights.h"
#include "recompute.h"
#include "redshift_spline.h"
#include "structs.h"
#include "IA.h"
#include "bias.h"

#include "log.c/src/log.h"

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// faster shear/ggl + IA routines
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

double int_for_C_ggl_IA_mpp(double a, void *params) {
  if (a >= 1.0) {
    log_fatal("a>=1");
    exit(1);
  }
  static double chi_a_min = 0;
  if (chi_a_min == 0) {
    chi_a_min = chi(limits.a_min);
  }
  double *ar = (double *) params;
  const double ell_prefactor1 = (ar[2])*(ar[2]+1.);
  double ell_prefactor2 = (ar[2]-1.)*ell_prefactor1*(ar[2]+2.);
  if (ell_prefactor2 <= 0.) {
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
  const double chi_0 = f_K(ell/k);
  const double chi_1 = f_K((ell+1.)/k);
  if (chi_1 > chi_a_min) {
    return 0;
  }
  const double a_0 = a_chi(chi_0);
  const double a_1 = a_chi(chi_1);
  const double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac(a)*
    nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);
  const double res = (W_gal(a, ar[0], chidchi.chi, hoverh0) +
    W_RSD(ell, a_0, a_1, ar[0]) +
    W_mag(a, fK, ar[0])*(ell_prefactor1/ell/ell - 1.))*
    (W_kappa(a, fK, ar[1]) - W_source(a, ar[1], hoverh0)*norm);
  return res*Pdelta(k,a)*chidchi.dchida/(fK*fK)*ell_prefactor2/(ell*ell);
}

double int_for_C_ggl_IA_mpp_b2(double a, void *params) {
  if (a >= 1.0) {
    log_fatal("a>=1");
    exit(1);
  }
  double *ar = (double *) params;
  const double b1 = gbias.b1_function(1./a-1.,(int)ar[0]);
  const double b2 = gbias.b2[(int)ar[0]];
  const double bs2 = gbias.bs2[(int)ar[0]];
  const double grow_fac = growfac(a);
  const double g4 = grow_fac*grow_fac*grow_fac*grow_fac;
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double ell = ar[2]+0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double norm =
    (cosmology.Omega_m*nuisance.c1rhocrit_ia/grow_fac)*
    nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia), nuisance.eta_ia);
  const double res = W_HOD(a, ar[0], hoverh0)*
    (W_kappa(a, fK, ar[1]) - W_source(a, ar[1], hoverh0)*norm);
  return res*(
      W_gal(a, ar[0], chidchi.chi,hoverh0)*Pdelta(k, a) +
      W_HOD(a, ar[0], hoverh0)*g4*(
        0.5*b2*PT_d1d2(k) +
        0.5*bs2*PT_d1s2(k) +
        0.5*b3nl_from_b1(b1)*PT_d1d3(k)
      )
    )*chidchi.dchida/(fK*fK);
}

double int_for_C_shear_shear_IA_mpp(double a, void *params) {
  if (a >= 1.0) {
    log_fatal("a>=1");
    exit(1);
  }
  double *ar = (double *) params;
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double ell = ar[2]+0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double ws1 = W_source(a, ar[0], hoverh0);
  const double ws2 = W_source(a, ar[1], hoverh0);
  const double wk1 = W_kappa(a, fK, ar[0]);
  const double wk2 = W_kappa(a, fK, ar[1]);
  const double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia/growfac(a)*
    nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);
  const double res = ws1*ws2*norm*norm - (ws1*wk2+ws2*wk1)*norm + wk1*wk2;
  return res*Pdelta(k,a)*chidchi.dchida/(fK*fK);
}

// COBAYA: - WE CHANGE SOME INTEGRATIONS FROM MEDIUM TO LOW PRECISION!!
double C_shear_shear_IA(double s, int ni, int nj) {
  double array[3] = {(double) ni, (double) nj,s};
  int j,k;
  if (ni <= nj) {
    j = nj;
    k = ni;
  }
  switch(like.IA){
    case 4:
      return int_gsl_integrate_low_precision(int_for_C_shear_shear_IA_mpp,
        (void*) array,amin_source(j), amax_source(k), NULL, 1000);
    default:
      log_fatal("like.IA = %d not supported", like.IA);
      exit(1);
  }
}

// VM - COBAYA - WE CHANGE SOME INTEGRATIONS FROM MEDIUM TO LOW PRECISION!!
double C_ggl_IA(double s, int nl, int ns) {
  double array[3] = {(double) nl, (double) ns,s};
  switch(like.IA){
    case 4:
      if (gbias.b2[nl]) {
        return int_gsl_integrate_medium_precision(
          int_for_C_ggl_IA_mpp_b2,
          (void*) array,
          amin_lens(nl),
          amax_lens(nl),
          NULL,
          1000
        );
      }
      return int_gsl_integrate_low_precision(
        int_for_C_ggl_IA_mpp,
        (void*) array,
        amin_lens(nl),
        0.9999,
        NULL,
        1000
      );
    default:
      log_fatal("like.IA = %d not supported", like.IA);
      exit(1);
  }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// tabulated versions of C_ggl_IA and C_shear_shear_IA
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// G-G lensing power spectrum, lens bin ni, source bin nj
double C_ggl_IA_tab(double l, int ni, int nj) {
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double **table, *sig;
  static int osc[100];
  static double ds = .0, logsmin = .0, logsmax = .0;
  if ( ni < 0 || ni >= tomo.clustering_Nbin ||
       nj < 0 || nj >= tomo.shear_Nbin ) {
    log_fatal("C_gl_tomo(l,%d,%d) outside tomo.X_Nbin range", ni, nj);
    exit(1);
  }
  if (recompute_ggl(C,G,N,ni)){
    if (table == 0) {
      table = create_double_matrix(0, tomo.ggl_Npowerspectra - 1, 0,
        Ntable.N_ell - 1);
      sig = create_double_vector(0,
        tomo.ggl_Npowerspectra - 1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell);
    }
    {
      const int k = 0;
      sig[k] = 1.;
      osc[k] = 0;
      const double res = C_ggl_IA(500., ZL(k), ZS(k));
      if (res < 0) {
        sig[k] = -1.;
      }
      #pragma omp parallel for
      for (int i=0; i<Ntable.N_ell; i++) {
        const double llog = logsmin + i*ds;
        table[k][i] = C_ggl_IA(exp(llog), ZL(k), ZS(k));
        if (res*sig[k] < 0.) {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0) {
        for(int i = 0; i < Ntable.N_ell; i++) {
          table[k][i] = log(sig[k]*table[k][i]);
        }
      }
    }
    #pragma omp parallel for
    for (int k = 1; k < tomo.ggl_Npowerspectra; k++) {
      sig[k] = 1.;
      osc[k] = 0;
      const double res = C_ggl_IA(500., ZL(k), ZS(k));
      if (res < 0) {
        sig[k] = -1.;
      }
      for (int i = 0; i < Ntable.N_ell; i++) {
        const double llog = logsmin + i*ds;
        table[k][i] = C_ggl_IA(exp(llog), ZL(k), ZS(k));
        if (res*sig[k] <0.) {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0) {
        for (int i = 0; i < Ntable.N_ell; i++) {
          table[k][i] = log(sig[k]*table[k][i]);
        }
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
    update_galpara(&G);
  }
  const int k = N_ggl(ni, nj);
  const int tmp = test_zoverlap(ni, nj);
  double f1 = 0.;
  if (tmp && osc[k] == 0) {
    f1 = sig[k]*exp(interpol_fitslope(table[k], Ntable.N_ell, logsmin,
      logsmax, ds, log(l), 1.));
  }
  else if (tmp && osc[k] ==1) {
    f1 = interpol_fitslope(table[k], Ntable.N_ell, logsmin, logsmax, ds,
      log(l), 1.);
  }
  if (isnan(f1)) {
    f1 = 0;
  }
  return f1;
}

//shear power spectrum of source galaxies in bins ni, nj
double C_shear_shear_IA_tab(double l, int ni, int nj) {
  static cosmopara C;
  static nuisancepara N;
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  if (ni < 0 || ni >= tomo.shear_Nbin || nj < 0 || nj >= tomo.shear_Nbin) {
    log_fatal("C_shear_tomo(l,%d,%d) outside tomo.shear_Nbin range", ni, nj);
    exit(1);
  }
  if (recompute_shear(C,N)) {
    if (table==0) {
      table = create_double_matrix(0, tomo.shear_Npowerspectra - 1, 0,
        Ntable.N_ell - 1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell);
    }
    {
      const int k = 0;
      table[k][0]= log(C_shear_shear_IA(exp(logsmin), Z1(k), Z2(k)));
      #pragma omp parallel for
      for (int i=1; i<Ntable.N_ell; i++) {
        const double llog = logsmin + i*ds;
        table[k][i] = log(C_shear_shear_IA(exp(llog), Z1(k), Z2(k)));
      }
    }
    #pragma omp parallel for
    for (int k=1; k<tomo.shear_Npowerspectra; k++) {
      for (int i=0; i<Ntable.N_ell; i++) {
        const double llog = logsmin + i*ds;
        table[k][i] = log(C_shear_shear_IA(exp(llog), Z1(k), Z2(k)));
      }
    }
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  double f1 = exp(interpol_fitslope(table[N_shear(ni,nj)], Ntable.N_ell,
    logsmin, logsmax, ds, log(l), 1.));
  if (isnan(f1)) {
    f1 = 0.;
  }
  return f1;
}