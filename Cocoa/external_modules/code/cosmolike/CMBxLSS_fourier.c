#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "basics.h"
#include "bias.h"
#include "cosmo3D.h"
#include "cosmo2D_fourier.h"
#include "CMBxLSS_fourier.h"
#include "pt_cfastpt.h"
#include "recompute.h"
#include "radial_weights.h"
#include "redshift_spline.h"
#include "structs.h"

// Naming convention:
// g = galaxy positions ("g" as in "galaxy")
// k = kappa CMB ("k" as in "kappa")
// s = kappa from source galaxies ("s" as in "shear")
// And alphabetical order

//lensing efficiency weightfunction for CMB lensing
double W_k(double a, double fK){
  return 1.5*cosmology.Omega_m*fK/a*g_cmb(a);
}

// galaxy position x kappa CMB
double int_for_C_gk(double a, void *params) {
  double *ar = (double *) params;
  const double growfac_a = growfac(a);
  struct chis chidchi = chi_all(a);
  const double hoverh0 = hoverh0v2(a, chidchi.dchida);
  const double ell = ar[1]+0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  //W_gal radial weight function for clustering, defined in cosmo2D_fourier.c
  const double res = W_gal(a, ar[0], chidchi.chi, hoverh0)/(fK*fK);
  return res*W_k(a,fK)*chidchi.dchida*Pdelta(k,a);
}

double int_for_C_gk_b2(double a, void *params) {
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

// shear x kappa CMB
double int_for_C_ks(double a, void *params) {
  double *ar = (double *) params;
  if (a >= 1.0) {
    error("a >1 in int_for_C_ks");
  }
  struct chis chidchi = chi_all(a);
  const double ell = ar[1]+0.5;
  const double fK = f_K(chidchi.chi);
  const double k = ell/fK;
  const double res = W_kappa(a, fK, ar[0])/(fK*fK);
  return res*W_k(a, fK)*chidchi.dchida*Pdelta(k,a);
}

// kappa CMB x kappaCMB
double int_for_C_kk(double a, void *params) {
   double *ar = (double *) params;
   if (a > 1.0) {
    error("a >1 in int_for_C_kk");
   }
   struct chis chidchi = chi_all(a);
   const double ell = ar[0]+0.5;
   const double fK = f_K(chidchi.chi);
   const double k = ell/fK;
   const double res = pow(W_k(a,fK), 2)*chidchi.dchida/(fK*fK);
   return res*Pdelta(k,a);
}

double int_for_C_ks_IA_mpp(double a, void *params) { // for like.IA==4
   double *ar = (double *) params;
   struct chis chidchi = chi_all(a);
   double hoverh0 = hoverh0v2(a, chidchi.dchida);
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

double C_ks_IA(double s, int ni) {
  double array[2] = {(double) ni,s};
  if (like.IA==4) {
    return int_gsl_integrate_medium_precision(int_for_C_ks_IA_mpp,
      (void*)array, amin_source(ni), amax_source(ni), NULL, 1000);
  }
  printf("CMBxLSS.c: C_ks_IA does not support like.IA = %d\nEXIT\n", like.IA);
  exit(1);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// power spectra - no look-up tables
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// galaxy position x kappa CMB, lens z-bin nl
double C_gk_nointerp(double l, int nl) {
  double array[2] = {(double)nl,l};
  if (gbias.b2[nl] || gbias.b2[nl]) {
    return int_gsl_integrate_medium_precision(int_for_C_gk_b2,
      (void*)array, amin_lens(nl), amax_lens(nl), NULL, 1000);
  }
  return int_gsl_integrate_medium_precision(int_for_C_gk, (void*)array,
    amin_lens(nl), 0.99999, NULL, 1000);
}

// shear x kappa CMB, for source z-bin ns
double C_ks_nointerp(double l, int ns) {
  if (like.IA) {
    return C_ks_IA(l,ns);
  }
  double array[2] = {(double) ns, l};
  return int_gsl_integrate_medium_precision(int_for_C_ks, (void*)array,
    amin_source(ns), amax_source(ns), NULL,1000);
}

// kappa CMB x kappa CMB
double C_kk_nointerp(double l) {
  double array[1] = {l};
  return int_gsl_integrate_medium_precision(int_for_C_kk, (void*)array,
    limits.a_min*(1.+1.e-5), 1.-1.e-5, NULL, 1000);
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// power spectra - look-up tables
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// galaxy position x kappa CMB, lens bin ni
double C_gk(double l, int ni) {
   static cosmopara C;
   static nuisancepara N;
   static galpara G;
   static double **table;
   static double ds = .0, logsmin = .0, logsmax = .0;
   if (ni < 0 || ni >= tomo.clustering_Nbin){
      printf("Bin %d outside tomo.clustering_Nbin range\nEXIT\n",ni); exit(1);
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
      table[k][0]= log(C_gk_nointerp(exp(logsmin),k));
      #pragma omp parallel for
      for (int i=1; i<Ntable.N_ell; i++) {
        const double llog = logsmin + i*ds;
        table[k][i]= log(C_gk_nointerp(exp(llog),k));
      }
    }
    #pragma omp parallel for
    for (int k=1; k<tomo.clustering_Nbin; k++) {
      for (int i=0; i<Ntable.N_ell; i++) {
        const double llog = logsmin + i*ds;
        table[k][i]= log(C_gk_nointerp(exp(llog),k));
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

// shear x kappa CMB, source bin ni
double C_ks(double l, int ni) {
   static cosmopara C;
   static nuisancepara N;
   static double **table, *sig;
   static int osc[100];
   static double ds = .0, logsmin = .0, logsmax = .0;
   if (ni < 0 || ni >= tomo.clustering_Nbin){
      printf("Bin %d outside tomo.clustering_Nbin range\nEXIT\n",ni); exit(1);
   }
   if (recompute_shear(C,N))
   {
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
      const double res = C_ks_nointerp(500.,k);
      if (res < 0){
        sig[k] = -1.;
      }
      #pragma omp parallel for
      for (int i=0; i<Ntable.N_ell; i++) {
        const double llog = logsmin + i*ds;
        table[k][i]= C_ks_nointerp(exp(llog), k);
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
      const double res = C_ks_nointerp(500.,k);
      if (res < 0){
        sig[k] = -1.;
      }
      for (int i=0; i<Ntable.N_ell; i++) {
        const double llog = logsmin + i*ds;
        table[k][i]= C_ks_nointerp(exp(llog), k);
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

// kappa CMB x kappa CMB
double C_kk(double l) {
  static cosmopara C;
  static double *table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  if (recompute_cosmo3D(C))
  {
    if (table == 0) {
       table   = create_double_vector(0, Ntable.N_ell-1);
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
       table[i]= log(C_kk_nointerp(exp(llog)));
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