#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cfastpt/cfastpt.h"
#include "basics.h"
#include "cosmo3D.h"
#include "pt_cfastpt.h"
#include "recompute.h"
#include "structs.h"

#include "log.c/src/log.h"

void FPT_input(double k[FPT.N], double P[FPT.N]) 
{
  if ((int) log10(FPT.k_max / FPT.k_min) * FPT.N_per_dec != FPT.N) 
  {
    log_fatal("inconsistent k-range and number of bins for FPT");
    log_fatal("FPT.k_min=%e, FPT.k_max=%e, FPT.N_per_dec=%d; FPT.N=%d", FPT.k_min, FPT.k_max, 
      FPT.N_per_dec, FPT.N);
    exit(1);
  }
  if (FPT.k_min < limits.k_min_cH0 || FPT.k_max > limits.k_max_cH0) 
  {
    log_fatal("k_min/k_max out of range FPT.k_min = %e, FPT.k_max =%e, "
      "limits.k_min = %e, and limits.k_max = %e", FPT.k_min, FPT.k_max,
      limits.k_min_cH0, limits.k_max_cH0);
    exit(1);
  }

  const double dlgk = log(10.) / (double)FPT.N_per_dec;
  const double tmp = log(FPT.k_min);
  {
    const int i = 0;
    k[i] = exp(tmp  + i*dlgk);
    P[i] = p_lin(k[i], 1.0);
  }
  #pragma omp parallel for
  for (int i=1; i<FPT.N; i++) 
  {
    k[i] = exp(tmp  + i*dlgk);
    P[i] = p_lin(k[i], 1.0);
  }
}

void get_FPT_bias(void) 
{
  static double cache[MAX_SIZE_ARRAYS];

  if (fdiff(cache[1], Ntable.random))
  {
    if ((int) log10(FPT.k_max / FPT.k_min) * FPT.N_per_dec != FPT.N) 
    {
      log_fatal("inconsistent k-range and number of bins for FPT");
      log_fatal("inconsistent k-range and number of bins for FPT: "
          "kmin=%e, kmax=%e, N_per_dec=%d, N=%d", FPT.k_min, 
          FPT.k_max, FPT.N_per_dec, FPT.N);
      exit(1);
    }

    if (FPT.tab_AB != NULL) free(FPT.tab_AB);
    FPT.tab_AB = (double**) malloc2d(7, FPT.N);
    FPT.k_min = 1.e-5*cosmology.coverH0;
    FPT.k_max = 1.e+3*cosmology.coverH0;
  }
  if (fdiff(cache[0], cosmology.random) || fdiff(cache[1], Ntable.random))
  {
    #pragma omp parallel for
    for (int i=0; i<FPT.N; i++) 
    {
      FPT.tab_AB[5][i] = exp(log(FPT.k_min) + i*(log(10.)/FPT.N_per_dec));
      FPT.tab_AB[6][i] = p_lin(FPT.tab_AB[5][i], 1.0);
    }

    double Pout[5][FPT.N];
    Pd1d2(FPT.tab_AB[5], FPT.tab_AB[6], FPT.N, Pout[0]);
    Pd2d2(FPT.tab_AB[5], FPT.tab_AB[6], FPT.N, Pout[1]);
    Pd1s2(FPT.tab_AB[5], FPT.tab_AB[6], FPT.N, Pout[2]);
    Pd2s2(FPT.tab_AB[5], FPT.tab_AB[6], FPT.N, Pout[3]);
    Ps2s2(FPT.tab_AB[5], FPT.tab_AB[6], FPT.N, Pout[4]);

    #pragma omp parallel for
    for (int i=0; i<FPT.N; i++) 
    {
      FPT.tab_AB[0][i] = Pout[0][i]; // Pd1d2
      FPT.tab_AB[1][i] = Pout[1][i]; // Pd2d2
      FPT.tab_AB[2][i] = Pout[2][i]; // Pd1s2
      FPT.tab_AB[3][i] = Pout[3][i]; // Pd2s2
      FPT.tab_AB[4][i] = Pout[3][i]; // Pd2s2
    }

    cache[0] = cosmology.random;
    cache[1] = Ntable.random;
  }
}

void get_FPT_IA(void) 
{
  static double cache[MAX_SIZE_ARRAYS];

  if (fdiff(cache[1], Ntable.random))
  {
    if ((int) log10(FPT.k_max / FPT.k_min) * FPT.N_per_dec != FPT.N) 
    {
      log_fatal("inconsistent k-range and number of bins for FPT");
      log_fatal("inconsistent k-range and number of bins for FPT: "
          "kmin=%e, kmax=%e, N_per_dec=%d, N=%d", FPT.k_min, 
          FPT.k_max, FPT.N_per_dec, FPT.N);
      exit(1);
    }

    if (FPT.tab_IA != NULL) free(FPT.tab_IA);
    FPT.tab_IA = (double**) malloc2d(12, FPT.N);
    FPT.k_min = 1.e-5*cosmology.coverH0;
    FPT.k_max = 1.e+3*cosmology.coverH0;
  }
  if (fdiff(cache[0], cosmology.random) || fdiff(cache[1], Ntable.random))
  {
    #pragma omp parallel for
    for (int i=0; i<FPT.N; i++) 
    {
      FPT.tab_IA[10][i] = exp(log(FPT.k_min) + i*(log(10.)/FPT.N_per_dec));
      FPT.tab_IA[11][i] = p_lin(FPT.tab_IA[10][i], 1.0);
    }

    IA_tt(FPT.tab_IA[10], FPT.tab_IA[11], FPT.N, FPT.tab_IA[0], FPT.tab_IA[1]);
    
    IA_ta(FPT.tab_IA[10], FPT.tab_IA[11], FPT.N, FPT.tab_IA[2], 
      FPT.tab_IA[3], FPT.tab_IA[4], FPT.tab_IA[5]);
    
    IA_mix(FPT.tab_IA[10], FPT.tab_IA[11], FPT.N,FPT.tab_IA[6], 
      FPT.tab_IA[7], FPT.tab_IA[8], FPT.tab_IA[9]);
    
    #pragma omp parallel for
    for (int i=0; i<FPT.N; i++) 
      FPT.tab_IA[7][i] *= 4.;
    
    cache[0] = cosmology.random;
    cache[1] = Ntable.random;
  }
}
