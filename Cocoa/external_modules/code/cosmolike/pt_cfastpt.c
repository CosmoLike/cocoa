#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cfastpt/cfastpt.h"
#include "basics.h"
#include "cosmo3D.h"
#include "pt_cfastpt.h"
#include "structs.h"

#include "log.c/src/log.h"

void get_FPT_bias(void) 
{
  static double cache[MAX_SIZE_ARRAYS];

  if (fdiff(cache[1], Ntable.random))
  {
    FPTIA.k_min     = 1.e-5;
    FPTIA.k_max     = 1.e+5;
    FPTbias.N       = 800 + 200 * Ntable.FPTboost;
    if (FPTbias.tab != NULL) {
      free(FPTbias.tab);
    }
    FPTbias.tab = (double**) malloc2d(7, FPTbias.N);
  }
  if (fdiff(cache[0], cosmology.random) || fdiff(cache[1], Ntable.random))
  {
    const double dlogk = (log(FPTbias.k_max) - log(FPTbias.k_min))/FPTbias.N;

    #pragma omp parallel for
    for (int i=0; i<FPTbias.N; i++) 
    {
      FPTbias.tab[5][i] = exp(log(FPTbias.k_min) + i*dlogk);
      FPTbias.tab[6][i] = p_lin(FPTbias.tab[5][i], 1.0);
    }

    double Pout[5][FPTbias.N];
    Pd1d2(FPTbias.tab[5], FPTbias.tab[6], FPTbias.N, Pout[0]);
    Pd2d2(FPTbias.tab[5], FPTbias.tab[6], FPTbias.N, Pout[1]);
    Pd1s2(FPTbias.tab[5], FPTbias.tab[6], FPTbias.N, Pout[2]);
    Pd2s2(FPTbias.tab[5], FPTbias.tab[6], FPTbias.N, Pout[3]);
    Ps2s2(FPTbias.tab[5], FPTbias.tab[6], FPTbias.N, Pout[4]);

    #pragma omp parallel for
    for (int i=0; i<FPTbias.N; i++) 
    {
      FPTbias.tab[0][i] = Pout[0][i]; // Pd1d2
      FPTbias.tab[1][i] = Pout[1][i]; // Pd2d2
      FPTbias.tab[2][i] = Pout[2][i]; // Pd1s2
      FPTbias.tab[3][i] = Pout[3][i]; // Pd2s2
      FPTbias.tab[4][i] = Pout[3][i]; // Pd2s2
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
    FPTIA.k_min     = 5.e-4;
    FPTIA.k_max     = 1.e+6;
    //FPTIA.N       = 1200 + 200 * Ntable.FPTboost;
    FPTIA.N         = 400 + 200 * Ntable.FPTboost;

    if (FPTIA.tab != NULL) {
      free(FPTIA.tab);
    }
    FPTIA.tab = (double**) malloc2d(12, FPTIA.N);
  }
  if (fdiff(cache[0], cosmology.random) || fdiff(cache[1], Ntable.random))
  {
    double lim[3];
    lim[0] = log(FPTIA.k_min);
    lim[1] = log(FPTIA.k_max);
    lim[2] = (lim[1] - lim[0])/FPTIA.N;
    
    #pragma omp parallel for
    for (int i=0; i<FPTIA.N; i++) 
    {
      FPTIA.tab[10][i] = exp(lim[0] + i*lim[2]);
      FPTIA.tab[11][i] = p_lin(FPTIA.tab[10][i], 1.0);
    }

    IA_tt(FPTIA.tab[10], FPTIA.tab[11], FPTIA.N, FPTIA.tab[0], FPTIA.tab[1]);
    
    IA_ta(FPTIA.tab[10], FPTIA.tab[11], FPTIA.N, FPTIA.tab[2], 
      FPTIA.tab[3], FPTIA.tab[4], FPTIA.tab[5]);
    
    IA_mix(FPTIA.tab[10], FPTIA.tab[11], FPTIA.N, FPTIA.tab[6], 
      FPTIA.tab[7], FPTIA.tab[8], FPTIA.tab[9]);
    
    #pragma omp parallel for
    for (int i=0; i<FPTIA.N; i++) {
      FPTIA.tab[7][i] *= 4.;
    }
    
    cache[0] = cosmology.random;
    cache[1] = Ntable.random;
  }
}
