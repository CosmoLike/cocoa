#include "clik_parametric.h"
#include "clik_parametric_addon.h"
#define PRM_NU0 143.

void gal_beta_dust_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err);
void gal_T_dust_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err);
void gal_alpha_non_thermal_derivative(parametric *egl, int iv, double *Rq, double *dRq, error **err);
double dust_spectrum(double nu, double T_dust, double beta_dust, double nu0);
double non_thermal_spectrum(double nu, double alpha_non_thermal, double nu0);