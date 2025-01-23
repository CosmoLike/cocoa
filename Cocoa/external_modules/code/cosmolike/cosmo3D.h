#include <stdbool.h>

#ifndef __COSMOLIKE_COSMO3D_H
#define __COSMOLIKE_COSMO3D_H
#ifdef __cplusplus
extern "C" {
#endif

struct chis {
  double chi;
  double dchida;
};

struct growths {
  double D;
  double f;
};

struct growths growfac_all(double io_a);

// return chi and dchi_da in a single spline call (faster)
struct chis chi_all(double io_a);

double norm_growfac(const double a, const bool normalize_z0);

double f_growth(const double z);

double growfac(double a);

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

double f_K(double chi);

double chi(double a);

double dchi_da(double a);

double hoverh0(double a);

double hoverh0v2(double io_a, double dchida);

double a_chi(double chi1);

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

double PkRatio_baryons(double kintern, double a); // k in h/Mpc

void set_run_mode(int ps_mode);

double p_lin(const double k, const double a);

double p_nonlin(const double k, const double a);

double Pdelta(double k_NL, double a); // k in coverH0 units

double MG_Sigma(double a);

double PkRatio_baryons(double k_NL, double a);

double sigma2(const double M);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD