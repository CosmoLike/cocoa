#ifndef __COSMOLIKE_CMBLSSFOURIER_H
#define __COSMOLIKE_CMBLSSFOURIER_H
#ifdef __cplusplus
extern "C" {
#endif

 //lensing efficiency weightfunction for CMB lensing
double W_k(double a, double fK);

double int_for_C_gk(double a, void *params);

double int_for_C_gk_b2(double a, void *params);

double int_for_C_ks(double a, void *params);

double int_for_C_kk(double a, void *params);

double int_for_C_ks_IA_mpp(double a, void *params);

double C_ks_IA(double s, int ni);

double C_gk_nointerp(double l, int nl);

double C_ks_nointerp(double l, int ns);

double C_kk_nointerp(double l);

double C_gk(double l, int ni);

double C_ks(double l, int ni);

double C_kk(double l);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD