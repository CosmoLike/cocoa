#ifndef __COSMOLIKE_IA_H
#define __COSMOLIKE_IA_H
#ifdef __cplusplus
extern "C" {
#endif

double int_for_C_II(double a, void *params);

double C_II_nointerp(double s, int ni, int nj);

double int_for_C_II_lin(double a, void *params);

double C_II_lin_nointerp(double s, int ni, int nj);

double int_for_C_gI(double a, void *params);

double C_gI_nointerp(double s, int ni, int nj);

double int_for_C_gI_lin(double a, void *params);

double C_gI_lin_nointerp(double s, int ni, int nj);

double int_for_C_GI(double a, void *params);

double C_GI_nointerp(double s, int ni, int nj);

double int_for_C_GI_lin(double a, void *params);

double C_GI_lin_nointerp(double s, int ni, int nj);

double int_for_C_ggl_IA_mpp(double a, void *params);

double int_for_C_ggl_IA_mpp_b2(double a, void *params);

double int_for_C_shear_shear_IA_mpp(double a, void *params);

double C_shear_shear_IA(double s, int ni, int nj);

double C_ggl_IA(double s, int nl, int ns);

double C_ggl_IA_tab(double l, int ni, int nj);

double C_shear_shear_IA_tab(double l, int ni, int nj);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD