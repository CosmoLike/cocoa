#ifndef __COSMOLIKE_COSMO2D_FULLSKY_TATT_H
#define __COSMOLIKE_COSMO2D_FULLSKY_TATT_H
#ifdef __cplusplus
extern "C" {
#endif

double C_EE_TATT(double l, int ni,int  nj);

double C_BB_TATT(double l, int ni, int nj);

double C_EE_tab(double l, int ni, int nj);

double C_BB_tab(double l, int ni, int nj);

double C_ggl_TATT_tab(double l, int ni, int nj);

double w_gamma_t_TATT(int nt, int ni, int nj); 

double xi_pm_TATT(int pm, int nt, int ni, int nj); 

//double C1_TA(double a, double nz);

double b_TA(double a, double nz);

//double C2_TT(double a, double nz);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD