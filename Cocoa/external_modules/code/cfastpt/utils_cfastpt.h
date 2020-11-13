#ifndef __CFASTPT_UTILS_CFASTPT_H
#define __CFASTPT_UTILS_CFASTPT_H
#ifdef __cplusplus
extern "C" {
#endif

void extrap_log_linear(double *fk, int N_origin, int N_extra, double *large_fk);

int J_table(int *alpha_ar, int *beta_ar, int *l1_ar, int *l2_ar, int *l_ar, double *coeff_A_ar, int Nterms, int *alpha_ar_new, int *beta_ar_new, int *J1_ar, int *J2_ar, int *Jk_ar, double *coeff_AB_ar);

double coeff_B(int l1, int l2, int l, int J1, int J2, int Jk);

long factorial(int n);

double wigner_3j_jjj_000(int j1, int j2, int j3);

double Delta_coef(int a, int b, int c);

double Racah(int a, int b, int c, int d, int e, int f);

double wigner_6j(int j1, int j2, int j3, int j4, int j5, int j6);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD