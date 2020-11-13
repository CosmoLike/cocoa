#ifndef __COSMOLIKE_COSMO2D_EXACT_FFT_H
#define __COSMOLIKE_COSMO2D_EXACT_FFT_H
#ifdef __cplusplus
extern "C" {
#endif

void f_chi_for_Psi_cl(double *chi_ar, int Nchi, double *f_chi_ar, int ni);

void f_chi_for_Psi_cl_RSD(double *chi_ar, int Nchi, double *f_chi_RSD_ar,
                          int ni);

void f_chi_for_Psi_cl_Mag(double *chi_ar, int Nchi, double *f_chi_Mag_ar,
                          int ni);

void C_cl_mixed(int L, int LMAX, int ni, int nj, double *Cl, double dev,
                double tolerance);

double w_tomo_nonLimber(int nt, int ni,
                        int nj); // w(theta) including non-Limber+RSD

void f_chi_for_Psi_sh(double *chi_ar, int Nchi, double *f_chi_ar, int nz);

void f_chi_for_Psi_sh_IA(double *chi_ar, int Nchi, double *f_chi_IA_ar, int nz);

void C_gl_mixed(int L, int LMAX, int ni, int nj, double *Cl, double dev,
                double tolerance);

double w_gamma_t_nonLimber(int nt, int ni, int nj);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD