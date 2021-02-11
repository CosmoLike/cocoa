#ifndef __COSMOLIKE_COSMO2D_FOURIER_H
#define __COSMOLIKE_COSMO2D_FOURIER_H
#ifdef __cplusplus
extern "C" {
#endif


// ------------------------------------------------------------------
// galaxy clustering power spectrum of galaxies in bins ni, nj
// ------------------------------------------------------------------
double C_cl_tomo(double l, int ni, int nj);

double C_cl_tomo_nointerp(double l, int ni, int nj);

// galaxy clustering power spectrum of galaxies in bin ni, using HOD model
double C_cl_HOD(double l, int ni);

// G-G lensing power spectrum, lens bin ni, source bin nj
double C_gl_tomo(double l, int ni, int nj);

double C_gl_tomo_nointerp(double l, int ni, int nj);

// G-G lensing power spectrum from HOD model, lens bin ni, source bin nj
double C_gl_HOD_tomo(double l, int ni, int nj);

// ------------------------------------------------------------------
// shear tomography power spectra
// ------------------------------------------------------------------
double C_shear_tomo(double l, int ni, int nj);

double C_shear_tomo_nointerp(double l, int ni, int nj);

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// CMB x LSS Functions
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

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