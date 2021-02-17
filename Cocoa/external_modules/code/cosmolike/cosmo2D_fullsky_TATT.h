#ifndef __COSMOLIKE_COSMO2D_FULLSKY_TATT_H
#define __COSMOLIKE_COSMO2D_FULLSKY_TATT_H
#ifdef __cplusplus
extern "C" {
#endif

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Angular Power Spectra
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double C_EE_TATT(double l, int ni, int nj);

double C_BB_TATT(double l, int ni, int nj);

double C_ggl_TATT(double l, int nl, int ns);

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Tabulated Angular Power Spectra
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double C_EE_tab(double l, int ni, int nj);

double C_BB_tab(double l, int ni, int nj);

double C_ggl_TATT_tab(double l, int ni, int nj);

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// Correlation Functions (real space)
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

double w_gamma_t_TATT(int nt, int ni, int nj);

double xi_pm_TATT(int pm, int nt, int ni, int nj);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD