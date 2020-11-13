#ifndef __COSMOLIKE_HALO_H
#define __COSMOLIKE_HALO_H
#ifdef __cplusplus
extern "C" {
#endif

// ----------------------------------------------------------------------
// relations for converting mass -> radius, using \Delta = 200 \rho_m for
// consistency with mass & bias function
// ----------------------------------------------------------------------
double m_Delta(double r, double a);

double r_Delta(double m, double a);

double bias_norm(double a);

// fitting functions for mass and bias
double massfunc(double m, double a);

double bias_norm(double a);

double mass_norm(double a);

double B1(double m, double a);

// ----------------------------------------------------------------------
// halo density profiles
// ----------------------------------------------------------------------

// mass-concentration relation
double conc(double m, double a);

// analytic expression for Fourier transform of NFW density profile
double u_nfw_c(double c, double k, double m, double a);

// Fourier transform of NFW density profile, truncated
// at r_Delta, through direct integration
double u_nfw(double c, double k, double m, double a);

// use this routine in inner I[0-1]j and modify int_rho_nfw if to work with
// different halo profiles, e.g. to include adiabatic contraction, AGN feedback,
// etc.

// halo model building blocks
double I0j(int j, double k1, double k2, double k3, double k4, double a);

double I_11(double k, double a);

double I1j(int j, double k1, double k2, double k3, double a);

/*halo model matter power spectrum, bispectrum, trispectrum*/
double p_1h(double k, double a);

double p_2h(double k, double a);

double Pdelta_halo(double k, double a);

/*look up table for 1-h halo sample variance term */
double I12_SSC(double k, double a);

// ----------------------------------------------------------------------
// simple abundance matching routines for source galaxies
// ----------------------------------------------------------------------

// comoving source galaxy density based on n_gal +n(z)
double n_s_cmv(double a);

// simple abundance matching routine to get bias for
// source galaxies from redshift distribution
double b_ngmatched(double a, double n_cmv);

double b_source(double a); // lookup table for b1 of source galaxies

double radius(double m);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD