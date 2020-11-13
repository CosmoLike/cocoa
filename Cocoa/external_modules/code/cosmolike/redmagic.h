#ifndef __COSMOLIKE_REDMAGIC_H
#define __COSMOLIKE_REDMAGIC_H
#ifdef __cplusplus
extern "C" {
#endif

// ------------------------------------------------------------------------
// HOD model routines
// ------------------------------------------------------------------------

// central galaxy occupation as a function of halo mass (in Msun/h) in bin nz at
// scale factor a
double n_c_rm(double mh);

// satellite galaxy occupation function
double n_s_rm(double mh);

// Fourier transformed radial galaxy density profile
double u_g_rm(double m, double k, double a);

// galaxy number density of HOD model  in bin nz at scale factor a (in c/H0 = 1
// units)
double ngal_rm(double a);

// mean galaxy bias of HOD model in bin nz at scale factor a
double bgal_rm(double a);

// mean host halo mass of HOD model in bin nz at scale factor a
double mmean_rm(double a);

// mean host halo mass of central HOD model in bin nz at scale factor a
double mc_rm(double a);

// satellite fraction of HOD model
double fsat_rm(double a);

double C_cl_HOD_rm_tomo(double l, int nz);

double C_gl_HOD_rm_tomo(double l, int nzl, int nzs);

// ------------------------------------------------------------------------
/// HOD model galaxy power spectrum + galaxy-matter cross spectrum
// ------------------------------------------------------------------------

// galaxy-galaxy power spectrum based on redmagic HOD model
double P_gg_rm(double k, double a);

// matter-galaxy power spectrum based on redmagic HOD model
double P_gm_rm(double k, double a);

// ------------------------------------------------------------------------
// ------------------------------------------------------------------------

double G02_rm(double k, double a);

double int_bgal_rm(double m, void *params);

double int_ngal_rm(double m, void *params);

double int_for_C_cl_HOD_rm(double a, void *params);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD