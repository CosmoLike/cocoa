#ifndef __COSMOLIKE_REDSHIFT_SPLINE_H
#define __COSMOLIKE_REDSHIFT_SPLINE_H
#ifdef __cplusplus
extern "C" {
#endif

// ------------------------------------------------------------------------
// integration boundary routines
// ------------------------------------------------------------------------

double amin_source(int i);

double amax_source(int i);

double amax_source_IA(int i);

double amin_lens(int i);

double amax_lens(int i);

// ------------------------------------------------------------------------
//redshift overlap tests, allowed tomography combinations
// ------------------------------------------------------------------------

// test whether the (l,zl) bin is in the linear
// clustering regime - return 1 if true, 0 otherwise
int test_kmax(double l, int zl);

// test whether source bin zs is behind lens bin zl
int test_zoverlap_cov(int zl, int zs);

// test whether source bin zs is behind lens bin zl
int test_zoverlap(int zl, int zs);

// test whether source bin zs is behind lens bin zl
int test_zoverlap_c(int zc, int zs);

// (z_l,z_s) -> N_tomo_ggl
int N_ggl(int zl, int zs);

void write_gglensing_zbins(char *surveyname);

int ZL(int Nbin);

int ZS(int Nbin);

int N_cgl(int zl, int zs);

int ZC(int Nbin);

int ZSC(int Nbin);

// find shear tomography bin number N_shear of tomography combination (z1,z2)
int N_shear(int z1, int z2);

// find z1 of tomography combination (z1,z2) constituting
// shear tomography bin Nbin
int Z1(int Nbin);

// find z2 of tomography combination (z1,z2) constituting
// shear tomography bin Nbin
int Z2(int Nbin);

// find zcl1 of tomography combination (zcl1,zcl2)
// constituting galaxy clustering tomography bin Nbin
int Zcl1(int Nbin);

// find zcl2 of tomography combination (zcl1,zcl2)
// constituting galaxy clustering tomography bin Nbin
int Zcl2(int Nbin);

// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// routines for redshift distributions, including photo-zs * (optional)
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------

// ------------------------------------------------------------------------
// begin shear routines
// ------------------------------------------------------------------------

// returns n(ztrue | j), works only with binned distributions;
// j =-1 -> no tomography; j>= 0 -> tomography bin j
double zdistr_photoz(double zz, int j);

// double zdistr_gaussian_ztrue (double zz, int j);

// ------------------------------------------------------------------------
// clustering routines
// ------------------------------------------------------------------------

// returns n(ztrue | j), works only with binned distributions;
// j =-1 -> no tomography; j>= 0 -> tomography bin j
double pf_photoz(double zz, int j);

// ------------------------------------------------------------------------
// specific photo-z models
// ------------------------------------------------------------------------

// compute photo-z scatter (sigma(z)) from nuisance.sigma_zphot_shear parameters
double sigma_zphot_shear(double z, int nz);

// compute photo-z bias (b_{zph}(z)) from nuisance.bias_zphot_shear parameters
double bias_zphot_shear(double z, int nz);

// ------------------------------------------------------------------------
// end begin shear routines
// ------------------------------------------------------------------------

double sigma_zphot_clustering(double z, int nz);

double bias_zphot_clustering(double z, int nz);

// ------------------------------------------------------------------------
// routines calculating the number of source and lens galaxies per bin
// ------------------------------------------------------------------------

double int_nsource(double z, void *param);

double int_nlens(double z, void *param);

// returns n_gal for shear tomography bin j,  works only with binned
// distributions; j =-1 -> no tomography; j>= 0 -> tomography bin j
double nsource(int j);

// returns n_gal for clustering tomography bin j, works only with binned
// distributions; j =-1 -> no tomography; j>= 0 -> tomography bin j
double nlens(int j);

// mean true redshift of (clustering/lens) galaxies within redshift bin j
double zmean(int j);

// mean true redshift of source galaxies in tomography bin j
double zmean_source(int j);

// ------------------------------------------------------------------------
// lens efficiencies
// ------------------------------------------------------------------------

// lens efficiency for CMB lensing
double g_cmb(double a);

// lens efficiency of source galaxies in tomography bin zbin
double g_tomo(double a, int zbin);

// lens efficiency of source galaxies in tomography bin zbin
double g2_tomo(double a, int zbin);

// lens efficiency of *lens* galaxies in tomography bin zbin -
// used for magnification calculations
double g_lens(double a, int zbin);

// ------------------------------------------------------------------------
// routines for redshift distributions, including photo-zs * (optional)
// ------------------------------------------------------------------------

// return nz(z) based on redshift file with one redshift distribution
double zdistr_histo_1(double z, void *params);

// ------------------------------------------------------------------------
// routines for association of a pair redshift bin numbers
// and power spectrum tomography bin
// ------------------------------------------------------------------------

double ggl_efficiency(int zl, int zs);

// no longer supported - declaration only to prevent compile errors
double g_bg(double a, int nzlens);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD