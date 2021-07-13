#ifndef __COSMOLIKE_REDSHIFT_SPLINE_H
#define __COSMOLIKE_REDSHIFT_SPLINE_H
#ifdef __cplusplus
extern "C" {
#endif

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// integration boundary routines
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

double amin_source(int i);

double amax_source(int i);

double amax_source_IA(int i);

double amin_lens(int i);

double amax_lens(int i);

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
//redshift overlap tests, allowed tomography combinations
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

int test_kmax(double l, int zl); // test whether the (l,zl) bin is in the linear clustering regime

int test_zoverlap(int zl, int zs); // test whether source bin zs is behind lens bin zl

int test_zoverlap_c(int zc, int zs); // test whether source bin zs is behind lens bin zl (clusters)

int ZL(int Nbin); // find zlens bin of tomo combination (galaxy-galaxy lensing)
int ZS(int Nbin); // find zsource bin of tomo combination (galaxy-galaxy lensing)
int N_ggl(int zl, int zs); // find tomo bin number tomography combination

int ZC(int Nbin); // find zlens bin of tomo combination (cluster-galaxy lensing)
int ZSC(int Nbin); // find zsource bin of tomo combination (cluster-galaxy lensing)
int N_cgl(int zl, int zs); // find tomo bin number tomography combination

int Z1(int Nbin); // find z1 bin of tomo combination 
int Z2(int Nbin); // find z2 bin of tomo combination 
int N_shear(int z1, int z2); // find tomo bin number tomography combination

int ZCL1(int Nbin); // find z1 bin of tomo combination (cluster-galaxy clustering)
int ZCL2(int Nbin); // find z2 bin of tomo combination (cluster-galaxy clustering)
int N_CL(int z1, int z2); // find tomo bin number tomography combination (cluster-galaxy clustering)

int ZCCL1(int Nbin); // find z1 bin of tomo combination (cluster-cluster clustering)
int ZCCL2(int Nbin); // find z2 bin of tomo combination (cluster-cluster clustering)
int N_CCL(int z1, int z2); // find tomo bin number tomography combination (cluster-cluster clustering)

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// shear routines (for redshift distributions, including photo-zs * (optional))
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

// returns n(ztrue | j), works only with binned distributions; 
double zdistr_photoz(double z, int j); // j =-1 -> no tomography; j>= 0 -> tomography bin j

double sigma_zphot_shear(double z, int nz); // compute photo-z scatter (sigma(z))

double bias_zphot_shear(double z, int nz); // compute photo-z bias

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// clustering routines (for redshift distributions, including photo-zs * (optional))
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

// returns n(ztrue | j), works only with binned distributions;
double pf_photoz(double zz, int j); // j =-1 -> no tomography; j>= 0 -> tomography bin j

double sigma_zphot_clustering(double z, int nz); // compute photo-z scatter (sigma(z))

double bias_zphot_clustering(double z, int nz); // compute photo-z bias

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// routines calculating the number of source and lens galaxies per bin
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

// returns n_gal for shear tomography bin j, works only with binned
double nsource(int j); // distributions; j =-1 -> no tomography; j>= 0 -> tomography bin j

// returns n_gal for clustering tomography bin j, works only with binned
double nlens(int j); // distributions; j =-1 -> no tomography; j>= 0 -> tomography bin j

double zmean(int j); // mean true redshift of (clustering/lens) galaxies within redshift bin j

double zmean_source(int j); // mean true redshift of source galaxies in tomography bin j

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// lens efficiencies
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------


double g_cmb(double a); // lens efficiency for CMB lensing

double g_tomo(double a, int zbin); // lens efficiency of source galaxies in tomography bin zbin

double g2_tomo(double a, int zbin); // lens efficiency of source galaxies in tomography bin zbin

// lens efficiency of *lens* galaxies in tomography bin zbin - used for magnification calculations
double g_lens(double a, int zbin);

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// routines for association of a pair redshift bin numbers and power spectrum tomography bin
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

double ggl_efficiency(int zl, int zs);

// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// SELECTION FUNCTION (Cluster)
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

// simplfied selection function, disregards evolution of N-M relation+mass function within z bin
double zdistr_cluster(int nz,  double z, double chi, double hoverh0) ;


#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD