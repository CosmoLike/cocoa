#include <gsl/gsl_interp2d.h>

#ifndef __COSMOLIKE_STRUCTS_H
#define __COSMOLIKE_STRUCTS_H
#ifdef __cplusplus
extern "C" {
#endif

#define CHAR_MAX_SIZE 1024
#define MAX_SIZE_ARRAYS 10

#define NO_IA 0
#define IA_NLA_LF 1
#define IA_REDSHIFT_BINNING 2
#define IA_REDSHIFT_EVOLUTION 3

#define IA_MODEL_NLA 0
#define IA_MODEL_TATT 1


typedef struct 
{
  double a_min;
  double a_min_hm;
  double k_min_cH0;
  double k_max_cH0;
  double M_min;
  double M_max;
  int LMIN_tab;
  int LMAX_NOLIMBER;
  int LMAX;
  double LMIN_hankel;
  double LMAX_hankel;
  double cluster_util_log_M_min;
  double cluster_util_log_M_max;
  double binned_P_lambda_obs_given_M_zmin_table;
  double binned_P_lambda_obs_given_M_zmax_table;
  double SDSS_P_lambda_obs_given_M_true_lambda_min;
  double SDSS_P_lambda_obs_given_M_true_lambda_max;
  double halo_exclusion_k_min;
  double halo_exclusion_k_max;
  double halo_exclusion_k_min_hankel;
  double halo_exclusion_k_max_hankel;
  double halo_exclusion_R_min;
  double halo_exclusion_R_max;
  double halo_uKS_cmin; // halo.c u_KS(double c, double k, double rv)
  double halo_uKS_cmax; // halo.c u_KS(double c, double k, double rv)
  double halo_uKS_xmin; // halo.c u_KS(double c, double k, double rv)
  double halo_uKS_xmax; // halo.c u_KS(double c, double k, double rv)
} lim;

typedef struct 
{
  double low;
  double medium;
  double high;
  double insane;
} pre;

typedef struct 
{
  int N_a;          
  int N_k_lin;
  int N_k_nlin;
  int N_ell;
  int Ntheta;
  int N_thetaH;
  int N_M;
  int N_ell_TATT;                  // Cosmo2D
  int NL_Nell_block;               // Cosmo2D - NL = NonLimber
  int NL_Nchi;                     // Cosmo2D - NL = NonLimber
  int N_a_halo_exclusion;          // N_a for binned_p_cc_incl_halo_exclusion (cluster_util.c)
  int N_k_halo_exclusion;          // N_k for binned_p_cc_incl_halo_exclusion (cluster_util.c)
  int N_k_hankel_halo_exclusion;   // N for 3D Hankel Transform (pk_to_xi and xi_to_pk)
  int N_R_halo_exclusion;
  int binned_P_lambda_obs_given_M_size_z_table;
  int binned_P_lambda_obs_given_M_size_M_table;
  int binned_p_cm_size_a_table;
  int halo_uKS_nc;                 // halo.c u_KS(double c, double k, double rv)
  int halo_uks_nx;                 // halo.c u_KS(double c, double k, double rv)
  int acc_boost_photoz_sampling;
  int photoz_interpolation_type;
  int high_def_integration;
  double vtmax;
  double vtmin;
} Ntab;

typedef struct 
{
  double xmin;
  double xmax;
  double Pmin;
  double Pmax;
  double dPmin;
  double dPmax;
} bin_avg;

typedef struct
{ // parameters for power spectrum passed to FASTPT
  double k_min;
  double k_max;
  int N;
  int N_per_dec;
  double** tab_AB; // parameters for table of bias terms
  int N_AB;        // parameters for table of bias terms
  double** tab_IA; // parameters for table of IA terms
  int N_IA;        // parameters for table of IA terms
} FPTpara;

typedef struct
{
  double Omega_b;  // baryon density paramter
  double Omega_m;  // matter density parameter
  double Omega_v;  // cosmogical constant parameter
  double h0;       // Hubble constant
  double Omega_nu; // Omega_m = Omega_cdm + Omega_nu (neutrinos) + omb
  double coverH0;  // units for comoving distances - speeds up code
  double rho_crit; // = 3 H_0^2/(8 pi G), critical comoving density
  double MGSigma;
  double MGmu;
  double sigma_8;
  double random; // Random number between zero and 1 - see interface.cpp
  //double w0;
  //double A_s;
  //double n_s;
  int is_cached;
} cosmopara;

typedef struct 
{
  int shear_Npowerspectra;       // number of shear-shear tomography power spectra
  int ggl_Npowerspectra;         // number of galaxy-galaxy lensing tomography power spectra
  int clustering_Npowerspectra;  // number of galaxy-galaxy clustering tomography power spectra
} tomopara;

typedef struct
{
  int shear_nbin;         // number of source tomography bins
  int shear_photoz;
  int shear_nzbins;
  double** shear_zdist_table;
  double shear_zdist_zmin_all;
  double shear_zdist_zmax_all;
  double shear_zdist_zmin[MAX_SIZE_ARRAYS];
  double shear_zdist_zmax[MAX_SIZE_ARRAYS];

  int clustering_nbin;    // number of lens galaxy bins
  int clustering_photoz;
  int clustering_nzbins;
  double** clustering_zdist_table;
  double clustering_zdist_zmin_all;
  double clustering_zdist_zmax_all;
  double clustering_zdist_zmin[MAX_SIZE_ARRAYS];
  double clustering_zdist_zmax[MAX_SIZE_ARRAYS];

  /*
  int cluster_Nbin;       // number of lens cluster redshift bins
  int clusters_photoz;
  char clusters_REDSHIFT_FILE[CHAR_MAX_SIZE];
  */
  /*
  int cgl_Npowerspectra;             // number of cluster-galaxy lensing tomography combinations
  int cg_clustering_Npowerspectra;   // number of cluster-galaxy clustering tomography combinations
  int cc_clustering_Npowerspectra;   // number of cluster-cluster clustering tomography combinations 
  double cluster_zmax[MAX_SIZE_ARRAYS];
  double cluster_zmin[MAX_SIZE_ARRAYS];
  // we assume cluster bin = galaxy bin (no cross)
  int external_selection_cg_clustering[MAX_SIZE_ARRAYS];
  */
} redshiftpara;

typedef struct
{
  int interpolate_survey_area;
  int bias_model;                 // Bias model
  int hmf_model;                  // HMF model 
  int nonlinear_bias;             // Do we include nonlinear bias in cluster analysis?

  int N_MOR;                      // Mass observable relation (number of nuisance params)
  int N_SF;                       // selection function (number of nuisance params)

  int halo_exclusion_model;
  double delta_exclusion;  // delta for exclusion radius (halo_exclusion) according to Baldauf 2013

  int N200_Nbin;                  // number of cluster bins in lambda_obs (observed richness)
  double N200_min;                // global lambda_obs_min (observed richness)
  double N200_max;                // global lambda_obs_max (observed richness)
  double N_min[MAX_SIZE_ARRAYS];  // lambda_obs_min in each bin in lambda_obs (observed richness)
  double N_max[MAX_SIZE_ARRAYS];  // lambda_obs_max in each bin in lambda_obs (observed richness)
  
  char model[CHAR_MAX_SIZE];
} clusterparams;

typedef struct
{
  // like.IA = 1 variables
  double LF_alpha;
  double LF_P;
  double LF_Q;
  double LF_red_alpha;
  double LF_red_P;
  double LF_red_Q;
  
  // like.IA = 3; NLA; = 4; NLA, = 5; TATT (per bin);  = 6: TATT (power law)
  double A_z[MAX_SIZE_ARRAYS];    // NLA normalization per source redshift bin
  double A2_z[MAX_SIZE_ARRAYS];   // NLA normalization per source redshift bin
  double b_ta_z[MAX_SIZE_ARRAYS]; // b_ta, per bin (like.IA = 6), or use b_ta_z[0] with like.IA = 5
  double A_ia;                    // A IA see Joachimi2012
  double A2_ia;                   // placeholder param for quadratic,etc IA
  double beta_ia;                 // beta IA see Joachimi2012
  double eta_ia;                  // eta_other IA see Joachimi2012
  double eta_ia_tt;               // same as eta_ia, for TT
  double eta_ia_highz;            // uncertainty in high z evolution
  double oneplusz0_ia;            // oneplusz0-ia MegaZ
  double c1rhocrit_ia;
  double fred[MAX_SIZE_ARRAYS];
  
  double shear_calibration_m[MAX_SIZE_ARRAYS];

  double sigma_zphot_shear[MAX_SIZE_ARRAYS];
  double bias_zphot_shear[MAX_SIZE_ARRAYS];

  double sigma_zphot_clustering[MAX_SIZE_ARRAYS];
  double bias_zphot_clustering[MAX_SIZE_ARRAYS];
  
  // Variables for the 4x2pt+N (see: 2008.10757 & 2010.01138)
  double cluster_MOR[MAX_SIZE_ARRAYS];
  double cluster_selection[MAX_SIZE_ARRAYS];

  double frac_lowz;   // photo-z = 5 model (see: 2108.00658)
  double frac_highz;  // photo-z = 5 model (see: 2108.00658)
} nuisancepara;

typedef struct
{
  // Compton-Y related variables
  double gas_Gamma_KS; // Gamma in K-S profile
  double gas_beta;     // beta: mass scaling index in bound gas fraction
  double gas_lgM0;     // critical halo mass, below which gas ejection is significant
  double gas_eps1;
  double gas_eps2;
  double gas_alpha;
  double gas_A_star;
  double gas_lgM_star;
  double gas_sigma_star;
  double gas_lgT_w;
  double gas_f_H;
} ynuisancepara;

typedef struct
{
  double area;                    // survey_area in deg^2.
  double n_gal;                   // galaxy density per arcmin^2
  double sigma_e;                 // rms inrinsic ellipticity noise
  double area_conversion_factor;  // factor from deg^2 to radian^2:
  double n_gal_conversion_factor; // factor from n_gal/arcmin^2 to n_gal/radian^2:
  double n_lens;                  // lens galaxy density per arcmin^2
  double m_lim;
  char name[CHAR_MAX_SIZE];
} sur;

typedef struct
{
  int is_cmb_bandpower;
  int is_cmb_kkkk_cov_from_sim;// if kkkk covmat is from sim, apply Hartlap factor
  double alpha_Hartlap_kkkk;
  int Ncl;
  int Ncos;
  int Ndata;
  int Nbp; // number of band-power bins for CMB lensing
  double lmin;
  double lmax;
  double* ell;
  double* theta;
  double cosmax;
  double Rmin_bias;
  double Rmin_shear;
  double lmax_shear;
  double lmin_bp;      // minimum l mode of CMB band power binning matrix (w/ corr)
  double lmax_bp;      // Note: for dv, only need binmat with corr
  double lmin_kappacmb;
  double lmax_kappacmb;
  int IA;
  int IA_MODEL;
  int bias;
  int shear_shear;
  int shear_pos;
  int pos_pos;
  int gk;
  int kk;
  int ks;
  int gy;
  int sy;
  int ky;
  int yy;
  int clusterN;
  int clusterWL;
  int clusterCG;
  int clusterCC;
  int use_full_sky_shear;
  int use_full_sky_ggl;
  int use_full_sky_clustering;
  int use_full_sky_gk;
  int use_full_sky_sk;
  int high_def_integration;
  int adopt_limber_gg;
  int adopt_limber_gammat;
  int use_ggl_efficiency_zoverlap;
} likepara;

typedef double (*B1_model)(double z, int nz);

typedef struct
{
  double b[MAX_SIZE_ARRAYS];      // linear galaxy bias paramter in clustering bin i
  double b2[MAX_SIZE_ARRAYS];     // quadratic bias parameter for redshift bin i
  double bs2[MAX_SIZE_ARRAYS];    // leading order tidal bias for redshift bin i
  double b_mag[MAX_SIZE_ARRAYS];  // amplitude of magnification bias, b_mag[i] = 5*s[i]+beta[i]-2
  
  double hod[MAX_SIZE_ARRAYS][MAX_SIZE_ARRAYS]; // HOD[i] contains HOD parameters of galaxies in 
                                                // clustering bin i, following 5 parameter model 
                                                // of Zehavi et al. 2011 + modification of 
                                                // concentration parameter
  double cg[MAX_SIZE_ARRAYS];                   // galaxy concentration parameter
  B1_model b1_function;
} galpara;

typedef struct
{
  char runmode[CHAR_MAX_SIZE];
} pdeltapara;

typedef struct
{
  int is_Pk_bary;
  int Na_bins;
  int Nk_bins;
  double* a_bins;
  double* logk_bins;
  double* log_PkR;
  gsl_interp2d_type* T;
  gsl_interp2d* interp2d;
} barypara;

typedef struct
{
  char name[CHAR_MAX_SIZE];
  double fwhm;                            // beam fwhm in rad
  double sensitivity;                     // white noise level in muK*rad
  char pathLensRecNoise[CHAR_MAX_SIZE];   // path to precomputed noise on reconstructed kappa
  char pathHealpixWinFunc[CHAR_MAX_SIZE]; // path to precomputed healpix window function
} Cmb;

typedef struct 
{
  int tinker_bias_ncosmo;
  int tinker_bias_nparam;
  int tinker_bias_nsamp;
  int tinker_bias_nparam_redshift;
  int tinker_bias_extrapolation_cut_in;
  int tinker_bias_extrapolation_cut_out;
  int tinker_hmf_ncosmo;
  int tinker_hmf_nparam;
  int tinker_hmf_nsamp;
  int tinker_hmf_nparam_redshift;
  int tinker_hmf_extrapolation_cut_in;
  int tinker_hmf_extrapolation_cut_out;
} TinkerEmuParameters;

double bgal_z(double z, int nz);

double b1_per_bin(double z, int nz);

extern likepara like;

extern cosmopara cosmology;

extern tomopara tomo;

extern redshiftpara redshift;

extern sur survey;

extern galpara gbias;

extern clusterparams Cluster;

extern pdeltapara pdeltaparams;

extern FPTpara FPT;

extern nuisancepara nuisance;

extern ynuisancepara ynuisance;

extern barypara bary;

extern Cmb cmb;

extern lim limits;

extern Ntab Ntable;

extern pre precision;

extern TinkerEmuParameters tinkerEmuParam;

// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------
//  RESET STRUCTS 
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------

void reset_bary_struct();

// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// UPDATE
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------

void update_cosmopara(cosmopara* C);

void update_galpara(galpara* G);

void update_nuisance(nuisancepara* N);

void update_ynuisance(ynuisancepara* N);

void update_table(Ntab* N);

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD
