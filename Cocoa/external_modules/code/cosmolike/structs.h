#include <gsl/gsl_interp2d.h>

#ifndef __COSMOLIKE_STRUCTS_H
#define __COSMOLIKE_STRUCTS_H
#ifdef __cplusplus
extern "C" {
#endif

#define CHAR_MAX_SIZE 1024
#define MAX_SIZE_ARRAYS 10

typedef struct 
{
  double a_min;
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
  int N_theta;
  int N_thetaH;
  int N_S2;
  int N_DS;
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
  double Omega_m;  // matter density parameter
  double Omega_v;  // cosmogical constant parameter
  double h0;       // Hubble constant
  double Omega_nu; // Omega_m = Omega_cdm + Omega_nu (neutrinos) + omb
  double coverH0;  // units for comoving distances - speeds up code
  double rho_crit; // = 3 H_0^2/(8 pi G), critical comoving density
  double MGSigma;
  double MGmu;
  double random; // Random number between zero and 1 - see interface.cpp
  //double w0;
  //double A_s;
  //double n_s;
  int is_cached;
} cosmopara;

typedef struct 
{
  int shear_Nbin;         // number of source tomography bins
  int clustering_Nbin;    // number of lens galaxy bins
  int cluster_Nbin;       // number of lens cluster redshift bins

  int shear_Npowerspectra; // number of shear-shear tomography power spectra
  int ggl_Npowerspectra;   // number of galaxy-galaxy lensing tomography power spectra
  int clustering_Npowerspectra;  // number of galaxy-galaxy clustering tomography power spectra
  
  int cgl_Npowerspectra;             // number of cluster-galaxy lensing tomography combinations
  int cg_clustering_Npowerspectra;   // number of cluster-galaxy clustering tomography combinations
  int cc_clustering_Npowerspectra;   // number of cluster-cluster clustering tomography combinations

  double shear_zmax[MAX_SIZE_ARRAYS];
  double shear_zmin[MAX_SIZE_ARRAYS];
  double n_source[MAX_SIZE_ARRAYS];

  double clustering_zmax[MAX_SIZE_ARRAYS];
  double clustering_zmin[MAX_SIZE_ARRAYS];
  double n_lens[MAX_SIZE_ARRAYS];
  
  double cluster_zmax[MAX_SIZE_ARRAYS];
  double cluster_zmin[MAX_SIZE_ARRAYS];

  // we assume cluster bin = galaxy bin (no cross)
  int external_selection_cg_clustering[MAX_SIZE_ARRAYS];
} tomopara;

typedef struct
{
  int shear_photoz;
  double shear_zdistrpar_zmin;
  double shear_zdistrpar_zmax;
  int shear_histogram_zbins;
  char shear_REDSHIFT_FILE[CHAR_MAX_SIZE];

  int clustering_photoz;
  double clustering_zdistrpar_zmin;
  double clustering_zdistrpar_zmax;
  int clustering_histogram_zbins;
  char clustering_REDSHIFT_FILE[CHAR_MAX_SIZE];

  int clusters_photoz;
  char clusters_REDSHIFT_FILE[CHAR_MAX_SIZE];
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
} nuisancepara;

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
  double ggl_overlap_cut;
} sur;

typedef struct
{
  int Ncl;
  int Ntheta;
  int Ncos;
  int Ndata;
  double lmin;
  double lmax;
  double vtmax;
  double vtmin;
  double* ell;
  double* theta;
  double cosmax;
  double Rmin_bias;
  double Rmin_shear;
  double lmax_shear;
  double lmin_kappacmb;
  double lmax_kappacmb;
  int IA;
  int bias;
  int shear_shear;
  int shear_pos;
  int pos_pos;
  int gk;
  int kk;
  int ks;
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
} likepara;

typedef double (*B1_model)(double z, int nz);

typedef struct
{
  double b[MAX_SIZE_ARRAYS];      // linear galaxy bias paramter in clustering bin i
  double b2[MAX_SIZE_ARRAYS];     // quadratic bias parameter for redshift bin i
  double bs2[MAX_SIZE_ARRAYS];    // leading order tidal bias for redshift bin i
  double rcorr[MAX_SIZE_ARRAYS];
  double cg[MAX_SIZE_ARRAYS];
  double n_hod[MAX_SIZE_ARRAYS];
  double b_mag[MAX_SIZE_ARRAYS];  // amplitude of magnification bias, b_mag[i] = 5*s[i]+beta[i]-2
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
  double fwhm;        // beam fwhm in rad
  double sensitivity; // white noise level in muK*rad
  char pathLensRecNoise[CHAR_MAX_SIZE]; // path to precomputed noise on reconstructed kappa
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

extern barypara bary;

extern Cmb cmb;

extern lim limits;

extern Ntab Ntable;

extern pre precision;

extern TinkerEmuParameters tinkerEmuParam;

//  ----------------------------------------------------------------------------------
//  ----------------------------------------------------------------------------------
//  ----------------------------------------------------------------------------------
//  RESET STRUCTS 
//  ----------------------------------------------------------------------------------
//  ----------------------------------------------------------------------------------
//  ----------------------------------------------------------------------------------

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

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD
