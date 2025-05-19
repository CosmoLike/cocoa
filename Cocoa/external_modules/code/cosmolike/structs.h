#include <gsl/gsl_interp2d.h>

#ifndef __COSMOLIKE_STRUCTS_H
#define __COSMOLIKE_STRUCTS_H
#ifdef __cplusplus
extern "C" {
#endif

#define CHAR_MAX_SIZE 1024
#define MAX_SIZE_ARRAYS 20

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
  // ---------------------------------------------------
  // CLUSTER ROUTINES (ALPHA STAGE)
  // ---------------------------------------------------
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
  // ---------------------------------------------------
  // ---------------------------------------------------
  // CACHE VARIABLES
  // ---------------------------------------------------
  // ---------------------------------------------------
  double random; 
  // ---------------------------------------------------
  // ---------------------------------------------------
  // CONTROL NUM POINTS ON COSMOLIKE TABLES
  // ---------------------------------------------------
  // ---------------------------------------------------
  int N_a;          
  int N_k_lin;
  int N_k_nlin;
  int N_ell;
  int Ntheta;
  int N_M;
  int NL_Nell_block;               // Cosmo2D - NL = NonLimber
  int NL_Nchi;                     // Cosmo2D - NL = NonLimber
  // ---------------------------------------------------
  // ---------------------------------------------------
  // INTERPOL TYPE ON SOURCE/LENS n(z) INTERPOLATION
  // ---------------------------------------------------
  // ---------------------------------------------------
  int photoz_interpolation_type;
  // ---------------------------------------------------
  // ---------------------------------------------------
  // THETA RANGE ON REAL SPACE CORRELATION FUNCTIONS
  // ---------------------------------------------------
  // ---------------------------------------------------
  double vtmax;
  double vtmin;
  // ---------------------------------------------------
  // ---------------------------------------------------
  // CONTROL NUM POINTS EVALUATED ON COSMOLIKE INTEGRALS
  // ---------------------------------------------------
  // ---------------------------------------------------
  int high_def_integration;
  // ---------------------------------------------------
  // ---------------------------------------------------
  // CONTROL NUM POINTS EVALUATED BY FASPT
  // ---------------------------------------------------
  // ---------------------------------------------------
  int FPTboost;
  // ---------------------------------------------------
  // ---------------------------------------------------
  // CLUSTER ROUTINES (ALPHA STAGE)
  // ---------------------------------------------------
  // ---------------------------------------------------
  int N_a_halo_exclusion;        // N_a for binned_p_cc_incl_halo_exclusion (cluster_util.c)
  int N_k_halo_exclusion;        // N_k for binned_p_cc_incl_halo_exclusion (cluster_util.c)
  int N_k_hankel_halo_exclusion; // N for 3D Hankel Transform (pk_to_xi and xi_to_pk)
  int N_R_halo_exclusion;
  int binned_P_lambda_obs_given_M_size_z_table;
  int binned_P_lambda_obs_given_M_size_M_table;
  int binned_p_cm_size_a_table;
  int halo_uKS_nc;               // halo.c u_KS(double c, double k, double rv)
  int halo_uks_nx;               // halo.c u_KS(double c, double k, double rv)
} Ntab;

typedef struct
{
  // ---------------------------------------------------
  // ---------------------------------------------------
  // CACHE VARIABLES
  // ---------------------------------------------------
  // ---------------------------------------------------
  double random;
  // ---------------------------------------------------
  // ---------------------------------------------------
  // COSMO PARAMETERS
  // ---------------------------------------------------
  // ---------------------------------------------------
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
  // ---------------------------------------------------
  // ---------------------------------------------------
  // MATTER POWER SPECTRUM
  // size = (lnP_nk+1,lnP_nz+1)
  // z = lnP[lnP_nk,j<lnP_nz]
  // k = lnP[i<lnP_nk,lnP_nz]
  // ---------------------------------------------------
  // ---------------------------------------------------
  int lnP_nk;
  int lnP_nz;
  double** lnP;
  // ---------------------------------------------------
  // ---------------------------------------------------
  // LINEAR MATTER POWER SPECTRUM
  // size = (lnP_nk+1,lnP_nz+1)
  // z = lnPL[lnP_nk,j<lnP_nz]
  // k = lnPL[i<lnP_nk,lnP_nz]
  // ---------------------------------------------------
  // ---------------------------------------------------
  int lnPL_nk;
  int lnPL_nz;
  double** lnPL;
  // ---------------------------------------------------
  // ---------------------------------------------------
  // DISTANCE
  // ---------------------------------------------------
  // ---------------------------------------------------
  // z   = chi[0,j<chi_nz]
  // chi = chi[1,j<chi_nz]
  int chi_nz;
  double** chi;
  // ---------------------------------------------------
  // ---------------------------------------------------
  // GROWTH FACTOR
  // ---------------------------------------------------
  // ---------------------------------------------------
  // z   = G[0,j<chi_nz]
  // chi = G[1,j<chi_nz]
  int G_nz;
  double** G;
} cosmopara;

typedef struct
{ // parameters for power spectrum passed to FASTPT
  int N;
  double k_min;
  double k_max;
  double** tab; 
} FPT;

typedef struct
{
  // ---------------------------------------------------
  // ---------------------------------------------------
  // CACHE VARIABLES
  // ---------------------------------------------------
  // ---------------------------------------------------
  double random_photoz_shear;
  double random_photoz_clustering;
  double random_ia;
  double random_galaxy_bias;
  double random_gas;
  // ---------------------------------------------------
  // ---------------------------------------------------
  // INTRINSIC ALIGMENT --------------------------------
  // ---------------------------------------------------
  // --------------------------------------------------- 
  // ia[0][0] = A_ia          if(IA_NLA_LF || IA_REDSHIFT_EVOLUTION)
  // ia[0][1] = eta_ia        if(IA_NLA_LF || IA_REDSHIFT_EVOLUTION)
  // ia[0][2] = eta_ia_highz  if(IA_NLA_LF, Joachimi2012)
  // ia[0][3] = beta_ia       if(IA_NLA_LF, Joachimi2012)
  // ia[0][4] = LF_alpha      if(IA_NLA_LF, Joachimi2012)
  // ia[0][5] = LF_P          if(IA_NLA_LF, Joachimi2012)
  // ia[0][6] = LF_Q          if(IA_NLA_LF, Joachimi2012)
  // ia[0][7] = LF_red_alpha  if(IA_NLA_LF, Joachimi2012)
  // ia[0][8] = LF_red_P      if(IA_NLA_LF, Joachimi2012)
  // ia[0][9] = LF_red_Q      if(IA_NLA_LF, Joachimi2012)
  // ------------------
  // ia[1][0] = A2_ia        if IA_REDSHIFT_EVOLUTION
  // ia[1][1] = eta_ia_tt    if IA_REDSHIFT_EVOLUTION
  // ------------------
  // ia[2][MAX_SIZE_ARRAYS] = b_ta_z[MAX_SIZE_ARRAYS]
  
  int IA;
  int IA_MODEL;
  double ia[MAX_SIZE_ARRAYS][MAX_SIZE_ARRAYS];
  double oneplusz0_ia;
  double c1rhocrit_ia;
  // ---------------------------------------------------
  // ---------------------------------------------------
  // PHOTOZ --------------------------------------------
  // ---------------------------------------------------
  // ---------------------------------------------------
  // 1st index: photoz[0][:][:] = SHEAR; photoz[1][:][:] = CLUSTERING
  // 2nd index: photoz[:][0][:] = bias; photoz[:][1][:] = stretch
  double photoz[MAX_SIZE_ARRAYS][MAX_SIZE_ARRAYS][MAX_SIZE_ARRAYS]; 
  // ---------------------------------------------------
  // ---------------------------------------------------
  // SHEAR CALIBRATION ---------------------------------
  // ---------------------------------------------------
  // ---------------------------------------------------
  double shear_calibration_m[MAX_SIZE_ARRAYS];
  // ---------------------------------------------------
  // ---------------------------------------------------
  // GALAXY BIAS ---------------------------------------
  // ---------------------------------------------------
  // ---------------------------------------------------
  // 1st index: b[0][i]: linear galaxy bias in clustering bin i
  //            b[1][i]: nonlinear b2 galaxy bias in clustering bin i
  //            b[2][i]: leading order tidal bs2 galaxy bias in clustering bin i
  //            b[3][i]: nonlinear b3 galaxy bias in clustering bin i 
  //            b[4][i]: amplitude of magnification bias in clustering bin i 
  double gb[MAX_SIZE_ARRAYS][MAX_SIZE_ARRAYS]; // galaxy bias
  // HOD[i] contains HOD parameters of galaxies in clustering bin i
  // 5 parameter model of Zehavi et al. 2011 + modification of concentration
  double hod[MAX_SIZE_ARRAYS][MAX_SIZE_ARRAYS]; 
  double gc[MAX_SIZE_ARRAYS];  // galaxy concentration parameter
  // ---------------------------------------------------
  // ---------------------------------------------------
  // GAS -----------------------------------------------
  // ---------------------------------------------------
  // ---------------------------------------------------
  //gas[0] = gas_Gamma_KS; // Gamma in K-S profile
  //gas[1] = gas_beta;     // beta: mass scaling index in bound gas fraction
  //gas[2] = gas_lgM0;     // critical halo mass, below which gas ejection is significant
  //gas[3] = gas_eps1;
  //gas[4] = gas_eps2;
  //gas[5] = gas_alpha;
  //gas[6] = gas_A_star;
  //gas[7] = gas_lgM_star;
  //gas[8] = gas_sigma_star;
  //gas[9] = gas_lgT_w;
  //gas[10] = gas_f_H;
  double gas[MAX_SIZE_ARRAYS]; // Compton-Y related variables
  // ---------------------------------------------------
  // ---------------------------------------------------
  // CLUSTER ROUTINES (ALPHA STAGE)
  // ---------------------------------------------------
  // ---------------------------------------------------
  /*
  // Variables for the 4x2pt+N (see: 2008.10757 & 2010.01138)
  double cluster_MOR[MAX_SIZE_ARRAYS];
  double cluster_selection[MAX_SIZE_ARRAYS];
  */
} nuisanceparams;

typedef struct
{
  int is_cmb_bandpower;
  int is_cmb_kkkk_cov_from_sim;// if kkkk covmat is from sim, apply Hartlap factor
  double alpha_Hartlap_kkkk;
  int Ncl;
  int Ncos;
  int Ndata;
  int Nbp; // number of band-power bins for CMB lensing
  int lmin;
  int lmax;
  double* ell;
  double cosmax;
  double Rmin_bias;
  double Rmin_shear;
  int lmax_shear;
  int lmin_bp;      // minimum l mode of CMB band power binning matrix (w/ corr)
  int lmax_bp;      // Note: for dv, only need binmat with corr
  int lmin_kappacmb;
  int lmax_kappacmb;
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
  int adopt_limber_gg;
  int use_ggl_efficiency_zoverlap;
  // ---------------------------------------------------
  // ---------------------------------------------------
  // HALO MODEL CHOICES
  // ---------------------------------------------------
  // ---------------------------------------------------
  int galaxy_bias_model[MAX_SIZE_ARRAYS]; // [0] = b1, 
                                          // [1] = b2, 
                                          // [2] = bs2, 
                                          // [3] = b3, 
                                          // [4] = bmag 
  int halo_model[MAX_SIZE_ARRAYS]; // [0] = HMF, 
                                   // [1] = BIAS, 
                                   // [2] = CONCENTRATION
                                   // [3] = HALO PROFILE
} likepara;

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
  char runmode[CHAR_MAX_SIZE];
} pdeltapara;

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
  int shear_Npowerspectra;       // number of shear-shear tomography power spectra
  int ggl_Npowerspectra;         // number of galaxy-galaxy lensing tomography power spectra
  int clustering_Npowerspectra;  // number of galaxy-galaxy clustering tomography power spectra
  int* ggl_exclude;              // l-s pairs that are excluded in ggl
  int N_ggl_exclude;             // number of l-s ggl pairs excluded
} tomopara;

typedef struct
{
  // ---------------------------------------------------
  // ---------------------------------------------------
  // CACHE VARIABLES
  // ---------------------------------------------------
  // ---------------------------------------------------
  double random_shear;
  double random_clustering;
  // ---------------------------------------------------
  // ---------------------------------------------------
  // SOURCE n(Z)
  // ---------------------------------------------------
  // ---------------------------------------------------
  int shear_nbin;         // number of source tomography bins
  int shear_photoz;
  int shear_nzbins;
  double** shear_zdist_table;
  double shear_zdist_zmin_all;
  double shear_zdist_zmax_all;
  double shear_zdist_zmin[MAX_SIZE_ARRAYS];
  double shear_zdist_zmax[MAX_SIZE_ARRAYS];
  // ---------------------------------------------------
  // ---------------------------------------------------
  // CLUSTERING n(z)
  // ---------------------------------------------------
  // ---------------------------------------------------
  int clustering_nbin;    // number of lens galaxy bins
  int clustering_photoz;
  int clustering_nzbins;
  double** clustering_zdist_table;
  double clustering_zdist_zmin_all;
  double clustering_zdist_zmax_all;
  double clustering_zdist_zmin[MAX_SIZE_ARRAYS];
  double clustering_zdist_zmax[MAX_SIZE_ARRAYS];
  double clustering_zdist_zmean[MAX_SIZE_ARRAYS];
  // ---------------------------------------------------
  // ---------------------------------------------------
  // CLUSTER ROUTINES (ALPHA STAGE)
  // ---------------------------------------------------
  // ---------------------------------------------------
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
} redshiftparams;

// ---------------------------------------------------
// ---------------------------------------------------
// PGM and PGG from HALO MODEL - USE W/ DARK EMULATOR
// ---------------------------------------------------
// ---------------------------------------------------
typedef struct
{
  // ---------------------------------------------------
  // ---------------------------------------------------
  // CACHE VARIABLES
  // ---------------------------------------------------
  // ---------------------------------------------------
  double random;

  // ---------------------------------------------------
  // ---------------------------------------------------
  // PGM POWER SPECTRUM (matter-galaxy)
  // size = (lnP_nk+1,lnP_nz+1)
  // z = lnP[lnP_nk,j<lnP_nz]
  // k = lnP[i<lnP_nk,lnP_nz]
  // ---------------------------------------------------
  // ---------------------------------------------------
  int lnPGM_nk;
  int lnPGM_nz;
  double** lnPGM;
  
  // ---------------------------------------------------
  // ---------------------------------------------------
  // PGM POWER SPECTRUM (matter-galaxy)
  // size = (lnP_nk+1,lnP_nz+1)
  // z = lnP[lnP_nk,j<lnP_nz]
  // k = lnP[i<lnP_nk,lnP_nz]
  // ---------------------------------------------------
  // ---------------------------------------------------
  int lnPGG_nk;
  int lnPGG_nz;
  double** lnPGG;
} halomodel;

/*
// ---------------------------------------------------
// CLUSTER ROUTINES (ALPHA STAGE)
// ---------------------------------------------------
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
*/

// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------
//  EXTERN STRUCTS 
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------

extern likepara like;

extern cosmopara cosmology;

extern tomopara tomo;

extern redshiftparams redshift;

extern sur survey;

extern pdeltapara pdeltaparams;

extern nuisanceparams nuisance;

extern Cmb cmb;

extern lim limits;

extern Ntab Ntable;

extern FPT FPTbias;

extern FPT FPTIA;

extern halomodel hmemu; // halo model emulator

//extern clusterparams Cluster;

// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------
//  RESET STRUCTS 
// --------------------------------------------------------------------
// --------------------------------------------------------------------
// --------------------------------------------------------------------

void reset_like_struct();
void reset_bary_struct();
void reset_nuisance_struct();
void reset_redshift_struct();
void reset_cosmology_struct();
void reset_tomo_struct();
void reset_Ntable_struct();
void reset_halomodel_struct();

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD
