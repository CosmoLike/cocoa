#ifndef __COSMOLIKE_STRUCTS_H
#define __COSMOLIKE_STRUCTS_H
#ifdef __cplusplus
extern "C" {
#endif

#define CHAR_MAX_SIZE 1024
#define MAX_SIZE_ARRAYS 10

typedef struct {
  int Ncl;
  int Ntheta;
  int Ncos;
  int Ndata;
  double lmin;
  double lmax;
  double vtmax;
  double vtmin;
  double *ell;
  double *theta;
  double cosmax;
  double Rmin_bias;
  double Rmin_shear;
  double lmax_shear;
  double lmax_kappacmb;
  int baryons;
  int IA;
  int bias;
  int wlphotoz;
  int clphotoz;
  int shearcalib;
  int clusterMobs;
  int Planck15_BAO_H070p6_JLA_w0wa;       // CH
  int Planck18_BAO_Riess18_Pantheon_w0wa; // CH
  int Planck18_BAO_w0wa;                  // CH
  int Planck18_w0;                        // CH
  int BAO;
  int SN_WFIRST;
  int GRS;
  int SRD;
  char DATA_FILE[CHAR_MAX_SIZE];
  char INV_FILE[CHAR_MAX_SIZE];
  char COV_FILE[CHAR_MAX_SIZE];
  char BARY_FILE[CHAR_MAX_SIZE];
  char MASK_FILE[CHAR_MAX_SIZE];
  int shear_shear;
  int shear_pos;
  int pos_pos;
  int clusterN;
  int clusterWL;
  int clusterCG;
  int clusterCC;
  // MANUWARNING: added "int gk, kk, ks;"
  int gk;
  int kk;
  int ks;
  char probes[CHAR_MAX_SIZE];
  char ext_data[CHAR_MAX_SIZE];
  int theta_s;
} likepara;

typedef struct {
  double Omega_m; // matter density parameter
  double Omega_v; // cosmogical constant parameter
  double h0;      // Hubble constant
  double Omega_nu; // density parameter of massive neutrinos;
                   // Omega_m = Omega_cdm+ Omega_nu + omb
  double coverH0;  // units for comoving distances - speeds up code
  double rho_crit; // = 3 H_0^2/(8 pi G), critical comoving density
  double MGSigma;
  double MGmu;
  double random; // Random number between zero and 1 - see interface.cpp
  int is_cached; 
} cosmopara;

typedef struct {
  int shear_Nbin;          // number of tomography bins
  int shear_Npowerspectra; // number of tomography power spectra+2+3+...+Nbin
  double shear_zmax[10];   // code needs modification if more than 10 zbins
  double shear_zmin[10];
  double n_source[10];
  int clustering_Nbin;          // number of tomography bins
  int clustering_Npowerspectra; // number of tomography power
                                // spectra+2+3+...+Nbin
  double clustering_zmax[10];
  double clustering_zmin[10];
  double n_lens[10];
  int cluster_Nbin; // number of cluster redshift bins
  double cluster_zmax[10];
  double cluster_zmin[10];
  int cluster_cg_Npowerspectra; // number of cluster-lensing tomography
                                // combinations
  int cgl_Npowerspectra;  // number of cluster-lensing tomography combinations
  int ggl_Npowerspectra;  // number of ggl tomography combinations
  int magnification_Nbin; // number of tomography bins
  int magnification_Npowerspectra; // number of tomography power
                                   // spectra+2+3+...+Nbin
  double magnification_zmax[10];
  double magnification_zmin[10];
} tomopara;

typedef struct {
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

  int magnification_photoz;
  double magnification_zdistrpar_zmin;
  double magnification_zdistrpar_zmax;
  int magnification_histogram_zbins;
  char magnification_REDSHIFT_FILE[CHAR_MAX_SIZE];
} redshiftpara;

typedef struct {
  double area;                   // survey_area in deg^2.
  double n_gal;                  // galaxy density per arcmin^2
  double sigma_e;                // rms inrinsic ellipticity noise
  double area_conversion_factor; // factor from deg^2 to radian^2:
                                 // 60*60*constants.arcmin*constants.arcmin
  double
      n_gal_conversion_factor; // factor from n_gal/arcmin^2 to n_gal/radian^2:
                               // 1.0/constants.arcmin/constants.arcmin
  double n_lens; // lens galaxy density per arcmin^2
  char Kcorrect_File[CHAR_MAX_SIZE];
  double m_lim;
  char name[CHAR_MAX_SIZE];
  int surveystage;
  char sourcephotoz[CHAR_MAX_SIZE];
  char lensphotoz[CHAR_MAX_SIZE];
  char galsample[CHAR_MAX_SIZE];
  double ggl_overlap_cut;
} sur;

// MANUWARNING: check this
typedef struct {
  char name[CHAR_MAX_SIZE];
  double fwhm;            // beam fwhm in rad
  double sensitivity;     // white noise level in muK*rad
  char pathLensRecNoise[CHAR_MAX_SIZE]; // path to precomputed noise on reconstructed kappa
} Cmb;

typedef double (*B1_model)(double z, int nz);
typedef struct {
  double b[10];   // linear galaxy bias paramter in clustering bin i
  double b2[10];  // quadratic bias parameter for redshift bin i
  double bs2[10]; // leading order tidal bias for redshift bin i
  double rcorr[10];
  double hod[10][6]; // HOD[i] contains HOD parameters of galaxies in clustering
                     // bin i, following 5 parameter model of Zehavi et al. 2011
                     // + modification of concentration parameter
  double cg[10];
  double n_hod[10];
  double b_mag[10]; // amplitude of magnification bias, b_mag[i] =
                    // 5*s[i]+beta[i]-2
  B1_model b1_function;
} galpara;

typedef struct {
  double hod[5];
  double cg;            // galaxy concentration
  double fc;            // central occupation
  double b_a;           // assembly bias
  int parameterization; // Zehavi, Reddick, histo
  double n0;
  char HOD_FILE[CHAR_MAX_SIZE];
} model_para_redmagic;

typedef struct {
  double N200_min;
  double N200_max;
  int N200_Nbin;
  double N_min[10];
  double N_max[10];
  int lbin;
  double l_min;
  double l_max;
  char model[CHAR_MAX_SIZE];
} clusterpara;

typedef struct {
  char runmode[CHAR_MAX_SIZE];
  char baryons[CHAR_MAX_SIZE];
  double DIFF_n; // difference fucntion describing the constant uncertainty in
                 // Pdelta for k>0.01
  double DIFF_A; // difference fucntion describing the scale dependent
                 // uncertainty in Pdelta for k>0.01
} pdeltapara;

typedef struct {
  double sglob;
  double aglob;
  double rglob;
  int reset_cov_tables;
} globalpara;

typedef struct { // parameters for power spectrum passed to FASTPT
  // general specifiers
  double k_min;
  double k_max;
  int N;
  int N_per_dec;
  char Plin_FILE[CHAR_MAX_SIZE];
  // parameters for table of bias terms
  double **tab_AB;
  int N_AB;
  // parameters for table of IA terms - note that N_IA needs to be initialized
  // below!!!
  double **tab_IA;
  int N_IA;
  char path[CHAR_MAX_SIZE];
  cosmopara C;
} FPTpara;

typedef struct {
  // like.IA = 3: NLA, per bin
  // like.IA = 4: NLA, power law
  // like.IA = 5: TATT, per bin
  // like.IA = 6: TATT, power law
  double A_z[10];  // NLA normalization per source redshift bin, for mpp analyis
                   // (activate with like.IA =3 or like.IA = 5)
  double A2_z[10]; // NLA normalization per source redshift bin, for mpp analyis
                   // (activate with like.IA = 5)
  double b_ta_z[10];   // b_ta, per bin (like.IA = 6), or use b_ta_z[0] with
                       // like.IA = 5
  double A_ia;         // A IA see Joachimi2012
  double A2_ia;        // placeholder param for quadratic,etc IA
  double beta_ia;      // beta IA see Joachimi2012
  double eta_ia;       // eta_other IA see Joachimi2012
  double eta_ia_tt;    // same as eta_ia, for TT
  double eta_ia_highz; // uncertainty in high z evolution
  double oneplusz0_ia; // oneplusz0-ia MegaZ
  double c1rhocrit_ia;
  double fred[10];
  double shear_calibration_m[10];
  double sigma_zphot_shear[10];
  double bias_zphot_shear[10];
  double sigma_zphot_clustering[10];
  double bias_zphot_clustering[10];
  double sigma_zphot_magnification[10];
  double bias_zphot_magnification[10];
  double LF_alpha;
  double LF_P;
  double LF_Q;
  double LF_red_alpha;
  double LF_red_P;
  double LF_red_Q;
  double fA_blue; // fractional IA amplitude of blue galaxies compared to red
  double cluster_Mobs_lgM0;
  double cluster_Mobs_sigma;
  double cluster_Mobs_alpha;
  double cluster_Mobs_beta;
  double cluster_Mobs_N_pivot;
  double cluster_Mobs_lgN0;
  double cluster_Mobs_sigma0;
  double cluster_Mobs_sigma_qm;
  double cluster_Mobs_sigma_qz;
  double cluster_completeness[10];
  double cluster_centering_f0;
  double cluster_centering_alpha;
  double cluster_centering_sigma;
  double cluster_centering_M_pivot;
  int N_cluster_MOR;
  double cluster_MOR[10];
  int N_cluster_selection;
  double cluster_selection[10];
} nuisancepara;

typedef struct {  // two parameters for each nuisance parameter: Center
                  // (prior.*[0]) + width of Gaussian (prior.*[1])
  double Omega_m; // matter density parameter
  double Omega_v; // cosmogical constant parameter
  double sigma_8; // power spectrum normalization
  double A_s;
  double n_spec;  // spectral index of initial power spectrum
  double alpha_s; // running of spectral index of initial power spectrum
  double w0;      // time dependent Dark energy parametrization zero order
  double wa;      // time dependent Dark energy parametrization first order
  double omb;     // Omega baryon
  double h0;      // Hubble constant
  double M_nu;
  double Omega_nu; // density parameter of massive neutrinos; Omega_m =
                   // Omega_cdm+ Omega_nu + omb
  double f_NL;
  double MGSigma;
  double MGmu;
  double A_ia[2]; // A IA see Joachimi2012
  double A2_ia[2];
  double beta_ia[2];      // beta IA see Joachimi2012
  double eta_ia[2];       // eta_other IA see Joachimi2012
  double eta_ia_highz[2]; // additional uncertainty at high z
  double LF_alpha[2];
  double LF_P[2];
  double LF_Q[2];
  double LF_red_alpha[2];
  double LF_red_P[2];
  double LF_red_Q[2];
  double shear_calibration_m[10][2];
  double sigma_zphot_shear[10][2];
  double bias_zphot_shear[10][2];
  double sigma_zphot_clustering[10][2];
  double bias_zphot_clustering[10][2];
  double sigma_zphot_magnification[10][2];
  double bias_zphot_magnification[10][2];
  double cluster_Mobs_lgM0[2];
  double cluster_Mobs_sigma[2];
  double cluster_Mobs_alpha[2];
  double cluster_Mobs_beta[2];
  double cluster_Mobs_N_pivot[2];
  double cluster_Mobs_lgN0[2];
  double cluster_Mobs_sigma0[2];
  double cluster_Mobs_sigma_qm[2];
  double cluster_Mobs_sigma_qz[2];
  double cluster_completeness[2];
  double cluster_centering_f0[2];
  double cluster_centering_alpha[2];
  double cluster_centering_sigma[2];
  double cluster_centering_M_pivot[2];
  double bary_Q1[2];
  double bary_Q2[2];
  double bary_Q3[2];
} priorpara;

typedef struct {
  double HOD_rm[5][2];
  double cg_rm[2];
  double fc_rm[2];
} flat_priorpara;

typedef struct input_cosmo_params_mpp {
  double omega_m;
  double sigma_8;
  double A_s;
  double n_s;
  double w0;
  double wa;
  double omega_b;
  double omega_nuh2;
  double h0;
  double MGSigma;
  double MGmu;
} input_cosmo_params_mpp;

typedef struct input_cosmo_params {
  double omega_m;
  double sigma_8;
  double n_s;
  double w0;
  double wa;
  double omega_b;
  double h0;
  double MGSigma;
  double MGmu;
} input_cosmo_params;

typedef struct input_nuisance_params_mpp {
  double bias[10];
  double bias2[10];
  double lens_z_bias[10];
  double source_z_bias[10];
  double shear_m[10];
  double A_z[10];
  double MOR[10];
  double selection[10];
  double b_mag[10];
} input_nuisance_params_mpp;

typedef struct input_nuisance_params {
  double bias[10];
  double source_z_bias[10];
  double source_z_s;
  double lens_z_bias[10];
  double lens_z_s;
  double shear_m[10];
  double A_ia;
  double beta_ia;
  double eta_ia;
  double eta_ia_highz;
  double lf[6];
  double m_lambda[6];
  double cluster_c[4];
  double bary[3];
  double b_mag[10];
} input_nuisance_params;

typedef struct input_nuisance_params_grs {
  double grsbias[7];
  double grssigmap[7];
  double grssigmaz;
  double grspshot;
  double grskstar;
} input_nuisance_params_grs;

typedef struct {
  double tmin;                // Theta min (arcmin) 
  double tmax;                // Theta max (arcmin) 
  int ntheta;                 // number of theta bins 
  double lmin;                // ell min  
  double lmax;                // ell max  
  int ncl;                    // number of ell bins 
  int ng;                     // ng covariance? 
  char outdir[CHAR_MAX_SIZE];           // output directory 
  char filename[CHAR_MAX_SIZE];         // output file name prefix 
  char C_FOOTPRINT_FILE[CHAR_MAX_SIZE]; // angular power spectrum of survey footprint, in
                              //   healpix format 
  char ss[8];                 // Calculate shear-shear components 
  char ls[8];                 // Calculate shear-position components 
  char ll[8];                 // Calculate position-position components 
} covpar;

typedef struct {
  int N_z;
  int N_t;
  int N_k;
  int N_mu;
  double k_star; // in h/Mpc
  double k_min;  // in h/Mpc
  double k_max;  // in h/Mpc
  double f_sky;
  double z[10];
  double V_z[10]; // in (Mpc/h)^3
  double H_ref[10];
  double DA_ref[10];
  double **datav;
  double **var;
  double *k;
  double *mu;
  double n_mt[10][2]; // in (h/Mpc)^3
  double b_mt[10][2];
  double sigma_z; // fractional accuracy
  double sigma_p; // in km/s
} GRSpara_mt;

typedef struct {
  char FILE_logPkR[CHAR_MAX_SIZE];
  char scenario[100]; // available options: mb2, illustris, eagle, HzAGN,
                      // TNG100, owls_AGN, owls_DBLIMFV1618, owls_NOSN,
                      // owls_NOSN_NOZCOOL, owls_NOZCOOL, owls_REF, owls_WDENS,
                      // owls_WML1V848, owls_WML4
  int Nabins;
  int Nkbins;
  double z_bins[50];
  int isPkbary; // if isPkbary=1
} barypara;

typedef struct {
  double ***S_integrands_cl; // galaxy density nonlimber integrand
  double ***S_integrands_sh; // galaxy shape nonlimber integrand
  int *recompute_cl;         // recompute the above or not
  int *recompute_sh;
  int Nell;
  int Nchi;
} fft_optimize;

typedef struct input_HOD_params {
  double lgMmin[10];
  double sigma_lgMin;
  double lgM1;
  double lgM0;
  double alpha;
  double f_c;
  double c_g;
} input_HOD_params;

double bgal_z(double z, int nz);

double b1_per_bin(double z, int nz);

extern likepara like;

extern cosmopara cosmology;

extern tomopara tomo;

extern redshiftpara redshift;

extern sur survey;

extern Cmb cmb;

extern galpara gbias;

extern model_para_redmagic redm;

extern clusterpara Cluster;

extern pdeltapara pdeltaparams;

extern globalpara global;

extern FPTpara FPT;

extern nuisancepara nuisance;

extern priorpara prior;

extern flat_priorpara flat_prior;

extern covpar covparams;

extern barypara bary;

extern fft_optimize fft_int;

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD