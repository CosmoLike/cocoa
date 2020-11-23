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
} pdeltapara;

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

double bgal_z(double z, int nz);

double b1_per_bin(double z, int nz);

extern likepara like;

extern cosmopara cosmology;

extern tomopara tomo;

extern redshiftpara redshift;

extern sur survey;

extern Cmb cmb;

extern galpara gbias;

extern clusterpara Cluster;

extern pdeltapara pdeltaparams;

extern FPTpara FPT;

extern nuisancepara nuisance;

extern barypara bary;

extern fft_optimize fft_int;

#ifdef __cplusplus
}
#endif
#endif // HEADER GUARD