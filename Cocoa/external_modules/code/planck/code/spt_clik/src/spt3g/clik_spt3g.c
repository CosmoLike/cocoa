#include "clik.h"
#include "clik_helper.h"
#include <errno.h>
#include <string.h>


void spt3g_free_();
void spt3g_lkl_(double*, double*);
void spt3g_only_one_(int*);
void spt3g_parameter_init_(int*,int*,
                           double*,
                           int*,int*,int*,
                           double*,double*,double*,
                           //int*, char*,
                           int*,int*,int*,
                           int*,double*,double*,
                           int*,double*,double*,
                           int*,double*,double*,double*,
                           double*,
                           int*,int*,
                           double*, double*,
                           int*,int*,
                           int*, double*,double*,double*,
                           int*,int*);





typedef struct {
  char tmpdir[800];
  } spt3g_teee_2020;

void free_spt3g(void **none) {
  //spt3g_free_();
}

double spt3g_lkl(void* none, double* pars, error **err) {
  double lkl;
  
  spt3g_lkl_(&lkl,pars);
  return -lkl;
}


cmblkl* clik_spt3g_init(cldf *df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  char *xnames_def[] = {"Kappa", "Dl_Poisson_90x90", "Dl_Poisson_90x150", 
                        "Dl_Poisson_90x220", "Dl_Poisson_150x150", 
                        "Dl_Poisson_150x220", "Dl_Poisson_220x220", 
                        "TDust", "ADust_TE_150", "BetaDust_TE", "AlphaDust_TE", 
                        "ADust_EE_150", "BetaDust_EE", "AlphaDust_EE", "MapTcal90", 
                        "MapTcal150", "MapTcal220", "MapPcal90", "MapPcal150", 
                        "MapPcal220"};

  int SPT3G_windows_lmin,SPT3G_windows_lmax;
  double aberration_coefficient,nu_0_radio, nu_0_dsfg, nu_0_dust,Kappa_prior_mean, Kappa_prior_sigma;
  double  AlphaDustEE_prior_mean,AlphaDustEE_prior_sigma, AlphaDustTE_prior_mean, AlphaDustTE_prior_sigma,BeamCovScaling;
  int ssl_switch, poisson_switch, dust_switch,ibin_min,ibin_max,Cal_prior_switch,Kappa_prior_switch;
  int AlphaDustEE_prior_switch, AlphaDustTE_prior_switch;
  int sz_spectra_to_fit_str;
  char* spectra_to_fit_str;
  double *full_bandpowers;
  int sz_full_bandpowers;
  int* spectra_to_fit_bandpower_indices;
  int N_s, N_b;
  double *bdp_covariance, *beam_covariance,*spectra_to_fit_nu_eff,*cal_inv_covariance,*full_windows;
  int *spectra_to_fit_tcal_indices, *spectra_to_fit_ecal_indices,*cal_row_use;
  int *spectra_to_fit_poisson_indices, *spectra_to_fit_kind;
  int bok;
  int xdim;
  cmblkl *cing;

  // make sure we have only one instance of the likelihood
  spt3g_only_one_(&bok);
  testErrorRet(bok!=0,-100,"spt3g already initialized",*err,__LINE__,NULL);

  SPT3G_windows_lmin = cldf_readint(df,"SPT3G_Y1_EETE_window_l_min",err);
  forwardError(*err,__LINE__,NULL);
  SPT3G_windows_lmax = cldf_readint(df,"SPT3G_Y1_EETE_window_l_max",err);
  forwardError(*err,__LINE__,NULL);
  
  aberration_coefficient = cldf_readfloat_default(df,"SPT3G_Y1_EETE_aberration_coefficient",0.0,err);
  forwardError(*err,__LINE__,NULL);

  ssl_switch = cldf_readint_default(df,"SPT3G_Y1_EETE_super_sample_lensing",1,err);
  forwardError(*err,__LINE__,NULL);
  poisson_switch = cldf_readint_default(df,"SPT3G_Y1_EETE_poisson_sources",1,err);
  forwardError(*err,__LINE__,NULL);
  dust_switch = cldf_readint_default(df,"SPT3G_Y1_EETE_dust",1,err);
  forwardError(*err,__LINE__,NULL);

  nu_0_radio = cldf_readfloat_default(df,"SPT3G_Y1_EETE_radio_galaxies_nu0",150.0,err);
  forwardError(*err,__LINE__,NULL);
  nu_0_dsfg = cldf_readfloat_default(df,"SPT3G_Y1_EETE_dsfg_nu0",150.0,err);
  forwardError(*err,__LINE__,NULL);
  nu_0_dust = cldf_readfloat_default(df,"SPT3G_Y1_EETE_dust_nu0",150.0,err);
  forwardError(*err,__LINE__,NULL);

  
  ibin_min = cldf_readint_default(df,"SPT3G_Y1_EETE_bin_min",1,err);
  forwardError(*err,__LINE__,NULL);
  testErrorRetVA(ibin_min<1 || ibin_min>44,-55555,"SPT3G_Y1_EETE_bin_min must be between 1 and 44, got %d",*err,__LINE__,NULL,ibin_min);
 
  ibin_max = cldf_readint_default(df,"SPT3G_Y1_EETE_bin_max",44,err);
  forwardError(*err,__LINE__,NULL);
  testErrorRetVA(ibin_max<ibin_min || ibin_max>44,-55555,"SPT3G_Y1_EETE_bin_max must be between 1 and 44, got %d",*err,__LINE__,NULL,ibin_min);
  
  N_b = ibin_max - ibin_min+1;

  Cal_prior_switch = cldf_readint_default(df,"SPT3G_Y1_EETE_cal_prior",0,err);
  forwardError(*err,__LINE__,NULL);

  Kappa_prior_switch = cldf_readint_default(df,"SPT3G_Y1_EETE_cal_prior",0,err);
  forwardError(*err,__LINE__,NULL);
  Kappa_prior_mean = cldf_readfloat_default(df,"SPT3G_Y1_EETE_Kappa_mean",0.0,err);
  forwardError(*err,__LINE__,NULL);
  Kappa_prior_sigma = cldf_readfloat_default(df,"SPT3G_Y1_EETE_Kappa_sigma",0.00045,err);
  forwardError(*err,__LINE__,NULL);

  AlphaDustEE_prior_switch =cldf_readint_default(df,"SPT3G_Y1_EETE_Alpha_Dust_EE_prior",0,err);
  forwardError(*err,__LINE__,NULL);
  AlphaDustEE_prior_mean = cldf_readfloat_default(df,"SPT3G_Y1_EETE_Alpha_Dust_EE_prior_mean",0.0,err);
  forwardError(*err,__LINE__,NULL);
  AlphaDustEE_prior_sigma = cldf_readfloat_default(df,"SPT3G_Y1_EETE_Alpha_Dust_EE_prior_sigma",0.00045,err);
  forwardError(*err,__LINE__,NULL);
  AlphaDustTE_prior_switch = cldf_readint_default(df,"SPT3G_Y1_EETE_Alpha_Dust_TE_prior",0,err);
  forwardError(*err,__LINE__,NULL);
  AlphaDustTE_prior_mean = cldf_readfloat_default(df,"SPT3G_Y1_EETE_Alpha_Dust_TE_prior_mean",0.0,err);
  forwardError(*err,__LINE__,NULL);
  AlphaDustTE_prior_sigma = cldf_readfloat_default(df,"SPT3G_Y1_EETE_Alpha_Dust_TE_prior_sigma",0.00045,err);
  forwardError(*err,__LINE__,NULL);

  BeamCovScaling = cldf_readfloat_default(df,"SPT3G_Y1_EETE_Beam_Covariance_Scaling",1.0,err);
  forwardError(*err,__LINE__,NULL);


  sz_full_bandpowers = 12*44;
  full_bandpowers = cldf_readfloatarray(df,"SPT3G_Y1_EETE_bandpower",&sz_full_bandpowers,err);
  forwardError(*err,__LINE__,NULL);
  
  N_s = -1;

  spectra_to_fit_bandpower_indices = cldf_readintarray(df,"SPT3G_Y1_EETE_spectra_to_fit",&N_s,err);
  forwardError(*err,__LINE__,NULL);
  spectra_to_fit_kind = cldf_readintarray(df,"SPT3G_Y1_EETE_spectra_to_fit_kind",&N_s,err);
  forwardError(*err,__LINE__,NULL);

  int covsz;
  covsz = N_s*N_b*N_s*N_b;

  bdp_covariance = cldf_readfloatarray(df,"SPT3G_Y1_EETE_covariance_matrix",&covsz,err);
  forwardError(*err,__LINE__,NULL);

  beam_covariance = cldf_readfloatarray(df,"SPT3G_Y1_EETE_beam_covariance_matrix",&covsz,err);
  forwardError(*err,__LINE__,NULL);

  int szcali;
  szcali = N_s*2;
  spectra_to_fit_tcal_indices = cldf_readintarray(df,"SPT3G_Y1_EETE_spectra_to_fit_tcal",&szcali,err);
  forwardError(*err,__LINE__,NULL);
  spectra_to_fit_ecal_indices = cldf_readintarray(df,"SPT3G_Y1_EETE_spectra_to_fit_ecal",&szcali,err);
  forwardError(*err,__LINE__,NULL);
  spectra_to_fit_nu_eff  = cldf_readfloatarray(df,"SPT3G_Y1_EETE_central_frequency",&szcali,err);
  forwardError(*err,__LINE__,NULL);
  
  szcali = 6;
  cal_row_use = cldf_readintarray(df,"SPT3G_Y1_EETE_cal_row_use",&szcali,err);
  forwardError(*err,__LINE__,NULL);
  
  szcali = cal_row_use[0] + cal_row_use[1] + cal_row_use[2] + cal_row_use[3] + cal_row_use[4] + cal_row_use[5];
  szcali = szcali*szcali;
  cal_inv_covariance = cldf_readfloatarray(df,"SPT3G_Y1_EETE_cal_inv_covariance_matrix",&szcali,err);
  forwardError(*err,__LINE__,NULL);

  int sz_windows;
  sz_windows = 44 * 12 * (1+SPT3G_windows_lmax-SPT3G_windows_lmin);
  full_windows =  cldf_readfloatarray(df,"SPT3G_Y1_EETE_windows",&sz_windows,err);
  forwardError(*err,__LINE__,NULL);

  spectra_to_fit_poisson_indices = cldf_readintarray(df,"SPT3G_Y1_EETE_spectra_to_fit_poisson_indices",&N_s,err);
  forwardError(*err,__LINE__,NULL);
  
  spt3g_parameter_init_(&SPT3G_windows_lmin,&SPT3G_windows_lmax,
                        &aberration_coefficient,
                        &ssl_switch, &poisson_switch, &dust_switch,
                        &nu_0_radio, &nu_0_dsfg, &nu_0_dust,
                        &ibin_min,&ibin_max,&Cal_prior_switch,
                        &Kappa_prior_switch,&Kappa_prior_mean, &Kappa_prior_sigma,
                        &AlphaDustEE_prior_switch,&AlphaDustEE_prior_mean, &AlphaDustEE_prior_sigma,
                        &AlphaDustTE_prior_switch,&AlphaDustTE_prior_mean, &AlphaDustTE_prior_sigma,&BeamCovScaling,
                        full_bandpowers, &N_s, spectra_to_fit_bandpower_indices,
                        bdp_covariance,beam_covariance,
                        spectra_to_fit_tcal_indices,spectra_to_fit_ecal_indices,
                        cal_row_use,cal_inv_covariance,full_windows,spectra_to_fit_nu_eff,spectra_to_fit_poisson_indices,spectra_to_fit_kind);

  free(full_bandpowers);
  free(spectra_to_fit_bandpower_indices);
  free(bdp_covariance);
  free(beam_covariance);
  free(spectra_to_fit_tcal_indices);
  free(spectra_to_fit_ecal_indices);
  free(cal_row_use);
  free(cal_inv_covariance);
  free(full_windows);
  free(spectra_to_fit_nu_eff);
  free(spectra_to_fit_kind);
  free(spectra_to_fit_poisson_indices);

  xdim = 20;
  
  cing = init_cmblkl(NULL, &spt3g_lkl, 
                     &free_spt3g,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,xdim,err);
  forwardError(*err,__LINE__,NULL);
  
  cmblkl_set_names(cing, xnames_def,err);
  forwardError(*err,__LINE__,NULL);
  
  return cing;
}


// spt3g ttteee 2018 yeah
void spt3g_ttteee2018_parameter_init_(int*, int*, double*, char*, int* , double* , 
                                          char*, int* , double*,  char*, int*  ,   
                                          double*, double*, char*, int* ,     
                                          double*,char*, int* , char*, int* ,      
                                          char*, int* ,  char*, int* ,char*, int* ,
                                          double*,double*,double*, 
                                          double*, double*, double*, double*, double*, double*, int*,double*,int*,double*,int*);



void spt3g_ttteee2018_only_one_(int*);
void spt3g_ttteee2018_free_();
void spt3g_ttteee2018_lkl_(double*, double*);

double spt3g_ttteee_2018_lkl(void* none, double* pars, error **err) {
  double lkl;
  
  spt3g_ttteee2018_lkl_(&lkl,pars);
  return -lkl;
}

void spt3g_ttteee2018_free(void **none) {
  spt3g_ttteee2018_free_();
}


cmblkl* clik_spt3g_ttteee_2018_options_init(cldf *df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, cdic *options,error **err) {
  char *xnames_def[] = {"Kappa" ,"Tcal90" ,"Tcal150" ,"Tcal220" ,"Ecal90" ,"Ecal150" ,
                        "Ecal220" ,"EE_Poisson_90x90" ,"EE_Poisson_90x150" ,"EE_Poisson_90x220", 
                        "EE_Poisson_150x150" ,"EE_Poisson_150x220" ,"EE_Poisson_220x220" ,
                        "EE_PolGalDust_Amp" ,"EE_PolGalDust_Alpha" ,"EE_PolGalDust_Beta" ,
                        "TE_PolGalDust_Amp" ,"TE_PolGalDust_Alpha" ,"TE_PolGalDust_Beta" ,
                        "TT_Poisson_90x90" ,"TT_Poisson_90x150" ,"TT_Poisson_90x220" ,
                        "TT_Poisson_150x150" ,"TT_Poisson_150x220" ,"TT_Poisson_220x220" ,
                        "TT_GalCirrus_Amp" ,"TT_GalCirrus_Alpha" ,"TT_GalCirrus_Beta" ,
                        "TT_CIBClustering_Amp" ,"TT_CIBClustering_Alpha" ,"TT_CIBClustering_Beta" ,
                        "TT_CIBClustering_decorr_90" ,"TT_CIBClustering_decorr_150" ,
                        "TT_CIBClustering_decorr_220" ,"TT_tSZ_Amp" ,"TT_tSZ_CIB_corr" ,"TT_kSZ_Amp", 
                        "H0",  "Omega_b", "sigma_8", "Omega_m", "ns", "tau"};

  int SPT3G_windows_lmin,SPT3G_windows_lmax;
  char *spectra_to_fit_list_string, *spec_bin_min_list_string, *spec_bin_max_list_string, *late_crop_msk_string,                \
        *full_bandpower_list_string, *full_covariance_list_string, *full_window_list_string, *full_beam_covariance_list_string, \
        *nu_eff_list_string;                       
  int l_spectra_to_fit_list_string, l_spec_bin_min_list_string, l_spec_bin_max_list_string, l_late_crop_msk_string,              \
      l_full_bandpower_list_string, l_full_covariance_list_string, l_full_window_list_string, l_full_beam_covariance_list_string,\
      l_nu_eff_list_string;

  double cov_eval_cut_threshold, cov_eval_large_number_replacement, beam_cov_scale, aberration_coefficient;
  int hk;
  int N_freq_0, N_b_0_TT, N_b_0_TE, N_b_0_EE, N_s_0, N_b_0_total;

  double *full_bandpowers, *full_covariance_matrix, *full_beam_covariance_matrix, *full_cal_covariance_matrix, *full_windows, *nu_eff_matrix;
  int sz_full_bandpowers, sz_full_covariance_matrix, sz_full_beam_covariance_matrix, sz_full_cal_covariance_matrix, sz_full_windows, sz_nu_eff_matrix;

  double nu_0_galdust, T_galdust, nu_0_CIB, T_CIB, nu_0_tSZ;
  int tSZCosmologyScalingEnabled, kSZCosmologyScalingEnabled;

  double *full_tSZ_template, *full_kSZ_template;
  int sz_full_tSZ_template, sz_full_kSZ_template;

  int include_logdet;

  int bok;
  int xdim;
  int version;

  cmblkl *cing;

  // make sure we have only one instance of the likelihood
  spt3g_only_one_(&bok);
  testErrorRet(bok!=0,-100,"spt3g ttteee 2018 already initialized",*err,__LINE__,NULL);

  version = cldf_readint(df,"version",err);
  forwardError(*err,__LINE__,NULL);
  testErrorRet((version!=1) && (version!=2),-100,"spt3g ttteee 2018 version unknown",*err,__LINE__,NULL);
  

  SPT3G_windows_lmin = cldf_readint(df,"SPT3G_2018_TTTEEE_window_l_min",err);
  forwardError(*err,__LINE__,NULL);
  SPT3G_windows_lmax = cldf_readint(df,"SPT3G_2018_TTTEEE_window_l_max",err);
  forwardError(*err,__LINE__,NULL);
  
  N_freq_0 = cldf_readint_default(df,"SPT3G_2018_TTTEEE_N_freq_0",3,err);
  forwardError(*err,__LINE__,NULL);
  
  N_b_0_TT = cldf_readint_default(df,"SPT3G_2018_TTTEEE_N_b_0_TT",44,err);
  forwardError(*err,__LINE__,NULL);
  N_b_0_TE = cldf_readint_default(df,"SPT3G_2018_TTTEEE_N_b_0_TE",44,err);
  forwardError(*err,__LINE__,NULL);
  N_b_0_EE = cldf_readint_default(df,"SPT3G_2018_TTTEEE_N_b_0_EE",44,err);
  forwardError(*err,__LINE__,NULL);
  N_b_0_total = cldf_readint_default(df,"SPT3G_2018_TTTEEE_N_b_0_total",6*N_b_0_TT+6*N_b_0_TE+6*N_b_0_EE,err);
  forwardError(*err,__LINE__,NULL);
  
  N_s_0 = cldf_readint_default(df,"SPT3G_2018_TTTEEE_N_s_0",18,err);
  forwardError(*err,__LINE__,NULL);

  l_spectra_to_fit_list_string = -1;
  spectra_to_fit_list_string = opdf_readstr(df,"SPT3G_2018_TTTEEE_spectra_to_fit",&l_spectra_to_fit_list_string,options,err);
  forwardError(*err,__LINE__,NULL);
  
  l_spec_bin_min_list_string = -1;
  spec_bin_min_list_string = opdf_readstr(df,"SPT3G_2018_TTTEEE_spectra_to_fit_bin_min",&l_spec_bin_min_list_string,options,err);
  forwardError(*err,__LINE__,NULL);
  //_DEBUGHERE_("%s",spec_bin_min_list_string);
  
  l_spec_bin_max_list_string = -1;
  spec_bin_max_list_string = opdf_readstr(df,"SPT3G_2018_TTTEEE_spectra_to_fit_bin_max",&l_spec_bin_max_list_string,options,err);
  forwardError(*err,__LINE__,NULL);
  
  hk = opdf_haskey(df,"SPT3G_2018_TTTEEE_late_crop_msk",options,err);
  forwardError(*err,__LINE__,NULL);
  if (hk==1) {
    l_late_crop_msk_string = -1;
    late_crop_msk_string = opdf_readstr(df,"SPT3G_2018_TTTEEE_late_crop_msk",&l_late_crop_msk_string,options,err);
    forwardError(*err,__LINE__,NULL);
  } else {
    late_crop_msk_string = malloc_err(sizeof(char),err);
    forwardError(*err,__LINE__,NULL);
    late_crop_msk_string[0] = ' ';
    l_late_crop_msk_string = 1;  
  }

  l_full_bandpower_list_string = -1;
  full_bandpower_list_string = cldf_readstr(df,"SPT3G_2018_TTTEEE_bandpower_file_order",&l_full_bandpower_list_string,err);
  forwardError(*err,__LINE__,NULL);
  
  l_full_covariance_list_string = -1;
  full_covariance_list_string = cldf_readstr(df,"SPT3G_2018_TTTEEE_covariance_matrix_order",&l_full_covariance_list_string,err);
  forwardError(*err,__LINE__,NULL);
  
  cov_eval_cut_threshold = cldf_readfloat_default(df,"SPT3G_2018_TTTEEE_cov_eval_cut_threshold",1e3,err);
  forwardError(*err,__LINE__,NULL);

  cov_eval_large_number_replacement = cldf_readfloat_default(df,"SPT3G_2018_TTTEEE_cov_eval_replacement",0.2,err);
  forwardError(*err,__LINE__,NULL);

  l_full_beam_covariance_list_string = -1;
  full_beam_covariance_list_string = cldf_readstr(df,"SPT3G_2018_TTTEEE_beam_covariance_matrix_order",&l_full_beam_covariance_list_string,err);
  forwardError(*err,__LINE__,NULL);

  beam_cov_scale  = opdf_readfloat_default(df,"SPT3G_2018_TTTEEE_beam_covariance_scale",1,options,err);
  forwardError(*err,__LINE__,NULL);

  l_full_window_list_string = -1;
  full_window_list_string = cldf_readstr(df,"SPT3G_2018_TTTEEE_window_folder_order",&l_full_window_list_string,err);
  forwardError(*err,__LINE__,NULL);

  l_nu_eff_list_string = -1;
  nu_eff_list_string = cldf_readstr(df,"SPT3G_2018_TTTEEE_central_frequency_file_order",&l_nu_eff_list_string,err);
  forwardError(*err,__LINE__,NULL);

  aberration_coefficient = opdf_readfloat_default(df,"SPT3G_2018_TTTEEE_aberration_coefficient",0.0,options,err);
  forwardError(*err,__LINE__,NULL);

  sz_full_bandpowers = N_b_0_total;
  full_bandpowers = cldf_readfloatarray(df,"SPT3G_2018_TTTEEE_bandpower",&sz_full_bandpowers,err);
  forwardError(*err,__LINE__,NULL);

  sz_full_covariance_matrix = N_b_0_total*N_b_0_total;
  full_covariance_matrix = cldf_readfloatarray(df,"SPT3G_2018_TTTEEE_covariance_matrix",&sz_full_covariance_matrix,err);
  forwardError(*err,__LINE__,NULL);

  sz_full_beam_covariance_matrix = N_b_0_total*N_b_0_total;
  full_beam_covariance_matrix = cldf_readfloatarray(df,"SPT3G_2018_TTTEEE_beam_covariance_matrix",&sz_full_beam_covariance_matrix,err);
  forwardError(*err,__LINE__,NULL);

  sz_full_cal_covariance_matrix = N_freq_0*2*N_freq_0*2;
  full_cal_covariance_matrix = cldf_readfloatarray(df,"SPT3G_2018_TTTEEE_cal_covariance_matrix",&sz_full_cal_covariance_matrix,err);
  forwardError(*err,__LINE__,NULL);

  sz_full_windows = N_b_0_EE * N_s_0 * (1+SPT3G_windows_lmax-SPT3G_windows_lmin);
  full_windows =  cldf_readfloatarray(df,"SPT3G_2018_TTTEEE_windows",&sz_full_windows,err);
  forwardError(*err,__LINE__,NULL);

  sz_nu_eff_matrix = 5*N_freq_0;
  nu_eff_matrix =  cldf_readfloatarray(df,"SPT3G_2018_TTTEEE_central_frequency",&sz_nu_eff_matrix,err);
  forwardError(*err,__LINE__,NULL);

  nu_0_galdust = cldf_readfloat_default(df,"SPT3G_2018_TTTEEE_galdust_nu0",150.0,err);
  forwardError(*err,__LINE__,NULL);
  
  T_galdust = opdf_readfloat_default(df,"SPT3G_2018_TTTEEE_galdust_T",19.6,options,err);
  forwardError(*err,__LINE__,NULL);

  nu_0_CIB = cldf_readfloat_default(df,"SPT3G_2018_TTTEEE_galdust_nu0",150.0,err);
  forwardError(*err,__LINE__,NULL);
  
  T_CIB = opdf_readfloat_default(df,"SPT3G_2018_TTTEEE_CIB_T",25.,options,err);
  forwardError(*err,__LINE__,NULL);

  nu_0_tSZ = cldf_readfloat_default(df,"SPT3G_2018_TTTEEE_tSZ_nu0",143.0,err);
  forwardError(*err,__LINE__,NULL);

  tSZCosmologyScalingEnabled = opdf_readint(df,"SPT3G_2018_TTTEEE_tSZ_cosmology_scaling",options,err);
  forwardError(*err,__LINE__,NULL);

  sz_full_tSZ_template = (1+SPT3G_windows_lmax-SPT3G_windows_lmin);
  full_tSZ_template = cldf_readfloatarray(df,"SPT3G_2018_TTTEEE_tSZ_template",&sz_full_tSZ_template,err);
  forwardError(*err,__LINE__,NULL);

  kSZCosmologyScalingEnabled = opdf_readint(df,"SPT3G_2018_TTTEEE_kSZ_cosmology_scaling",options,err);
  forwardError(*err,__LINE__,NULL);

  sz_full_kSZ_template = (1+SPT3G_windows_lmax-SPT3G_windows_lmin);
  full_kSZ_template = cldf_readfloatarray(df,"SPT3G_2018_TTTEEE_kSZ_template",&sz_full_kSZ_template,err);
  forwardError(*err,__LINE__,NULL);

  include_logdet = opdf_readint_default(df,"SPT3G_2018_TTTEEE_include_logdet",1,options,err);
  forwardError(*err,__LINE__,NULL);


  spt3g_ttteee2018_parameter_init_(&SPT3G_windows_lmin, &SPT3G_windows_lmax,                                                              
                                   full_bandpowers, full_bandpower_list_string, &l_full_bandpower_list_string,                            
                                   full_covariance_matrix, full_covariance_list_string, &l_full_covariance_list_string,                  
                                   full_beam_covariance_matrix,  full_beam_covariance_list_string, &l_full_beam_covariance_list_string , 
                                   full_cal_covariance_matrix, full_windows, full_window_list_string, &l_full_window_list_string,        
                                   nu_eff_matrix,nu_eff_list_string, &l_nu_eff_list_string,                                              
                                   spectra_to_fit_list_string, &l_spectra_to_fit_list_string,                                            
                                   spec_bin_min_list_string, &l_spec_bin_min_list_string,                                                
                                   spec_bin_max_list_string, &l_spec_bin_max_list_string,                                                
                                   late_crop_msk_string, &l_late_crop_msk_string,                                                        
                                   &cov_eval_cut_threshold, &cov_eval_large_number_replacement,                                           
                                   &beam_cov_scale, &aberration_coefficient,                                                              
                                   &nu_0_galdust, &T_galdust, &nu_0_CIB, &T_CIB, &nu_0_tSZ,                                                  
                                   &tSZCosmologyScalingEnabled,full_tSZ_template,                                                        
                                   &kSZCosmologyScalingEnabled,full_kSZ_template,
                                   &include_logdet);                                                               
  
  free(spectra_to_fit_list_string);
  free(spec_bin_min_list_string);
  free(spec_bin_max_list_string);
  free(late_crop_msk_string);
  free(full_bandpower_list_string);
  free(full_covariance_list_string);
  free(full_window_list_string);
  free(full_beam_covariance_list_string);
  free(nu_eff_list_string );
  free(full_bandpowers);
  free(full_covariance_matrix);
  free(full_beam_covariance_matrix);
  free(full_cal_covariance_matrix);
  free(full_windows);
  free(nu_eff_matrix);
  free(full_tSZ_template);
  free(full_kSZ_template);

  if (tSZCosmologyScalingEnabled==1 || kSZCosmologyScalingEnabled ==1) {
    xdim = 43;
  } else {
    xdim =37;
  }
  
  
  cing = init_cmblkl(NULL, &spt3g_ttteee_2018_lkl, 
                   &spt3g_ttteee2018_free,
                     nell,ell,
                     has_cl,ell[nell-1],unit,wl,0,bins,nbins,xdim,err);
  forwardError(*err,__LINE__,NULL);
  
  cmblkl_set_names(cing, xnames_def,err);
  forwardError(*err,__LINE__,NULL);
  
  return cing;
}

cmblkl* clik_spt3g_ttteee_2018_init(cldf *df, int nell, int* ell, int* has_cl, double unit,double* wl, double *bins, int nbins, error **err) {
  cdic  *options;
  cmblkl * res;

  options = cdic_init(err);
  forwardError(*err,__LINE__,NULL);
  
  res = clik_spt3g_ttteee_2018_options_init(df, nell, ell, has_cl, unit,wl, bins, nbins, options,err);
  forwardError(*err,__LINE__,NULL);

  return res;  

}
