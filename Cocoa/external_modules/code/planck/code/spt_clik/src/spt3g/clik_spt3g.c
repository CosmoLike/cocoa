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
