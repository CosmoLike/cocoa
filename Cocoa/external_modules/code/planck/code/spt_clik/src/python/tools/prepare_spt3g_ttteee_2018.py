#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import numpy as nm
import clik.parobject as php
import clik
import os.path as osp
import clik.hpy as hpy

overridekeys = [
  "SPT3G_Y1_EETE_spectra_to_fit",
  "SPT3G_Y1_EETE_bin_min",
  "SPT3G_Y1_EETE_bin_max",
  "SPT3G_Y1_EETE_super_sample_lensing",
  "SPT3G_Y1_EETE_radio_galaxies",
  "SPT3G_Y1_EETE_dsfg",
  "SPT3G_Y1_EETE_dust",
  "SPT3G_Y1_EETE_radio_galaxies_nu0",
  "SPT3G_Y1_EETE_dsfg_nu0",
  "SPT3G_Y1_EETE_dust_nu0",
  "SPT3G_Y1_EETE_params_file",
  "SPT3G_Y1_EETE_bandpower_file",
  "SPT3G_Y1_EETE_bandpower_file_order",
  "SPT3G_Y1_EETE_covariance_matrix",
  "SPT3G_Y1_EETE_covariance_matrix_order",
  "SPT3G_Y1_EETE_beam_covariance_matrix",
  "SPT3G_Y1_EETE_beam_covariance_matrix_order",
  "SPT3G_Y1_EETE_window_folder",
  "SPT3G_Y1_EETE_window_folder_order",
  "SPT3G_Y1_EETE_window_l_min",
  "SPT3G_Y1_EETE_window_l_max",
  "SPT3G_Y1_EETE_central_frequency_file",
  "SPT3G_Y1_EETE_central_frequency_file_order",
  "SPT3G_Y1_EETE_cal_covariance_matrix",
  "SPT3G_Y1_EETE_cal_prior",
  "SPT3G_Y1_EETE_kappa_prior",
  "SPT3G_Y1_EETE_kappa_mean",
  "SPT3G_Y1_EETE_kappa_sigma",
  "SPT3G_Y1_EETE_Alpha_Dust_EE_prior",
  "SPT3G_Y1_EETE_Alpha_Dust_EE_mean",
  "SPT3G_Y1_EETE_Alpha_Dust_EE_sigma",
  "SPT3G_Y1_EETE_Alpha_Dust_TE_prior",
  "SPT3G_Y1_EETE_Alpha_Dust_TE_mean",
  "SPT3G_Y1_EETE_Alpha_Dust_TE_sigma",
  "SPT3G_Y1_EETE_Beam_Covariance_Scaling",
  ]

def add_data_path(data_path,path):
  if osp.exists(path):
    return path
  else:
    return osp.join(data_path,path)

option_list = [
  "SPT3G_2018_TTTEEE_late_crop_msk",
  "SPT3G_2018_TTTEEE_beam_covariance_scale",
  "SPT3G_2018_TTTEEE_galdust_T",
  "SPT3G_2018_TTTEEE_CIB_T",
  "SPT3G_2018_TTTEEE_tSZ_cosmology_scaling",
  "SPT3G_2018_TTTEEE_kSZ_cosmology_scaling",
  "SPT3G_2018_TTTEEE_spectra_to_fit_bin_min",
  "SPT3G_2018_TTTEEE_spectra_to_fit_bin_max"
]

def main(argv):
  pars = clik.miniparse(argv[1])
  
  root_grp,hf = php.baseCreateParobject(pars.res_object)
  hascl = [0]*6
  hascl[0] = 1
  hascl[1] = 1
  hascl[3] = 1
  hascl = nm.array(hascl,dtype=nm.int)
  data_path = pars.data_path
  override = {}
  for k in overridekeys:
    if k in pars:
      override[k]=getattr(pars,k)
  
  print(override)

  dataset = clik.miniparse(pars.dataset_path,**override)

  lmin = dataset.int.SPT3G_2018_TTTEEE_window_l_min
  lmax = dataset.int.SPT3G_2018_TTTEEE_window_l_max

  lkl_grp = php.add_lkl_generic(root_grp,"spt3g_ttteee_2018",1,hascl,lmax,0)
  lkl_grp.attrs["version"]=2

  php.add_option_list(lkl_grp,option_list)

  N_freq_0 = pars.int(default=3).SPT3G_2018_TTTEEE_N_freq_0
  lkl_grp["SPT3G_2018_TTTEEE_N_freq_0"] = N_freq_0
  
  N_b_0_TT = pars.int(default=44).SPT3G_2018_TTTEEE_N_b_0_TT
  lkl_grp["SPT3G_2018_TTTEEE_N_b_0_TT"] = N_b_0_TT
  
  N_b_0_TE = pars.int(default=44).SPT3G_2018_TTTEEE_N_b_0_TE
  lkl_grp["SPT3G_2018_TTTEEE_N_b_0_TE"] = N_b_0_TE
  
  N_b_0_EE = pars.int(default=44).SPT3G_2018_TTTEEE_N_b_0_EE
  lkl_grp["SPT3G_2018_TTTEEE_N_b_0_EE"] = N_b_0_EE
  
  N_b_0_total = pars.int(default=6*N_b_0_TT+6*N_b_0_TE+6*N_b_0_EE).SPT3G_2018_TTTEEE_N_b_0_total
  lkl_grp["SPT3G_2018_TTTEEE_N_b_0_total"] = N_b_0_total
  
  N_s_0 = pars.int(default=18).SPT3G_2018_TTTEEE_N_s_0
  lkl_grp["SPT3G_2018_TTTEEE_N_s_0"] = N_s_0
  
  lkl_grp["SPT3G_2018_TTTEEE_window_l_min"] = dataset.int.SPT3G_2018_TTTEEE_window_l_min
  lkl_grp["SPT3G_2018_TTTEEE_window_l_max"] = dataset.int.SPT3G_2018_TTTEEE_window_l_max

  lkl_grp["SPT3G_2018_TTTEEE_aberration_coefficient"] = pars.float(default=-0.0004826).SPT3G_2018_TTTEEE_aberration_coefficient

  lkl_grp["SPT3G_2018_TTTEEE_spectra_to_fit"] = dataset.str.SPT3G_2018_TTTEEE_spectra_to_fit

  lkl_grp["SPT3G_2018_TTTEEE_spectra_to_fit_bin_min"] = dataset.str.SPT3G_2018_TTTEEE_spectra_to_fit_bin_min
  lkl_grp["SPT3G_2018_TTTEEE_spectra_to_fit_bin_max"] = dataset.str.SPT3G_2018_TTTEEE_spectra_to_fit_bin_max

  lkl_grp["SPT3G_2018_TTTEEE_bandpower_file_order"] = dataset.str.SPT3G_2018_TTTEEE_bandpower_file_order
  lkl_grp["SPT3G_2018_TTTEEE_covariance_matrix_order"] = dataset.str.SPT3G_2018_TTTEEE_covariance_matrix_order
  lkl_grp["SPT3G_2018_TTTEEE_beam_covariance_matrix_order"] = dataset.str.SPT3G_2018_TTTEEE_beam_covariance_matrix_order
  lkl_grp["SPT3G_2018_TTTEEE_window_folder_order"] = dataset.str.SPT3G_2018_TTTEEE_window_folder_order
  lkl_grp["SPT3G_2018_TTTEEE_central_frequency_file_order"] = dataset.str.SPT3G_2018_TTTEEE_central_frequency_file_order
  lkl_grp["SPT3G_2018_TTTEEE_late_crop_msk"] = dataset.str.SPT3G_2018_TTTEEE_late_crop_msk
  lkl_grp["SPT3G_2018_TTTEEE_late_crop"] = dataset.str(default="F").SPT3G_2018_TTTEEE_late_crop.lower() in ("1","t","true",".true.")
  
  lkl_grp["SPT3G_2018_TTTEEE_cov_eval_cut_threshold"] = dataset.float.SPT3G_2018_TTTEEE_cov_eval_cut_threshold
  lkl_grp["SPT3G_2018_TTTEEE_cov_eval_replacement"] = dataset.float.SPT3G_2018_TTTEEE_cov_eval_replacement

  lkl_grp["SPT3G_2018_TTTEEE_beam_covariance_scale"] = dataset.float.SPT3G_2018_TTTEEE_beam_covariance_scale

  lkl_grp["SPT3G_2018_TTTEEE_galdust_T"] = dataset.float.SPT3G_2018_TTTEEE_galdust_T
  lkl_grp["SPT3G_2018_TTTEEE_galdust_nu0"] = dataset.float.SPT3G_2018_TTTEEE_galdust_nu0
  lkl_grp["SPT3G_2018_TTTEEE_CIB_T"] = dataset.float.SPT3G_2018_TTTEEE_CIB_T
  lkl_grp["SPT3G_2018_TTTEEE_CIB_nu0"] = dataset.float.SPT3G_2018_TTTEEE_CIB_nu0
  lkl_grp["SPT3G_2018_TTTEEE_tSZ_nu0"] = dataset.float.SPT3G_2018_TTTEEE_tSZ_nu0

  bandpowers_data = nm.loadtxt(add_data_path(data_path,dataset.SPT3G_2018_TTTEEE_bandpower_file))
  #bandpowers_data.shape=(N_b_0_EE,-1)
  lkl_grp["SPT3G_2018_TTTEEE_bandpower"] = bandpowers_data.flat[:]

  cov = nm.loadtxt(add_data_path(data_path,dataset.SPT3G_2018_TTTEEE_covariance_matrix))
  lkl_grp["SPT3G_2018_TTTEEE_covariance_matrix"] = cov.flat[:]

  beam_cov = nm.loadtxt(add_data_path(data_path,dataset.SPT3G_2018_TTTEEE_beam_covariance_matrix))
  lkl_grp["SPT3G_2018_TTTEEE_beam_covariance_matrix"] = beam_cov.flat[:]

  full_windows = nm.zeros((N_s_0,1+lmax-lmin,N_b_0_EE))
  #print(full_windows.shape)
  for i in range(N_b_0_EE):
    full_windows[:,:,i] = (nm.loadtxt(add_data_path(data_path+"../../",dataset.SPT3G_2018_TTTEEE_window_folder)+"window_%d.txt"%(i+1))[:,1:]).T
  lkl_grp["SPT3G_2018_TTTEEE_windows"] = full_windows.flat[:]

  central_freq = nm.loadtxt(add_data_path(data_path,dataset.SPT3G_2018_TTTEEE_central_frequency_file)).T
  lkl_grp["SPT3G_2018_TTTEEE_central_frequency"] = central_freq.flat[:]

  tSZ = nm.loadtxt(add_data_path(data_path,dataset.SPT3G_2018_TTTEEE_tSZ_template_file))[:,1:]
  kSZ = nm.loadtxt(add_data_path(data_path,dataset.SPT3G_2018_TTTEEE_kSZ_template_file))[:,1:]

  lkl_grp["SPT3G_2018_TTTEEE_tSZ_template"] = tSZ.flat[:]
  lkl_grp["SPT3G_2018_TTTEEE_kSZ_template"] = kSZ.flat[:]

  lkl_grp["SPT3G_2018_TTTEEE_tSZ_cosmology_scaling"] = dataset.str(default="F").SPT3G_2018_TTTEEE_tSZ_cosmology_scaling.lower() in ("1","t","true",".true.")
  lkl_grp["SPT3G_2018_TTTEEE_kSZ_cosmology_scaling"] = dataset.str(default="F").SPT3G_2018_TTTEEE_kSZ_cosmology_scaling.lower() in ("1","t","true",".true.")

  cal_cov = nm.loadtxt(add_data_path(data_path,dataset.SPT3G_2018_TTTEEE_cal_covariance_matrix))
  lkl_grp["SPT3G_2018_TTTEEE_cal_covariance_matrix"] = cal_cov.flat[:]

  
  root_grp.close()

import sys
if __name__=="__main__":
  main(sys.argv)