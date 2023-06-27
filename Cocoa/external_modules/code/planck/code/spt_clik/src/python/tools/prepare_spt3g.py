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

def main(argv):
  pars = clik.miniparse(argv[1])
  
  root_grp,hf = php.baseCreateParobject(pars.res_object)
  hascl = [0]*6
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

  
  lmin = dataset.int.SPT3G_Y1_EETE_window_l_min
  lmax = dataset.int.SPT3G_Y1_EETE_window_l_max

  lkl_grp = php.add_lkl_generic(root_grp,"spt3g",1,hascl,lmax+1,0)
  lkl_grp.attrs["version"]=1

  lkl_grp["SPT3G_Y1_EETE_window_l_min"] = dataset.int.SPT3G_Y1_EETE_window_l_min
  lkl_grp["SPT3G_Y1_EETE_window_l_max"] = dataset.int.SPT3G_Y1_EETE_window_l_max

  with open(pars.ini_path) as ini :
    for l in ini:
      if "SPT3G_Y1_EETE_aberration_coefficient" in l:
        lkl_grp["SPT3G_Y1_EETE_aberration_coefficient"] = float(l.split("=")[-1].strip())
        break
  lkl_grp["SPT3G_Y1_EETE_super_sample_lensing"] = int(dataset.bool(default=True).SPT3G_Y1_EETE_super_sample_lensing)
  lkl_grp["SPT3G_Y1_EETE_poisson_sources"] = dataset.bool(default=True).SPT3G_Y1_EETE_poisson_sources
  lkl_grp["SPT3G_Y1_EETE_dust"] = int(dataset.bool(default=True).SPT3G_Y1_EETE_dust)
  lkl_grp["SPT3G_Y1_EETE_radio_galaxies_nu0"] = dataset.float(default=150).SPT3G_Y1_EETE_radio_galaxies_nu0
  lkl_grp["SPT3G_Y1_EETE_dsfg_nu0"] = dataset.float(default=150).SPT3G_Y1_EETE_dsfg_nu0
  lkl_grp["SPT3G_Y1_EETE_dust_nu0"] = dataset.float(default=150).SPT3G_Y1_EETE_dust_nu0
  lkl_grp["SPT3G_Y1_EETE_bin_min"] = dataset.int(default=0).SPT3G_Y1_EETE_bin_min
  lkl_grp["SPT3G_Y1_EETE_bin_max"] = dataset.int(default=44).SPT3G_Y1_EETE_bin_max
  lkl_grp["SPT3G_Y1_EETE_cal_prior"] = int(dataset.bool(default=False).SPT3G_Y1_EETE_cal_prior)
  lkl_grp["SPT3G_Y1_EETE_Kappa_mean"] = dataset.float(default=0).SPT3G_Y1_EETE_Kappa_mean
  lkl_grp["SPT3G_Y1_EETE_Kappa_sigma"] = dataset.float(default=0.00045).SPT3G_Y1_EETE_Kappa_sigma
  lkl_grp["SPT3G_Y1_EETE_Alpha_Dust_EE_prior"] = int(dataset.bool(default=False).SPT3G_Y1_EETE_Alpha_Dust_EE_prior)
  lkl_grp["SPT3G_Y1_EETE_Alpha_Dust_EE_prior_mean"] = dataset.float(default=-2.42).SPT3G_Y1_EETE_Alpha_Dust_EE_prior_mean
  lkl_grp["SPT3G_Y1_EETE_Alpha_Dust_EE_prior_sigma"] = dataset.float(default=0.02).SPT3G_Y1_EETE_Alpha_Dust_EE_prior_sigma
  lkl_grp["SPT3G_Y1_EETE_Alpha_Dust_TE_prior"] = int(dataset.bool(default=False).SPT3G_Y1_EETE_Alpha_Dust_TE_prior)
  lkl_grp["SPT3G_Y1_EETE_Alpha_Dust_TE_prior_mean"] = dataset.float(default=-2.42).SPT3G_Y1_EETE_Alpha_Dust_TE_prior_mean
  lkl_grp["SPT3G_Y1_EETE_Alpha_Dust_TE_prior_sigma"] = dataset.float(default=0.02).SPT3G_Y1_EETE_Alpha_Dust_TE_prior_sigma
  lkl_grp["SPT3G_Y1_EETE_Beam_Covariance_Scaling"] = dataset.float(default=1.0).SPT3G_Y1_EETE_Beam_Covariance_Scaling

  bandpowers = nm.loadtxt(add_data_path(data_path,dataset.SPT3G_Y1_EETE_bandpower_file))[:,1:]
  lkl_grp["SPT3G_Y1_EETE_bandpower"] = bandpowers.flat[:]

  bdp_order = dict([(b,a) for a,b in enumerate(dataset.str_array.SPT3G_Y1_EETE_bandpower_file_order)])
  bdp_select = nm.array([bdp_order[v] for v in dataset.str_array.SPT3G_Y1_EETE_spectra_to_fit])
  lkl_grp["SPT3G_Y1_EETE_spectra_to_fit"] = bdp_select+1
  lkl_grp["SPT3G_Y1_EETE_spectra_to_fit_kind"] = [3 if "T" in v else 1 for v in dataset.str_array.SPT3G_Y1_EETE_spectra_to_fit]

  N_b_0 = 44
  N_s_0 = 12
  bmin = dataset.int(default=0).SPT3G_Y1_EETE_bin_min-1
  bmax = dataset.int(default=44).SPT3G_Y1_EETE_bin_max
  N_b = bmax-bmin
  N_s = len(bdp_select)
  
  cov = nm.loadtxt(add_data_path(data_path,dataset.SPT3G_Y1_EETE_covariance_matrix))
  cov.shape = (N_b_0*N_s_0,N_b_0*N_s_0)
  cov_order = dict([(b,a) for a,b in enumerate(dataset.str_array.SPT3G_Y1_EETE_covariance_matrix_order)])
  cov_select = [cov_order[v] for v in dataset.str_array.SPT3G_Y1_EETE_spectra_to_fit]
  ncov = nm.zeros((N_b*N_s,N_b*N_s))
  for i in range(N_s):
    for j in range(i,N_s):
      # select block
      fromblock = cov[cov_select[i]*N_b_0:(cov_select[i]+1)*N_b_0,cov_select[j]*N_b_0:(cov_select[j]+1)*N_b_0]
      ncov[i*N_b:(i+1)*N_b,j*N_b:(j+1)*N_b] = fromblock[bmin:bmax,:][:,bmin:bmax]
      ncov[j*N_b:(j+1)*N_b,i*N_b:(i+1)*N_b] = ncov[i*N_b:(i+1)*N_b,j*N_b:(j+1)*N_b].T
  lkl_grp["SPT3G_Y1_EETE_covariance_matrix"] = ncov.flat[:]

  beamcov = nm.loadtxt(add_data_path(data_path,dataset.SPT3G_Y1_EETE_beam_covariance_matrix))
  beamcov.shape = (N_b_0*N_s_0,N_b_0*N_s_0)
  beamcov_order = dict([(b,a) for a,b in enumerate(dataset.str_array.SPT3G_Y1_EETE_beam_covariance_matrix_order)])
  beamcov_select = [beamcov_order[v] for v in dataset.str_array.SPT3G_Y1_EETE_spectra_to_fit]
  nbeamcov = nm.zeros((N_b*N_s,N_b*N_s))
  for i in range(N_s):
    for j in range(i,N_s):
      # select block
      fromblock = beamcov[beamcov_select[i]*N_b_0:(beamcov_select[i]+1)*N_b_0,beamcov_select[j]*N_b_0:(beamcov_select[j]+1)*N_b_0]
      nbeamcov[i*N_b:(i+1)*N_b,j*N_b:(j+1)*N_b] = fromblock[bmin:bmax,:][:,bmin:bmax]
      nbeamcov[j*N_b:(j+1)*N_b,i*N_b:(i+1)*N_b] = nbeamcov[i*N_b:(i+1)*N_b,j*N_b:(j+1)*N_b].T
  lkl_grp["SPT3G_Y1_EETE_beam_covariance_matrix"] = nbeamcov.flat[:]

  mapset = set()
  for sp  in dataset.str_array.SPT3G_Y1_EETE_spectra_to_fit:
    mapset.update(sp.split("x"))

  calrow = nm.array([int(v in mapset) for v in ["90_T","150_T","220_T","90_E","150_E","220_E"]])

  lkl_grp["SPT3G_Y1_EETE_cal_row_use"] = calrow

  idscalT = {"90":15,"150":16,"220":17}
  idscalE = {"90":18,"150":19,"220":20}
  
  spectra_to_fit_tcal_indices = [[idscalT[v.split("_")[0]] for v in sp.split("x")] for sp in dataset.str_array.SPT3G_Y1_EETE_spectra_to_fit]
  spectra_to_fit_ecal_indices = [[idscalE[v.split("_")[0]] for v in sp.split("x")] for sp in dataset.str_array.SPT3G_Y1_EETE_spectra_to_fit]

  lkl_grp["SPT3G_Y1_EETE_spectra_to_fit_tcal"] = nm.array(spectra_to_fit_tcal_indices).flat[:]
  lkl_grp["SPT3G_Y1_EETE_spectra_to_fit_ecal"] = nm.array(spectra_to_fit_ecal_indices).flat[:]

  calcov = nm.loadtxt(add_data_path(data_path,dataset.SPT3G_Y1_EETE_cal_covariance_matrix))
  ncalcov = calcov[calrow==1][:,calrow==1]
  icalcov = nm.linalg.inv(ncalcov)

  lkl_grp["SPT3G_Y1_EETE_cal_inv_covariance_matrix"] = icalcov.flat[:]

  central_freq = dict([[k,v] for k,v in zip(dataset.str_array.SPT3G_Y1_EETE_central_frequency_file_order, nm.loadtxt(add_data_path(data_path,dataset.SPT3G_Y1_EETE_central_frequency_file)))])
  
  spectra_to_fit_nu_eff = [[central_freq[v.split("_")[0]] for v in sp.split("x")] for sp in dataset.str_array.SPT3G_Y1_EETE_spectra_to_fit]
  #spectra_to_fit_nu_eff = nm.arange(24)*1.
  lkl_grp["SPT3G_Y1_EETE_central_frequency"] = nm.array(spectra_to_fit_nu_eff).flat[:]
  
  import re
  idspoi = {"90x90":2, "90x150":3, "90x220":4, "150x150":5, "150x220":6, "220x220":7}
  spectra_to_fit_poisson_indices = [ idspoi[re.compile("_[T,E]").sub("",sp)] for sp in dataset.str_array.SPT3G_Y1_EETE_spectra_to_fit]
  lkl_grp["SPT3G_Y1_EETE_spectra_to_fit_poisson_indices"] = nm.array(spectra_to_fit_poisson_indices).flat[:]

  full_windows = nm.zeros((N_b_0,1+lmax-lmin,N_s_0))
  win_order = dict([(b,a) for a,b in enumerate(dataset.str_array.SPT3G_Y1_EETE_window_folder_order)])
  win_select = [win_order[v] for v in dataset.str_array.SPT3G_Y1_EETE_spectra_to_fit]
  for i in range(N_b_0):
    sel_windows = nm.loadtxt(add_data_path(data_path+"../../",dataset.SPT3G_Y1_EETE_window_folder)+"window_%d.txt"%(i+1))[:,1:]
    for j in range(N_s_0):
      full_windows[i,:,j] = sel_windows[:,win_select[j]]
  lkl_grp["SPT3G_Y1_EETE_windows"] = nm.array(full_windows).flat[:]


  root_grp.close()

import sys
if __name__=="__main__":
  main(sys.argv)