#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import numpy as nm
import clik.parobject as php
import clik
import os.path as osp
import clik.hpy as hpy

def read_binning(pars,name,fl):
  
  assert len(set(getattr(pars.str_array,"%s_out_order"%name)))==1,"bad out_order"

  lmax = pars.int.cl_lmax
  nb = pars.int.nbins
  bmin = pars.int.use_min -1
  bmax = pars.int.use_max 
  rnb = bmax - bmin
  n_in = len(getattr(pars.str_array,"%s_in_order"%name))

  rord = list(range(n_in))
  if n_in!=1:
    iord = getattr(pars.str_array,"%s_in_order"%name)
    rord = [iord.index(v) for v in ["PP","TT","EE","BB","TE","TB","EB"] if v in iord]

  mat = nm.zeros((rnb,(lmax+1)*n_in))

  for i in range(bmin,bmax):
    print(osp.join(osp.dirname(fl),getattr(pars.str,"%s_files"%name)%(i+1)))
    bi = nm.loadtxt(osp.join(osp.dirname(fl),getattr(pars.str,"%s_files"%name)%(i+1)))
    for ioor,oor in enumerate(rord):
      mat[i,(bi[:,0]+(lmax+1)*ioor).astype(nm.int)] = bi[:,oor+1]
  return mat

def main(argv):
  pars = clik.miniparse(argv[1])
  
  root = hpy.File(pars.res_object,"w")
  root.create_group("clik_lensing")
  root["clik_lensing/itype"]=4
  
  dars = clik.miniparse(pars.str.dataset)
  
  bins = read_binning(dars,"bin_window",pars.str.dataset)
  cors = read_binning(dars,"linear_correction_bin_window",pars.str.dataset)
  lmax = dars.int.cl_lmax
  nb = dars.int.nbins
  bmin = dars.int.use_min -1
  bmax = dars.int.use_max 
  hascl = nm.array([int(v in dars.str.linear_correction_bin_window_in_order.lower()) for v in ["tt","ee","bb","te","tb","eb"]])

  root["clik_lensing/nbins"] = bmax-bmin
  root["clik_lensing/lmax"] = lmax
  root["clik_lensing/hascl"] = hascl.astype(nm.int)

  root["clik_lensing/pp_hat"] = nm.loadtxt(osp.join(osp.dirname(pars.str.dataset),dars.cl_hat_file))[bmin:bmax,4]

  root["clik_lensing/bins"] = bins.flat[:]

  root["clik_lensing/siginv"] = nm.linalg.inv(nm.loadtxt(osp.join(osp.dirname(pars.str.dataset),dars.covmat_fiducial))[bmin:bmax,bmin:bmax]).flat[:]

  root["clik_lensing/cors"] = cors.flat[:]  

  if "linear_correction_fiducial_file" in dars:
    root["clik_lensing/cor0"] = nm.loadtxt(osp.join(osp.dirname(pars.str.dataset),dars.linear_correction_fiducial_file))[bmin:bmax,1]  
  else:
    root["clik_lensing/cor0"] = nm.loadtxt(osp.join(osp.dirname(pars.str.dataset),dars.linear_correction_fiducial))[bmin:bmax,1]  


  root["clik_lensing/renorm"] = pars.int(default=1).renorm
  root["clik_lensing/ren1"] = pars.int(default=1).ren1
  root["clik_lensing/has_calib"] = pars.int(default=1).has_calib

  clf = nm.loadtxt(pars.str.cl_fid)

  clrf = nm.zeros((1+hascl.sum(),lmax+1))
  ell = clf[:,0][:lmax+1-2]

  llmax = int(ell[-1]+1)

  clrf[0,2:llmax] = (clf[:,-1][:lmax+1-2] /  ((ell*(ell+1))*(ell*(ell+1)))) *2*nm.pi
  ii = 1
  print(hascl)
  if hascl[0]:
    clrf[ii,2:llmax] = clf[:,1][:lmax+1-2] /ell/(ell+1)*2*nm.pi
    ii+=1
    print(clrf[:,-3:])
  if hascl[1]:
    clrf[ii,2:llmax] = clf[:,3][:lmax+1-2] /ell/(ell+1)*2*nm.pi
    ii+=1
    print(clrf[:,-3:])
  if hascl[2]:
    clrf[ii,2:llmax] = clf[:,4][:lmax+1-2] /ell/(ell+1)*2*nm.pi
    ii+=1
    print(clrf[:,-3:])
  if hascl[3]:
    clrf[ii,2:llmax] = clf[:,2][:lmax+1-2] /ell/(ell+1)*2*nm.pi
    ii+=1
    print(clrf[:,-3:])

  print(clrf.shape)
  print(clrf[:,-3:])
  root["clik_lensing/cl_fid"] = clrf.flat[:]

  root.close()

import sys
if __name__=="__main__":
  main(sys.argv)