#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import numpy as nm
import numpy.random as ra
import numpy.linalg as la
import clik.parobject as php
import clik
import re
import os.path as osp
import os
import shutil
import clik.hpy as hpy

def change_simlow(inhf,lklfile,outfile,lmin,lmax):
  hpy.copyfile(lklfile,outfile)
  outhf = hpy.File(outfile,"r+")
  outhf["clik/lmax"] = [-1,lmax,-1,-1,-1,-1]
  outhf["clik/lkl_0/lmin"] = lmin
  outhf["clik/lkl_0/lmax"] = lmax
  php.remove_selfcheck(root_grp=outhf["clik"])
  outhf.close()


def change_gibbs_gauss(inhf,lklfile,outfile,lmin,lmax):
  olmin = inhf["clik/lkl_0/lmin"]
  if lmin == -1 or lmin>lmax:
    lmin = inhf["clik/lkl_0/lmin"]
  if lmax == -1:
    lmax = inhf["clik/lkl_0/lmax"]
  if lmax>249:
    print("not possible")
    sys.exit(-1)
  hpy.copyfile(lklfile,outfile)
  outhf = hpy.File(outfile,"r+")
  outhf["clik/lmax"] = [lmax,-1,-1,-1,-1,-1]
  outhf["clik/lkl_0/lmin"] = lmin
  outhf["clik/lkl_0/lmax"] = lmax
  outhf["clik/lkl_0/delta_l"] = lmax
  php.remove_selfcheck(root_grp=outhf["clik"])
  outhf.close()

def cutarr(inhf,outhf,olmin,olmax,lmin,lmax,nn):
  arr = inhf[nn]
  arr.shape=(olmax+1-olmin,-1)
  narr = arr[lmin-olmin:lmax+1-olmin]
  outhf[nn] = narr.flat[:]

def change_plik_cmbonly(inhf,lklfile,outfile,lmin,lmax):
  olmin = 30
  olmax = 2508
  if lmin == -1 or lmin>lmax or lmin<olmin:
    lmin = olmin
  if lmax == -1 or lmax>olmax:
    lmax = olmax
  nb = 215
  blmin = nm.loadtxt(inhf._name+"/clik/lkl_0/_external/blmin.dat")
  blmax = nm.loadtxt(inhf._name+"/clik/lkl_0/_external/blmax.dat")
  bmin = nm.argmin((blmin+olmin-lmin)**2)
  bmax = nm.argmin((blmax+olmin-lmax)**2)

  lmin = blmin[bmin]+olmin
  lmax = blmax[bmax]+olmin

  print("restrict to %d %d [%d %d]"%(lmin,lmax,bmin,bmax))

  hpy.copyfile(lklfile,outfile)
  outhf = hpy.File(outfile,"r+")
  outhf["clik/lkl_0/bin_min_tt"] = bmin+1
  outhf["clik/lkl_0/bin_max_tt"] = bmax+1
  outhf["clik/lkl_0/bin_min_te"] = 1
  outhf["clik/lkl_0/bin_max_te"] = 199
  outhf["clik/lkl_0/bin_min_ee"] = 1
  outhf["clik/lkl_0/bin_max_ee"] = 199

  php.remove_selfcheck(root_grp=outhf["clik"])
  outhf.close()

def change_smica(inhf,lklfile,outfile,lmin,lmax):
  olmin = inhf["clik/lkl_0/lmin"]
  olmax = inhf["clik/lkl_0/lmax"]
  if lmin == -1 or lmin>lmax or lmin<olmin:
    lmin = olmin
  if lmax == -1 or lmax>olmax:
    lmax = olmax
  
  nb = inhf["clik/lkl_0/nbins"]

  blmin = inhf["clik/lkl_0/bin_lmin"]
  blmax = inhf["clik/lkl_0/bin_lmax"]

  bmin = nm.argmin((blmin+olmin-lmin)**2)
  bmax = nm.argmin((blmax+olmin-lmax)**2)

  lmin = blmin[bmin]+olmin
  lmax = blmax[bmax]+olmin

  print("restrict to %d %d"%(lmin,lmax))

  hascl = inhf["clik/lkl_0/has_cl"]
  assert hascl.sum()==1, "do not work for joint T+P (incuding TE) yet"

  mT = inhf["clik/lkl_0/m_channel_T"]
  mP = inhf["clik/lkl_0/m_channel_P"]
    
  ord = inhf["clik/lkl_0/criterion_gauss_ordering"]
  ord.shape=(-1,2)
  nmsk = inhf["clik/lkl_0/criterion_gauss_mask"]
  nmsk.shape=(nb,mT+mP,mT+mP)

  kp = []
  mx = 0
  for i,j in ord:
    cur = nm.arange(nb)[nmsk[:,i,j]==1]
    kp += [nm.where((cur<bmax+1) * (cur>bmin-1))[0]+mx]
    mx += len(cur)
  kp = nm.concatenate(kp)
  
  siginv = inhf["clik/lkl_0/criterion_gauss_mat"]
  siginv.shape=(siginv.shape[0]**.5,-1)
  sig = nm.linalg.inv(siginv)
  del(siginv)
  nsig = ((sig[kp])[:,kp])*1.
  del(sig)
  nsiginv = nm.linalg.inv(nsig)
  del(nsig)

  hpy.copyfile(lklfile,outfile)
  outhf = hpy.File(outfile,"r+")
  
  lmaxo = hascl*(lmax+1)-1
  outhf["clik/lmax"] = lmaxo
  
  outhf["clik/lkl_0/lmin"] = lmin
  outhf["clik/lkl_0/lmax"] = lmax
  
  outhf["clik/lkl_0/nbins"] = bmax + 1 - bmin
  cutarr(inhf,outhf,0,nb-1,bmin,bmax,"clik/lkl_0/bin_lmin")
  outhf["clik/lkl_0/bin_lmin"] = outhf["clik/lkl_0/bin_lmin"]-(lmin-olmin)
  cutarr(inhf,outhf,0,nb-1,bmin,bmax,"clik/lkl_0/bin_lmax")
  outhf["clik/lkl_0/bin_lmax"] = outhf["clik/lkl_0/bin_lmax"]-(lmin-olmin)
  cutarr(inhf,outhf,olmin,olmax,lmin,lmax,"clik/lkl_0/bin_ws")

  cutarr(inhf,outhf,0,nb-1,bmin,bmax,"clik/lkl_0/Rq_hat")

  cutarr(inhf,outhf,0,nb-1,bmin,bmax,"clik/lkl_0/wq")
  cutarr(inhf,outhf,0,nb-1,bmin,bmax,"clik/lkl_0/criterion_gauss_mask")
  outhf["clik/lkl_0/criterion_gauss_mat"] = nsiginv.flat[:]

  for i in range(1,inhf["clik/lkl_0/n_component"]):
    if "lmin" not in inhf["clik/lkl_0/component_%d"%i]:
      continue
    outhf["clik/lkl_0/component_%d/lmin"%i]=lmin
    outhf["clik/lkl_0/component_%d/lmax"%i]=lmax

  php.remove_selfcheck(root_grp=outhf["clik"])
  outhf.close()


def main(argv):

  if len(sys.argv)!=5:
    print("usage : %s input_clik lmin lmax output_clik\n  copy input_clik to output_clik, change lmin and lmax of output_clik.\n  lmin and lmax can be set to -1 to keep the input_clik values.\n  input_clik, must be a plik or commander file.\n  if input_clik is a binned plik file, the effective lmin and lmax will be the set to the closest bins."%osp.basename(sys.argv[0]))
    sys.exit(0)
  lklfile = sys.argv[1]
  lmin = int(sys.argv[2])
  lmax = int(sys.argv[3])
  outfile = sys.argv[4]

  inhf = hpy.File(lklfile)
  ty = inhf["clik/lkl_0/lkl_type"]
  if ty not in ("smica","gibbs_gauss","plik_cmbonly","simlow"):
    print("can only change lmin and lmax for plik, plik_lite, and commander TT likelihoods")
    sys.exit(-1)
  assert ty in ["smica","gibbs_gauss","plik_cmbonly","simlow"],"Cannot change lrange for likelihood type %s"%ty
  fnc = globals()["change_%s"%ty]
  fnc(inhf,lklfile,outfile,lmin,lmax)
  
    
import sys
if __name__=="__main__":
  main(sys.argv)