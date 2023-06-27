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


def cutarr(inhf,outhf,olmin,olmax,lmin,lmax,nn,inarr=None):
  if inarr is None:
    arr = inhf[nn]
  else:
    arr = inarr
  arr.shape=(olmax+1-olmin,-1)
  narr = arr[lmin-olmin:lmax+1-olmin]
  outhf[nn] = narr.flat[:]

def change_smica(inhf,lklfile,outfile,lmins,lmaxs,beg,nd):
  olmin = inhf["clik/lkl_0/lmin"]
  olmax = inhf["clik/lkl_0/lmax"]
  lmaxs = nm.array([olmax if (lm == -1 or lm>olmax) else lm for lm in lmaxs])
  lmins = nm.array([olmin if (lm == -1 or lm<olmin) else lm for lm in lmins])
  beg_notch = nm.array([lmx if (lm == -1 or lm>lmx) else lm for lm,lmx in zip(beg,lmaxs)])
  end_notch = nm.array([lmx if (lm == -1 or lm<lmx) else lm for lm,lmx in zip(nd,lmins)])
  
  kpp =  lmaxs>lmins
  
  nb = inhf["clik/lkl_0/nbins"]

  hascl = inhf["clik/lkl_0/has_cl"]
  mT = inhf["clik/lkl_0/m_channel_T"]
  mP = inhf["clik/lkl_0/m_channel_P"]  
  mE = mT*hascl[0]
  mB = mE + mP*hascl[1]
  m =  mB + mP*hascl[2]
  
  blmin = inhf["clik/lkl_0/bin_lmin"][:int(nb/nm.sum(hascl))]
  blmax = inhf["clik/lkl_0/bin_lmax"][:int(nb/nm.sum(hascl))]
  
  
  bmins = nm.array([nm.argmin((blmin+olmin-lm)**2) for lm in lmins])
  bmaxs = nm.array([nm.argmin((blmax+olmin-lm)**2) for lm in lmaxs])

  
  bbeg_notch = nm.array([nm.argmin((blmin+olmin-lm)**2) for lm in beg_notch])
  bend_notch = nm.array([nm.argmin((blmin+olmin-lm)**2) for lm in end_notch])

  lmins = nm.array([blmin[bm]+olmin for bm in bmins])
  lmaxs = nm.array([blmax[bm]+olmin for bm in bmaxs])

  
  beg_notch = nm.array([blmin[bm]+olmin for bm in bbeg_notch])
  end_notch = nm.array([blmax[bm]+olmin for bm in bend_notch])

  
  dnames = inhf["clik/lkl_0/dnames"]
  dnames = [dnames[i*256:(i+1)*256].strip("\0") for i in range(len(dnames)//256)]
  print("restrict to")
  cc = 0
  for a in range(3):
    aT = "TEB"[a]
    for ii in range(mT if a==0 else mP*hascl[a]):
      for b in range(3):
        bT = "TEB"[b]
        for jj in range(mT if b==0 else mP*hascl[b]):
          print("  %s%s %s%s lmin = %d, lmax = %d"%(aT,dnames[mT*(a!=0)+ii][:-1],bT,dnames[mT*(b!=0)+jj][:-1],lmins[cc]*kpp[cc],lmaxs[cc]*kpp[cc]),end="")
          if beg_notch[cc]*kpp[cc]<end_notch[cc]*kpp[cc]:
            print("[notch %d -> %d]"%(beg_notch[cc]*kpp[cc],end_notch[cc]*kpp[cc],),end="")
          print("")
          cc+=1


  lmin = min(lmins[kpp])
  lmax = max(lmaxs[kpp])
  bmin = min(bmins[kpp])
  bmax = max(bmaxs[kpp])
  
  #sys.exit(-1)

  ord = inhf["clik/lkl_0/criterion_gauss_ordering"]
  ord.shape=(-1,2)
  nmsk = inhf["clik/lkl_0/criterion_gauss_mask"]
  nmsk.shape=(int(nb/nm.sum(hascl)),m,m)

  kp = []
  mx = 0
  bmins.shape=(m,m)
  bmaxs.shape=(m,m)
  bbeg_notch.shape=(m,m)
  bend_notch.shape=(m,m)
  kpp.shape=(m,m)
  nnmsk = nmsk*1

  for i,j in ord:
    cur = nm.arange(int(nb/nm.sum(hascl)))[nmsk[:,i,j]==1]
    nnmsk[:bmins[i,j],i,j] = 0
    nnmsk[bmaxs[i,j]+1:,i,j] = 0
    nnmsk[bbeg_notch[i,j]:bend_notch[i,j]+1,i,j]=0
    if kpp[i,j]:
      kp += [nm.where((cur<bmaxs[i,j]+1) * (cur>bmins[i,j]-1) * (((cur>bbeg_notch[i,j]-1) * (cur<bend_notch[i,j]+1))==False))[0]+mx]
      #print "cur",cur+mx
      #print "kp",kp[-1]
      #print i,j,len(cur),len(kp[-1])
      #print (cur<bmaxs[i,j]+1)
      #print (cur>bmins[i,j]-1)
      #print (cur<bbeg_notch[i,j]+1),bbeg_notch[i,j]
      #print (cur>bend_notch[i,j]-1),bend_notch[i,j]
    else:
      nnmsk[:,i,j] = 0
    mx += len(cur)
    
  kp = nm.concatenate(kp)
  
  siginv = inhf["clik/lkl_0/criterion_gauss_mat"]
  siginv.shape=(int(siginv.shape[0]**.5),-1)
  sig = nm.linalg.inv(siginv)
  del(siginv)
  
  nsig = ((sig[kp])[:,kp])*1.
  del(sig)
  nsiginv = nm.linalg.inv(nsig)
  del(nsig)

  hpy.copyfile(lklfile,outfile)
  outhf = hpy.File(outfile,"r+")
  
  #outhf["clik/lmax"] = [lmax,-1,-1,-1,-1,-1]
  
  #outhf["clik/lkl_0/lmin"] = lmin
  #outhf["clik/lkl_0/lmax"] = lmax
  
  #outhf["clik/lkl_0/nbins"] = bmax + 1 - bmin
  #cutarr(inhf,outhf,0,nb-1,bmin,bmax,"clik/lkl_0/bin_lmin")
  #outhf["clik/lkl_0/bin_lmin"] = outhf["clik/lkl_0/bin_lmin"]-(lmin-olmin)
  #cutarr(inhf,outhf,0,nb-1,bmin,bmax,"clik/lkl_0/bin_lmax")
  #outhf["clik/lkl_0/bin_lmax"] = outhf["clik/lkl_0/bin_lmax"]-(lmin-olmin)
  #cutarr(inhf,outhf,olmin,olmax,lmin,lmax,"clik/lkl_0/bin_ws")

  #cutarr(inhf,outhf,0,nb-1,bmin,bmax,"clik/lkl_0/Rq_hat")

  #cutarr(inhf,outhf,0,nb-1,bmin,bmax,"clik/lkl_0/wq")
  #cutarr(inhf,outhf,0,nb-1,bmin,bmax,"clik/lkl_0/criterion_gauss_mask",nnmsk)
  outhf["clik/lkl_0/criterion_gauss_mask"] = nnmsk.flat[:]
  outhf["clik/lkl_0/criterion_gauss_mat"] = nsiginv.flat[:]

  #for i in range(1,inhf["clik/lkl_0/n_component"]):
  #  if "lmin" not in inhf["clik/lkl_0/component_%d"%i]:
  #    continue
  #  outhf["clik/lkl_0/component_%d/lmin"%i]=lmin
  #  outhf["clik/lkl_0/component_%d/lmax"%i]=lmax

  php.remove_selfcheck(root_grp=outhf["clik"])
  outhf.close()


def main(argv):

  if len(sys.argv)!=2:
    print("usage : %s parfile\n  copy input_clik to output_clik, change lmin and lmax of output_clik.\n  lmin and lmax can be set to -1 to keep the input_clik values.\n  input_clik, must be a plik or commander file.\n  if input_clik is a binned plik file, the effective lmin and lmax will be the set to the closest bins.\n  Here is an example parfile\n    #input and output clik files\n    input_clik = plik_dx11dr2_HM_v18_TT.clik\n\n    output_clik = plik_dx11dr2_HM_v18_TT_no143x143.clik\n\n    #set to -1 to use default value\n    #set to anything else to change lmin\n    lmin = -1    -1    -1    &\n           -1    -1    -1    &\n           -1    -1    -1   \n\n \n    #set to -1 to use default value\n    #set to anything else to change lmax\n    #set any non negative value inferior than lmin (0 is a good example) ro remove a cross spectra\n    lmax = -1    -1    -1    &\n           -1     0    -1    &\n           -1    -1    -1   \n\n "%osp.basename(sys.argv[0]))
    sys.exit(0)
  pars = clik.miniparse(argv[1])
  lklfile = pars.str.input_clik
  lmin = pars.int_array.lmin
  lmax = pars.int_array.lmax
  if "beg_notch" in pars:
    beg = pars.int_array.beg_notch
    nd = pars.int_array.end_notch
  else:
    beg = [-1]*len(lmin)
    nd = [-1]*len(lmin)

  outfile = pars.str.output_clik
  
  inhf = hpy.File(lklfile)
  ty = inhf["clik/lkl_0/lkl_type"]
  if ty not in ("smica"):
    print("can only change lmin and lmax for plik likelihoods")
    sys.exit(-1)
  assert ty in ["smica"],"Cannot change lrange for likelihood type %s"%ty
  fnc = globals()["change_%s"%ty]
  fnc(inhf,lklfile,outfile,lmin,lmax,beg,nd)
  
    
import sys
if __name__=="__main__":
  main(sys.argv)