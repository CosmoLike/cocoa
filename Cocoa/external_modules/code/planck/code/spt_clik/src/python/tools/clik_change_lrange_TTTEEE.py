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

# clik_change_lrange_TTTEEE.py 35 650 

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

#GO: start
def cutarr2(inhf, outhf, olmin, olmax, lmin, lmax, nn, hascl):
	ns = hascl.sum()
	arr = inhf[nn]
	if (nn == "clik/lkl_0/Rq_hat" or nn == "clik/lkl_0/criterion_gauss_mask" or nn == "clik/lkl_0/wq"):
		arr.shape = (int(olmax+1-olmin), -1)
		narr = arr[lmin-olmin:lmax+1-olmin]
		outhf[nn] = narr.flat[:]
	else:
		arr.shape = (ns*(olmax+1-olmin),-1)
		narr = []
		for i in range(ns):
			narr += [arr[(i*(olmax+1-olmin)) + lmin-olmin:(i*(olmax+1-olmin)) + lmax+1-olmin]]
		outhf[nn] = nm.concatenate(narr).flat[:]
#GO: end

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
  #GO assert hascl[1:].sum()==0,"do not work yet on polar data"
  if (hascl[1:].sum() == 0):  #GO this is the TT case that was already programmed in 
	  mT = inhf["clik/lkl_0/m_channel_T"]
	    
	  ord = inhf["clik/lkl_0/criterion_gauss_ordering"]
	  ord.shape=(-1,2)
	  nmsk = inhf["clik/lkl_0/criterion_gauss_mask"]
	  nmsk.shape=(nb,mT,mT)

	  kp = []
	  mx = 0
	  for i,j in ord:
	    cur = nm.arange(nb)[nmsk[:,i,j]==1]
	    kp += [nm.where((cur<bmax+1) * (cur>bmin-1))[0]+mx]
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
	  
	  outhf["clik/lmax"] = [lmax,-1,-1,-1,-1,-1]
	  
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
  #GO: start
  elif (hascl[1:].sum() != 0):
    print("Starting new code #GO")
    mT = inhf["clik/lkl_0/m_channel_T"] #number of T channels : 100, 143, 217 GHz
    mP = inhf["clik/lkl_0/m_channel_P"] #number of P channels 
    ord = inhf["clik/lkl_0/criterion_gauss_ordering"] #what order to read the mask matrix below
    ord.shape = (-1,2) 
    nmsk = inhf["clik/lkl_0/criterion_gauss_mask"] #mask to show which channels are correlated at each bin
    #print(type(nb/(hascl.sum())))
    #print(type(mT+mP))
    #print(nb/(hascl.sum()))
    nmsk.shape = (int(nb/(hascl.sum())), mT+mP, mT+mP) #rearrange into nb/hascl.sum() matrices each of size mT+mP. Each matrix has a 0 or a 1 to indicated which channels are correlated at that bin

    kp = []
    mx = 0
    for i,j in ord:
        cur = nm.arange(nb/hascl.sum())[nmsk[:,i,j]==1] #indices of bins (out of the first nb/hascl.sum() only) that have the i,j correlation
        kp += [nm.where((cur<bmax+1) * (cur>bmin-1))[0] + mx] #keep only the ones between bmin and bmax; this limits the ell range
        mx += len(cur) #keep count of how many you kept to ensure monotonicity 
    kp = nm.concatenate(kp)

    siginv = inhf["clik/lkl_0/criterion_gauss_mat"]
    siginv.shape = (int(siginv.shape[0]**.5),-1)
    sig = nm.linalg.inv(siginv)
    del(siginv)
    nsig = ((sig[kp])[:,kp])*1.
    del(sig)
    nsiginv = nm.linalg.inv(nsig)
    del(nsig)

	#WRITE THE NEW LIKELIHOOD
    hpy.copyfile(lklfile,outfile)
    outhf = hpy.File(outfile,"r+")

    outhf["clik/lmax"] = [lmax,lmax,-1,lmax,-1,-1]

    outhf["clik/lkl_0/lmin"] = lmin
    outhf["clik/lkl_0/lmax"] = lmax

    outhf["clik/lkl_0/nbins"] = hascl.sum()*(bmax + 1 - bmin)

    #We need to reset the l value to ensure that blmin[i+1]-blmax[i]=1
    ns = hascl.sum()
    blmin.shape = (nb, -1)
    blmax.shape = (nb, -1)
    nblmin = []
    nblmax = []
    for i in range(ns):
        if len(nblmin) == 0:
            offset = lmin-olmin
        else:
            offset = blmin[i*int(nb/ns) + bmin] - 1 - nblmax[-1][-1]
        #print(type(nb/ns))
        #print(nb/ns)
        nblmin += [blmin[i*int(nb/ns) + bmin:i*int(nb/ns) + bmax+1] - offset]
        nblmax += [blmax[i*int(nb/ns) + bmin:i*int(nb/ns) + bmax+1] - offset]

    outhf["clik/lkl_0/bin_lmin"] = nm.concatenate(nblmin).flat[:]
    outhf["clik/lkl_0/bin_lmax"] = nm.concatenate(nblmax).flat[:]

    cutarr2(inhf,outhf,olmin,olmax,lmin,lmax,"clik/lkl_0/bin_ws", hascl.sum())

    cutarr2(inhf,outhf,0,nb/hascl.sum()-1,bmin,bmax,"clik/lkl_0/Rq_hat", hascl.sum())

    cutarr2(inhf,outhf,0,nb/hascl.sum()-1,bmin,bmax,"clik/lkl_0/wq", hascl.sum())
    cutarr2(inhf,outhf,0,nb/hascl.sum()-1,bmin,bmax,"clik/lkl_0/criterion_gauss_mask", hascl.sum())
    outhf["clik/lkl_0/criterion_gauss_mat"] = nsiginv.flat[:]

    for i in range(1,inhf["clik/lkl_0/n_component"]):
        if "lmin" not in inhf["clik/lkl_0/component_%d"%i]:
            continue
        outhf["clik/lkl_0/component_%d/lmin"%i]=lmin
        outhf["clik/lkl_0/component_%d/lmax"%i]=lmax

    php.remove_selfcheck(root_grp=outhf["clik"])
    outhf.close()


  #GO: end

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
  #if ty not in ("smica","gibbs_gauss"):
  #  print("can only change lmin and lmax for plik and commander TT likelihoods")
  #  sys.exit(-1)
  #assert ty in ["smica","gibbs_gauss"],"Cannot change lrange for likelihood type %s"%ty
  fnc = globals()["change_%s"%ty]
  fnc(inhf,lklfile,outfile,lmin,lmax)
  
    
import sys
if __name__=="__main__":
  main(sys.argv)
