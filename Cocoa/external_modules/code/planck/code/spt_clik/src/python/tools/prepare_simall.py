#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import numpy as nm
import clik.parobject as php
import clik

def main(argv):
  pars = clik.miniparse(argv[1])
  
  root_grp,hf = php.baseCreateParobject(pars.res_object)
  hascl = nm.array(pars.bool_array.has_cl,dtype=nm.int)
  hascl = nm.array(hascl,dtype=nm.int)
  lmax,lmin = pars.int(default=30).lmax,pars.int(default=2).lmin
  print lmin,lmax
  lkl_grp = php.add_lkl_generic(root_grp,"simall",1,hascl,pars.int(default=30).lmax,pars.int(default=2).lmin)
  
  php.add_pid(lkl_grp,"simall_EE_BB_TE")

  lkl_grp["lmin"] = 2
  lkl_grp["nell"] = 29
  lkl_grp["lmax"] = 30
  if (hascl[1]==1): lkl_grp["nstepsEE"] = pars.int(default=3000).nstepsEE
  if (hascl[2]==1): lkl_grp["nstepsBB"] = pars.int(default=1000).nstepsBB
  if (hascl[3]==1): lkl_grp["nstepsTE"] = pars.int(default=10000).nstepsTE
  if (hascl[1]==1): lkl_grp["stepEE"] = pars.float(default=0.0001).stepEE
  if (hascl[2]==1): lkl_grp["stepBB"] = pars.float(default=0.0001).stepBB
  if (hascl[3]==1): lkl_grp["stepTE"] = pars.float(default=0.001).stepTE
  if (hascl[1]==1): lkl_grp["probEE"] = nm.fromfile(pars.fileEE)
  if (hascl[2]==1): lkl_grp["probBB"] = nm.fromfile(pars.fileBB)
  if (hascl[3]==1): lkl_grp["probTE"] = nm.fromfile(pars.fileTE)

  if "input_cl" in pars:
    cls = nm.loadtxt(pars.input_cl)
    a=nm.zeros(31*hascl.sum())
    ind=0
    for isp in nm.arange(4):
	if (hascl[isp]==1):
#		print isp
		ind+=2
		for il in nm.arange(29):
			l=il+2.0
			lnorm=l*(l+1.0)/2./nm.pi
			a[ind]=cls[il,isp+1]/lnorm
#                        print ind,a[ind],cls[il,isp+1]
			ind+=1

    res = php.add_selfcheck(pars.res_object,a)
    print res
  
import sys
if __name__=="__main__":
  main(sys.argv)
