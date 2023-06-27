#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import numpy as nm
import clik.parobject as php
import clik
import os.path as osp

def main(argv):
  pars = clik.miniparse(argv[1])
  test_cl = nm.loadtxt(osp.join(pars.lowlike_data,"data/test_cls_v5.dat"))
  
  mcl = nm.zeros((4,1201),dtype=nm.double)
  llp1s2pi = nm.arange(1201)*nm.arange(1,1202)/2./nm.pi
  mcl[:,2:] = (test_cl[:1201-2,1:5].T)/llp1s2pi[2:]
  
  #print has_tt,has_te,ttmin,ttmax,temin,temax
  
  root_grp,hf = php.baseCreateParobject(pars.res_object)
  hascl = [0]*6
  hascl[0] = 1
  hascl[1:4] = [1]*3
  hascl = nm.array(hascl,dtype=nm.int)
  #print hascl
  
  lmin = 0
  lmax = 32
  
  mcl = (nm.compress(hascl[:4],mcl[:,:lmax+1],0)).flat[:]
  
  lkl_grp = php.add_lkl_generic(root_grp,"lowlike",1,hascl,lmax,lmin)
  
  lkl_grp.attrs["ttmin"] = 2
  lkl_grp.attrs["temin"] = 2
  lkl_grp.attrs["ttmax"] = 32
  lkl_grp.attrs["temax"] = 32
  
  lkl_grp.attrs["use_gibbs"] = 1
  lkl_grp.attrs["use_lowl_pol"] = 1
  
  lkl_grp.attrs["use_wmap_pol"] = pars.int(default=0).use_wmap_pol
  
  #lkl_grp.attrs["external_dir"] = osp.realpath(pars.lowlike_data)
  php.add_external_data(osp.realpath(pars.lowlike_data),lkl_grp,tar=bool(pars.int(default=1).include))

  php.add_pid(lkl_grp,pars.str(default="").pid)

  hf.close()
  
  if hasattr(clik,"clik"):
    res = php.add_selfcheck(pars.res_object,mcl)
    print("lkl for init cl %g"%res)
  
  if "cl_save" in pars:
    f=open(pars.cl_save,"w")
    for ci in mcl:
      print(ci, file=f)
    f.close()

import sys
if __name__=="__main__":
  main(sys.argv)