#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import numpy as nm
import clik.parobject as php
import clik
import os.path as osp

def main(argv):
  pars = clik.miniparse(argv[1])
  
  root_grp,hf = php.baseCreateParobject(pars.res_object)
  hascl = [0]*6
  hascl[1] = 1
  #hascl[0] = 1
  hascl = nm.array(hascl,dtype=nm.int)
  lmax,lmin = pars.int(default=30).lmax,pars.int(default=2).lmin
  print lmin,lmax
  lkl_grp = php.add_lkl_generic(root_grp,"simlow",1,hascl,pars.int(default=30).lmax,pars.int(default=2).lmin)
  
  php.add_pid(lkl_grp,"simlow_ee_2_32_2016_03_31")

  lkl_grp["lmin"] = 2
  lkl_grp["nell"] = 30
  lkl_grp["nstep"] = 4000
  lkl_grp["step"] = 1e-4
  lkl_grp["tab"] = nm.fromfile(pars.tabname)

  hf.close()

  if "input_cl" in pars:
    cls = nm.loadtxt(pars.input_cl)
    res = php.add_selfcheck(pars.res_object,cls[:lmax+1])


import sys
if __name__=="__main__":
  main(sys.argv)