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
  #hascl[0] = 1
  hascl = nm.array(hascl,dtype=nm.int)
  
  lkl_grp = php.add_lkl_generic(root_grp,"ftau",1,hascl,1,0)
  
  php.add_pid(lkl_grp,"")

  lkl_grp["tau_name"] = php.pack256(pars.str(default="tau").tau_name)
  tauP = nm.loadtxt(pars.str.tau_table,skiprows=pars.int(default=1).skip)
  d,f=tauP[:,0].searchsorted([pars.float(default=.015).tau_min,pars.float(default=.1).tau_max])
  tv = tauP[d:f+1,0]
  Pv = nm.log(tauP[d:f+1,1])
  lkl_grp["tau"] = tv
  lkl_grp["Ptau"] = Pv
  lkl_grp["n_tab"] = len(tv)

  
  hf.close()



import sys
if __name__=="__main__":
  main(sys.argv)