#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import clik
import clik.hpy as hpy
import clik.parobject as php
import numpy as nm

def main(argv):
  #pars = clik.miniparse(argv[1])

  hpy.copyfile(argv[1],argv[2])
  outhf = hpy.File(argv[2],"r+")
  lkl_grp = outhf["clik/lkl_0"]
  php.add_free_calib(lkl_grp,"A_planck")

  if "clik/check_param" in outhf:
    outhf["clik/check_param"] = nm.concatenate((outhf["clik/check_param"],[1]))
  #php.remove_selfcheck(root_grp=outhf["clik"])

    
  outhf.close()

    
import sys
if __name__=="__main__":
  main(sys.argv)