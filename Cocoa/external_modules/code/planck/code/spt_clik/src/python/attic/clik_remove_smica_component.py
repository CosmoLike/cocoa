#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import clik.parobject as php
import clik
import clik.hpy as hpy
import clik.smicahlp as slp

def main(argv):
  pars = clik.miniparse(argv[1])
  
  inhf = hpy.File(pars.input_object)
  assert inhf["clik/lkl_0"].attrs["lkl_type"].lower() == "smica"
  inhf.close()
  
  hpy.copyfile(pars.input_object,pars.res_object)
  outhf = hpy.File(pars.res_object,"r+")
  lkl_grp = outhf["clik/lkl_0"]
  del lkl_grp["component_%d"%component_index]
  
  php.remove_selfcheck(root_grp = outhf["clik"])
  outhf.close()

import sys
if __name__=="__main__":
  main(sys.argv)