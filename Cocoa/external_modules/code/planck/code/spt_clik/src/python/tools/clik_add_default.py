#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import numpy as nm
import clik.parobject as php
import clik
import clik.hpy as hpy



def main(argv):
  pars = clik.miniparse(argv[1])
  
  clikl = clik.clik(pars.input_object)
  extn = clikl.extra_parameter_names
  
  name = pars.str_array.name
  loc = pars.float_array.loc
  assert len(name)==len(loc),"name and loc have different sizes"

  hpy.copyfile(pars.input_object,pars.res_object)
  outhf = hpy.File(pars.res_object,"r+")
    
  php.add_default(outhf["clik"],name,loc,extn)
  php.remove_selfcheck(root_grp=outhf["clik"])
  outhf.close()

    
import sys
if __name__=="__main__":
  main(sys.argv)