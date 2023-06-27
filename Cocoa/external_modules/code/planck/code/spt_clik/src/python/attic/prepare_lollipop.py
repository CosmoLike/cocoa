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

  hascl = nm.array(pars.bool_array.hascl,dtype=nm.int)
  lkl_grp = php.add_lkl_generic(root_grp,"lollipop",1,hascl,pars.int.lmax,pars.int.lmin)
  
  import tempfile
  dr = tempfile.mkdtemp()
  import os
  print(dr)
  import shutil
  shutil.copy(pars.data.strip(),dr+"/lollipop_data")    
  shutil.copy(pars.fiducial.strip(),dr+"/lollipop_fid")    
  shutil.copy(pars.invcov.strip(),dr+"/lollipop_invcov")    

  php.add_external_data(dr,lkl_grp,tar=True)

  hf.close()

  shutil.rmtree(dr)


import sys
if __name__=="__main__":
  main(sys.argv)