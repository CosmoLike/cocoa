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
  hascl[2] = 1
  hascl = nm.array(hascl,dtype=nm.int)
  
  lkl_grp = php.add_lkl_generic(root_grp,"bicep",1,hascl,600,0)
  
  php.add_pid(lkl_grp,"")

  import tempfile
  dr = tempfile.mkdtemp()
  import os
  print(dr)
  import shutil
  shutil.copytree(pars.cliksrc.strip()+"/src/bicep/windows",dr+"/windows")    
  shutil.copy(pars.cliksrc.strip()+"/src/bicep/B2_3yr_bpcm_no-sysuncer_20140314.txt",dr+"/")
  shutil.copy(pars.cliksrc.strip()+"/src/bicep/B2_3yr_likelihood_bandpowers_20140314.txt",dr+"/")
  
  php.add_external_data(dr,lkl_grp,tar=True)

  hf.close()

  shutil.rmtree(dr)


import sys
if __name__=="__main__":
  main(sys.argv)