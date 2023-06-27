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
  
  mspec_folder = pars.mspec_folder
  import cosmoslik as K
  import mspec
  print(mspec.__file__)
  script = K.load_script(osp.join(mspec_folder,'script.py'),cls_set_externally=True)
  lmax = [max([-1]+[lmax for ((x1,_),(x2,_)),(_,lmax) in list(script.params.mspec.use.items()) if x1+x2==x]) for x in ['TT','EE','BB','TE','TB','EB']]
  nuis = list(script.get_sampled().keys())

  hascl = [int(lm!=-1) for lm in lmax]
  hascl = nm.array(hascl,dtype=nm.int)
  lm = max(lmax)
  
  lkl_grp = php.add_lkl_generic(root_grp,"mspec",1,hascl,lm,0)
  
  lkl_grp.attrs["n_nuisance"] = len(nuis)
  lkl_grp.attrs["nuisance"] = php.pack256(*nuis)

  php.add_pid(lkl_grp,"")

  import tempfile
  dr = tempfile.mkdtemp()
  import os
  print(dr)
  import shutil
  shutil.copytree(mspec_folder.strip(),dr+"/mspec_data")    
  
  php.add_external_data(dr,lkl_grp,tar=True)

  hf.close()

  shutil.rmtree(dr)


import sys
if __name__=="__main__":
  main(sys.argv)