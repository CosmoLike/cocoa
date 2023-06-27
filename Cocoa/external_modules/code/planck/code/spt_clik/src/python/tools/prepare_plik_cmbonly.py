#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import numpy as nm
import clik.parobject as php
import clik
import os.path as osp

def main(argv):
  pars = clik.miniparse(argv[1])
  #test_cl = nm.loadtxt(osp.join(pars.wmap_data,"data/v2a1s_best_lcdm_6000.txt"))
  
  #mcl = nm.zeros((4,1201),dtype=nm.double)
  #llp1s2pi = nm.arange(1201)*nm.arange(1,1202)/2./nm.pi
  #mcl[:,2:] = (test_cl[:1201-2,1:].T)/llp1s2pi[2:]
  
  root_grp,hf = php.baseCreateParobject(pars.res_object)
  
  hascl = nm.array(pars.bool_array.hascl,dtype=nm.int)
  
  lkl_grp = php.add_lkl_generic(root_grp,"plik_cmbonly",1,hascl,2508,30)
  lkl_grp["cmbonly_version"]=pars.int(default=18).version
  php.add_external_data(osp.realpath(pars.data),lkl_grp,tar=bool(pars.int(default=1).include))
  #assert os.system("cd %s;tar cvf data.tar *"%dr)==0
  #f=open(osp.join(dr,"data.tar"),"r")
  #dts = f.read()
  #f.close()
  #lkl_grp.create_dataset("external_data",data=nm.fromstring(dts,dtype=nm.uint8))
  hf.close()


    #if hasattr(clik,"clik"):
  #  res = php.add_selfcheck(pars.res_object,mcl)
  #  print "lkl for init cl %g"%res
  #
  #if "cl_save" in pars:
  #  f=open(pars.cl_save,"w")
  #  for ci in mcl:
  #    print >>f,ci
  #  f.close()

import sys
if __name__=="__main__":
  main(sys.argv)