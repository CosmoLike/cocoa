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
  parfile = pars.comm_lowl_parfile.strip()
  clw_pars = clik.miniparse(parfile)
  lmin = clw_pars.int.LMIN
  lmax = clw_pars.int.LMAX
  datafile = clw_pars.DATAFILE.strip()
  fid_cl = clw_pars.FIDUCIAL_CL_FILE.strip()

  root_grp,hf = php.baseCreateParobject(pars.res_object)
  hascl = [1]*6
  hascl = nm.array(hascl,dtype=nm.int)
  
  lkl_grp = php.add_lkl_generic(root_grp,"comm_lowl",1,hascl,lmax,lmin)
  
  php.add_pid(lkl_grp,pars.str(default="").pid)

  import tempfile
  dr = tempfile.mkdtemp()
  print(dr)
  import shutil
  shutil.copy(parfile,dr+"/comm_lowl.par")
  shutil.copy(datafile,dr+"/"+datafile)
  shutil.copy(fid_cl,dr+"/"+fid_cl)
  
  php.add_external_data(dr,lkl_grp,tar=True)

  #assert os.system("cd %s;tar cvf data.tar *"%dr)==0
  #f=open(osp.join(dr,"data.tar"),"r")
  #dts = f.read()
  #f.close()
  #lkl_grp.create_dataset("external_data",data=nm.fromstring(dts,dtype=nm.uint8))
  hf.close()

  shutil.rmtree(dr)

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