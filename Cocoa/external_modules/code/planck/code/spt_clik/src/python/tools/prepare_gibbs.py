#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import numpy as nm
import clik.parobject as php
import clik
import os.path as osp

def main(argv):
  pars = clik.miniparse(argv[1])
  if "version" in pars:
    if "gauss" in pars.str.version.lower():
      main_v2(argv)
      return
  main_v1(argv)

def main_v1(argv):
  pars = clik.miniparse(argv[1])
  #test_cl = nm.loadtxt(osp.join(pars.wmap_data,"data/v2a1s_best_lcdm_6000.txt"))
  
  #mcl = nm.zeros((4,1201),dtype=nm.double)
  #llp1s2pi = nm.arange(1201)*nm.arange(1,1202)/2./nm.pi
  #mcl[:,2:] = (test_cl[:1201-2,1:].T)/llp1s2pi[2:]
  
  lmin = pars.int(default=2).lmin
  lmax = pars.int(default=32).lmax
  firstchain = pars.int.firstchain
  lastchain = pars.int.lastchain
  firstsample = pars.int.firstsample
  lastsample = pars.int.lastsample
  step = pars.int.step
  
  root_grp,hf = php.baseCreateParobject(pars.res_object)
  hascl = [0]*6
  hascl[0] = 1
  hascl = nm.array(hascl,dtype=nm.int)
  
  lkl_grp = php.add_lkl_generic(root_grp,"gibbs",1,hascl,lmax,lmin)
  
  lkl_grp.attrs["firstsample"] = firstsample
  lkl_grp.attrs["lastsample"] = lastsample
  lkl_grp.attrs["firstchain"] = firstchain
  lkl_grp.attrs["lastchain"] = lastchain
  lkl_grp.attrs["step"] = step
  if "approx_chi2" in pars:
    lkl_grp.attrs["approx_chi2"] = pars.int.approx_chi2
  if "ltrans" in pars:
    lkl_grp.attrs["ltrans"] = pars.int.ltrans
    print(pars.ltrans)
  php.add_pid(lkl_grp,pars.str(default="").pid)

  import tempfile
  dr = tempfile.mkdtemp()
  import os
  os.mkdir(dr+"/data")
  print(dr)
  import shutil
  shutil.copy(pars.sigma_file.strip(),dr+"/data/sigma.fits")
  shutil.copy(pars.cl_file.strip(),dr+"/data/cl.dat")

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

def main_v2(argv):
  pars = clik.miniparse(argv[1])
  
  lmin = pars.int(default=2).lmin
  lmax = pars.int(default=32).lmax
  delta_l = pars.int.delta_l
  
  root_grp,hf = php.baseCreateParobject(pars.res_object)
  hascl = [0]*6
  hascl[0] = 1
  hascl = nm.array(hascl,dtype=nm.int)
  
  lkl_grp = php.add_lkl_generic(root_grp,"gibbs_gauss",1,hascl,lmax,lmin)
  
  if pars.str.version.lower()=="gauss3":
    lkl_grp.attrs["version"] = 3
  lkl_grp.attrs["delta_l"] = delta_l
  php.add_pid(lkl_grp,pars.str(default="").pid)

  import tempfile
  dr = tempfile.mkdtemp()
  import os
  print(dr)
  import shutil
  shutil.copy(pars.sigma_file.strip(),dr+"/sigma.fits")
  
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