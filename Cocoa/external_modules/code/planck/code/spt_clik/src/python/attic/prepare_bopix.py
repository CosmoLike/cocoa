#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import numpy as nm
import numpy.random as ra
import numpy.linalg as la
import clik.parobject as php
import clik
import re
import os.path as osp
import os

def exists(partxt,data_dir,ENTRY,nofail=False):
  rr = re.findall(ENTRY+"\s*=\s*(.+)",partxt)
  if rr:
    ent = rr[0]
    #print rr
    return osp.exists(osp.join(data_dir,ent))
  assert bool(rr) or nofail
  
def main(argv):
  pars = clik.miniparse(argv[1])

  bopix_dir = osp.realpath(pars.bopix_data)
  fpartxt=open(osp.join(bopix_dir,"input_file.txt"))
  partxt = "\n".join([l.split("#")[0].strip() for l in fpartxt if l.split("#")[0
].strip()])
  #print partxt
  
  # get data directory
  data_dir = osp.join(bopix_dir,re.findall("BOPIX_DATA_DIR\s*=\s*(.+)",partxt)[0
])
  assert osp.exists(data_dir)
  
  for fi in ["BOPIX_CL_FILENAME","COV_FILE","WINDOW_FILE"]:
    assert exists(partxt,data_dir,fi)
  
  if not exists(partxt,data_dir,"MAP_FILE",True):
    assert exists(partxt,data_dir,"MASKED_MAP_P")
    assert exists(partxt,data_dir,"MASKED_MAP_T")

  if not exists(partxt,data_dir,"MASKFILE",True):
    assert exists(partxt,data_dir,"MASKFILE_P")
    assert exists(partxt,data_dir,"MASKFILE_T")
      
  
  lmax = int(re.findall("BOPIX_CL_LMAX\s*=\s*(.+)",partxt)[0])
  lmin = 0

  root_grp,hf = php.baseCreateParobject(pars.res_object)
  hascl = nm.array([1,1,0,1,0,0],dtype=nm.int)
  
  lkl_grp = php.add_lkl_generic(root_grp,"bopix",1,hascl,lmax,lmin)
  if "external" in pars and pars.int.external == 1:
    lkl_grp.attrs["external_dir"] = bopix_dir
  else:
    assert os.system("cd %s;tar cvf /tmp/bopix_data.tar *"%bopix_dir)==0
    dts = open("/tmp/bopix_data.tar").read()
    assert os.system("rm /tmp/bopix_data.tar")==0
  
    lkl_grp.create_dataset("external_data",data=nm.frombuffer(dts,dtype=nm.uint8))
  
  php.add_pid(lkl_grp,pars.str(default="").pid)
  
  hf.close()
  
  cls = nm.loadtxt(osp.join(data_dir,re.findall("BOPIX_CL_FILENAME\s*=\s*(.+)",partxt)[0]))
  mcl = nm.zeros((3,lmax+1))
  llp1s2pi = nm.arange(2,lmax+1)*nm.arange(3,lmax+2)/2./nm.pi
  #mcl[0,2:] = cls[:lmax-1,1]/llp1s2pi
  #mcl[1,2:] = cls[:lmax-1,3]/llp1s2pi
  #mcl[2,2:] = cls[:lmax-1,2]/llp1s2pi
  mcl[0,2:] = cls[:lmax-1,1]/llp1s2pi
  mcl[1,2:] = cls[:lmax-1,2]/llp1s2pi
  mcl[2,2:] = cls[:lmax-1,4]/llp1s2pi
  #print "TT wow",mcl[0,2:6]
  #print "EE wow",mcl[1,2:6]
  #print "TE wow",mcl[2,2:6]
  
  if hasattr(clik,"clik"):
    res = php.add_selfcheck(pars.res_object,mcl.flat[:])
    print "lkl for init cl %g"%res
  
  if "cl_save" in pars:
    f=open(pars.cl_save,"w")
    for ci in mcl.flat:
      print >>f,ci
    f.close()
  
import sys
if __name__=="__main__":
  main(sys.argv)
