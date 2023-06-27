#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import clik.parobject as php
import clik
import os.path as osp
import os
import clik.hpy as hpy

def main(argv):
  pars = clik.miniparse(argv[1])
  
  inhf = hpy.File(pars.input_object)
  install_path = inhf["clik/lkl_%d"%pars.int(default=0).lkl_id].attrs["external_dir"]
  
  hpy.copyfile(pars.input_object,pars.res_object)
  outhf = hpy.File(pars.res_object,"r+")
  lkl_grp = outhf["clik/lkl_%d"%pars.int(default=0).lkl_id]
  if hpy.is_h5py_object(inhf):
    assert os.system("cd %s;tar cvf data.tar *"%install_path)==0
    f=open(osp.join(install_path,"data.tar"),"r")
    dts = f.read()
    f.close()
    del lkl_grp.attrs["external_dir"]
    lkl_grp.create_dataset("external_data",data=nm.fromstring(dts,dtype=nm.uint8))
  else:
    local_path = osp.join(lkl_grp._name,"_external")
    os.mkdir(local_path)
    shutil.copytree(install_path,local_path)
    lkl_grp["external_dir"]="."    

  outhf.close()
  

import sys
if __name__=="__main__":
  main(sys.argv)