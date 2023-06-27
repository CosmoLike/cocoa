#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import numpy as nm
import clik.parobject as php
import clik
import clik.hpy as hpy


def read_gcal_data(pars,lkl_grp):
  # returns the dtemplate data
  datacal = pars.str_array.datacal
  ngcal = pars.float_array.ngcal
  lmin = lkl_grp.attrs["lmin"]
  lmax = lkl_grp.attrs["lmax"]
  if len(datacal) == 1:
    dat = nm.loadtxt(datacal[0]).flat[:]
  else:
    assert len(datacal)==len(ngcal)
    dat = ()
    i=0
    lmax_tpl = pars.int_array(default=[-1]).lmax_tpl
    if len(lmax_tpl)==1:
      lmax_tpl = list(lmax_tpl)*len(ngcal)
    assert len(ngcal)==len(lmax_tpl)
    for ff,lm in zip(datacal,lmax_tpl):
      assert lm>=lmax
      idat = nm.loadtxt(ff).flat[:]
      if lm!=-1:
        idat.shape = (-1,lm+1)
        idat = idat[:ngcal[i],lmin:lmax+1]
        idat = idat.flat[:]
      dat = nm.concatenate((dat,idat))  
  return dat
  
def main(argv):
  pars = clik.miniparse(argv[1])
  inhf = hpy.File(pars.input_object)
  inhf.close()
  assert inhf["clik/lkl_0"].attrs["lkl_type"].lower() == "smica"
  hpy.copyfile(pars.input_object,pars.res_object)
  outhf = hpy.File(pars.res_object,"r+")
  lkl_grp = outhf["clik/lkl_0"]
  nc = lkl_grp.attrs["n_component"]
  position = pars.int(default=nc).position
  assert position <=nc
  for ic in range(nc,position,-1):
    lkl_grp.copy("component_%d"%ic-1,"component_%d"%ic)
    del lkl_grp["component_%d"%ic]
  
  agrp = outhf.create_group("clik/lkl_0/component_%d"%(position))
  nc+=1
  lkl_grp.attrs["n_component"] = nc

  ngcal = pars.float_array.ngcal
  
  dat = read_gcal_data(pars,lkl_grp)

  agrp.attrs["ngcal"] = ngcal
  agrp.create_dataset("gcaltpl",data=dat.flat[:])

  if "name" in pars:
    name = pars.str_array.name
    assert len(name) == nm.sum(ngcal)
    agrp.attrs["name"] = php.pack256(*name)
  
  agrp.attrs["binned"] = int(pars.int(default=0).binned!=0)
  outhf.close()  
    
import sys
if __name__=="__main__":
  main(sys.argv)