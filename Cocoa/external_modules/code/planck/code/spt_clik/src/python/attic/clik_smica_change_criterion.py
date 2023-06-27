#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import numpy as nm
import clik.parobject as php
import clik
import clik.hpy as hpy
import clik.smicahlp as smh


  
def main(argv):
  pars = clik.miniparse(argv[1])
  inhf = hpy.File(pars.input_object)
  assert inhf["clik/lkl_0"].attrs["lkl_type"].lower() == "smica"
  inhf.close()

  hpy.copyfile(pars.input_object,pars.res_object)
  outhf = hpy.File(pars.res_object,"r+")
  lkl_grp = outhf["clik/lkl_0"]

  crit = pars.criterion
  extra={}
  if crit=="quad":
    if "rq" in pars:
      extra["fid"]=nm.fromfile(pars.rq)
    if "mask" in pars:

      extra["mask"] = nm.asarray(nm.reshape(nm.loadtxt(pars.mask),(lkl_grp.attrs["m_channel_T"],lkl_grp.attrs["m_channel_T"])),nm.int)
      print(extra["mask"])
  smh.set_criterion(lkl_grp,crit,**extra)
  outhf.close()  
    
import sys
if __name__=="__main__":
  main(sys.argv)