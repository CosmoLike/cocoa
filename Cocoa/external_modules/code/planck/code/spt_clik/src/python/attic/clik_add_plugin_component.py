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
import shutil
import clik.hpy as h5py
import clik.smicahlp as slp

def main(argv):
  pars = clik.miniparse(argv[1])
  
  inhf = h5py.File(pars.input_object)
  assert inhf["clik/lkl_0"].attrs["lkl_type"].lower() == "smica"
  inhf.close()
  
  shutil.copyfile(pars.input_object,pars.res_object)
  outhf = h5py.File(pars.res_object,"r+")
  lkl_grp = outhf["clik/lkl_0"]

  if "lmax" in lkl_grp.attrs:
    lmax = lkl_grp.attrs["lmax"]
    lmin = 0
    if "lmin" in lkl_grp.attrs:
      lmin = lkl_grp.attrs["lmin"]
    ellh = False
  else:
    ell = lkl_grp.attrs["ell"]
    lmax = nm.max(ell)      
    lmin = nm.min(ell)

  name = pars.name
  dets = pars.int_array.freq
  vpars = list(pars.str_array.parameter)
  defaults = dict(zip(pars.str_array(default="").default_par,pars.str_array(default="").default_value))
  position = -1
  
  if pars.str(default="sky").component_type.lower() == "sky":
    position =1

  data = None
  if pars.str(default="").template_file:
    data = nm.loadtxt(pars.template_file).flat[:]
  color = None
  if pars.str(default="").color_file:
    color = nm.loadtxt(pars.color_file).flat[:]
  if "color" in pars:
    color = nm.array(pars.float_array.color)

  rename = {}
  if "rename_prefix" in pars:
    prefix = pars.rename_prefix
    lp = len(prefix)
    rename = dict([(pp,pp[lp:]) for pp in vpars + defaults.keys() if pp.startwith(prefix)])
  if "rename_from" in pars:
    rename = dict([vv for vv in zip(pars.str_array.rename_from, pars.str_array.rename_to)])

  slp.add_parametric_component(lkl_grp,name,dets,vpars,lmin,lmax,defaults,color=color, rename = rename, data=data,position=position)

  php.remove_selfcheck(root_grp = outhf["clik"])
  
  outhf.close()

import sys
if __name__=="__main__":
  main(sys.argv)