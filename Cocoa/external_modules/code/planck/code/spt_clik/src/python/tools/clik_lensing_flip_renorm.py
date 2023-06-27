#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import numpy as nm
import clik.parobject as php
import clik
import clik.hpy as hpy
import os.path as osp

def main(argv):
  lkl = clik.clik_lensing(argv[1])
  del(lkl)

  if hpy.cldf.is_cldf(argv[1]):
    hpy.copyfile(argv[1],argv[2])
    root = hpy.File(argv[2],"r+")
    root["clik_lensing/renorm"] = int( not bool(root["clik_lensing/renorm"]))
    root.close()
  else:
    rr = open(argv[1]).read()
    if ("# format: mono") in rr:
      itype = 1
    elif ("# format: qecl") in rr:
      itype = 2
    else:
      itype = 3
    root = hpy.File(argv[2],"w")
    root.create_group("clik_lensing")
    root["clik_lensing/renorm"]=0
    root["clik_lensing/ren1"]=1
    root["clik_lensing/itype"]=itype
    hpy.copyfile(argv[1],argv[2]+"/clik_lensing/"+osp.basename(argv[1]))
    root["clik_lensing/filename"]=osp.basename(argv[1])
    
import sys
if __name__=="__main__":
  main(sys.argv)