#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import numpy as nm
import clik.cldf as cldf
import os.path as osp

def main(argv):

  base = osp.dirname(argv[1])
  f = cldf.File(base)
  print(f[osp.basename(argv[1])])

import sys
if __name__=="__main__":
  main(sys.argv)