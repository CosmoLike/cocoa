#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import numpy as nm
import clik.cldf as cldf
import os.path as osp
import os

def main(argv):

  base = (argv[1])
  try:
    f = cldf.File(base)
    kk = list(f.keys())
  except IOError as e:
    print("Can't cldf_ls %s"%argv[1])
    return
  kk.sort()
  res = []
  for k in kk:
    if isinstance(f[k],cldf.File):
      res +=[k+"/"]
    else:
      res +=[k]
  
  rows, columns = os.popen('stty size', 'r').read().split()    
  col = int(columns)
  sz = max([len(k)+1 for k in res])
  sz = max(20,sz)
  cc = col/sz
  cnt = 0
  txt = ""
  for k in res:
    txt += k.ljust(sz)
    cnt +=1
    if cnt == cc:
      cnt = 0
      txt +="\n"
  print(txt)

import sys
if __name__=="__main__":
  main(sys.argv)