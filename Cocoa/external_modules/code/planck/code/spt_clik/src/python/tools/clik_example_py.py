#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import numpy as nm
import clik

def main(argv):
  if clik.try_lensing(argv[1]):
    main_lensing(argv)
    return 
  main_CMB(argv)

def main_CMB(argv):
  from optparse import OptionParser
  parser = OptionParser()
  parser.add_option("-p", "--par", dest="options",
                  help="add option 'TOTO = TATA'",action="append")
  opt,argv = parser.parse_args(argv)
  options = {}
  if opt.options is not None:
    for kv in opt.options:
      kvs = kv.split("=")
      if len(kvs)!=2:
        continue
      k = kvs[0].strip()
      v = kvs[1].strip()
      options[k]=v
  lklfile = argv[1]
  #print(options)
  lkl = clik.clik(lklfile,**options)
  for clfile in argv[2:]:
    cls = nm.loadtxt(clfile)
    nres = lkl(cls.flat[:])
    print(nres)

def main_lensing(argv):
  lklfile = argv[1]
  lkl = clik.clik_lensing(lklfile)
  for clfile in argv[2:]:
    cls = nm.loadtxt(clfile)
    nres = lkl(cls.flat[:])
    print(nres)

if __name__=="__main__":
  main(sys.argv)