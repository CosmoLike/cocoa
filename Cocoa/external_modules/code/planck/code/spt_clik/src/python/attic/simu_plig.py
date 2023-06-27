#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import sys
import os
import healpy as hp
import numpy as nm
import numpy.random as ra
import numpy.linalg as la
import clik.parobject as php
import clik
import re
import clik.smicahlp as smh
import clik.cldf as cldf

def do_one_freq(cmb_alms, nT,nP,beam,noises,name):
  bm = hp.read_cl(beam)
  cmb = hp.alm2map([hp.almxfl(alm,bm) for alm in cmb_alms],nside=2048)
  for i,ns in enumerate(noises):
    noi = hp.read_map(ns,field=(0,1,2))
    if i<nT:
      mp = cmb[0]+noi[0]/nT
      hp.write_map(name+"_T_%d.fits"%i,mp)
    if i<nP:
      mp = cmb[1]+noi[1]/nP
      hp.write_map(name+"_Q_%d.fits"%i,mp)
      mp = cmb[2]+noi[2]/nP
      hp.write_map(name+"_U_%d.fits"%i,mp)
    
def gen_noise(noise_dir,freq,n):
  return nm.random.permutation([l for l in os.listdir(noise_dir+"/%d"%freq) if "hr" not in l])[:n]

def main(argv):
  pars = clik.miniparse(argv[1])

  cmb_alms = [hp.read_alm(f) for f in pars.str_array.cmb_alm)
]  noise_dir = pars.noise_dir
  res_dir = pars.res_dir

  f=open(res_dir+"/description.txt")
  print >>f,"cmb:",pars.str.cmb_alms

  for beam, freq, nT, nP in zip(pars.str_array.beams, pars.int_array.freq, pars.int_array.nT, pars.int_array.nP):
    noises = gen_noise(noise_dir,freq,max(nT,nP))
    print >>f, "%d beam"%freq,beam
    print >>f, "%d noise"%freq," ".join(noises)
    name = res_dir+"/map_%d"%freq
    do_one_freq(cmb_alms, nT,nP,beam,noises,name)




import sys
if __name__=="__main__":
  main(sys.argv)