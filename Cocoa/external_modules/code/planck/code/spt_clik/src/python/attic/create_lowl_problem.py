#! PYTHONEXE

import sys
sys.path = ["REPLACEPATH"]+sys.path

import healpy as hp
import numpy as nm
import re
import numpy.random as ra
import clik


def read_somearray(somepath):
  # for now only ascii arrays
  return nm.loadtxt(somepath)
      
def main(argv):
  pars = clik.miniparse(argv[1])
  
  hascl = nm.array([int(ss) for ss in pars.int_array.has_cl])
  if hascl[3]:
    if not (hascl[0] and hascl[1]):
      raise Exception("bad 1")
  if hascl[4]:
    if not (hascl[0] and hascl[2]):
      raise Exception("bad 2")
  if hascl[5]:
    if not (hascl[1] and hascl[2]):
      raise Exception("bad 3")
  ncl = nm.sum(hascl)
  if ncl == 0:
    raise Exception("argl")

  
  cl = read_somearray(pars.cl)
  if len(cl.shape)==2:
    # could I be dealing with a camb cl ?
    if nm.sum(cl[0:10,0]-nm.arange(2,12) )==0:
      # that's it !
      ell = cl[:,0]
      llp1 = ell*(ell+1)/2./nm.pi
      print "I think this is a cl from WMAP !"
      mcl = []
      for i in range(4):
        if hascl[i]:
          mcl+=[nm.concatenate(([0.,0.],cl[:,i+1]/llp1))]
      for i in range(4,6):
        if hascl[i]:
          mcl+=[0*nm.concatenate(([0.,0.],cl[:,0]))]
      cl = nm.concatenate(mcl)

        
  cl.shape=(ncl,-1)
  
  lmax = cl.shape[-1]-1
  if hasattr(pars,"lmax"):
    lmax = pars.int.lmax
    
  cl = cl[:,:lmax+1]
  ell = nm.arange(lmax+1)
  
  nside = pars.int.nside
  
  if pars.str.beamtype.lower() == "gauss":
    fwhm_deg = pars.float.fwhm
    sigma = fwhm_deg / (2.*nm.sqrt(2.*nm.log(2.)))*(nm.pi/180.)
    beam = nm.exp(-0.5*ell*(ell+1)*sigma**2)

  if pars.str.beamtype.lower() == "cos":
    c_start = pars.int.cos_start
    c_end = pars.int.cos_end
    beam = nm.ones((lmax+1))
    c_span = -c_start+c_end+1
    beam[c_start:c_end+1] = nm.cos(nm.linspace(0,nm.pi/2,c_span))
    beam[c_end+1:lmax+1] = 0

  hhcl = cl * beam[nm.newaxis,:]**2
  
  if nm.sum(hascl[1:])==0:
    hpcl = hhcl[0]
  else:
    hpcl = nm.zeros((4,lmax+1))
    hpcl[0] = hhcl[0]
    hpcl[2] = hhcl[1]
    off = 0
    if hascl[2]:
      hpcl[3] = hhcl[2]
      off = 1
    if hascl[3]:
      hpcl[1] = hhcl[off+2]
  
  hmap = nm.array(hp.synfast(hpcl,nside,lmax))
  
  if pars.str.noiseunit.strip() == "muK^2.pixel":
    noise = pars.float.noise
  else:
    noise = pars.float.noise / (4*nm.pi/(12*nside*nside))

  print "noise per pixel",noise
  
  nmap = ra.standard_normal(size=hmap.shape)*nm.sqrt(noise)
  
  rmap = hmap + nmap
  
  print rmap.shape
  hp.write_map(pars.str.mapfile,rmap)

  f=open(pars.str.beamfile,"w")
  for v in beam:
    print >>f,v**2
  f.close()

  f=open(pars.str.clbffile,"w")
  for v in cl.flat[:]:
    print >>f,v
  f.close()
  
if __name__=="__main__":
  import sys
  print 1
  main(sys.argv)