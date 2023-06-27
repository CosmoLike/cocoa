#! PYTHONEXE

import sys
sys.path = ["REPLACEPATH"]+sys.path

import numpy as nm
import numpy.random as ra
import numpy.linalg as la
import clik.parobject as php
import clik
import re
import clik.hpy as h5py
import clik.smicahlp as smh

    
    
def do_nlbl(fwhm,ns,Acmb,lmin,lmax,bins=None):
  # fwhm is in arcmin
  sigma2 = (fwhm/60./180.*nm.pi)**2 / (8.*nm.log(2.))
  ell = nm.arange(lmin,lmax+1)
  bl = nm.exp(-ell*(ell+1)*sigma2/2.)
  nl = ns * Acmb**2/(bl**2)*4*nm.pi
  if bins is not None:
    return nm.dot(bins,nl)
  return nl

def uniformbin (lstart, lstop, deltal):
  # from maude's spherelib
  L    = lstop - lstart + 1
  nbin = nm.ceil(L / deltal).astype(nm.int32)
  if nbin==0:
      nbin = 1
  bin  = nm.zeros((nbin, 2))
  bin[0,0] = lstart
  for q in range(1,nbin):
      bin[q  , 0] = bin[q-1, 0] + deltal
      bin[q-1, 1] = bin[q  , 0] - 1
  bin[-1,1] = lstop
  return bin
  
def bin_cl (cl, bin):
  # from maude's spherelib
  nbin = nm.shape(bin)[0]
  cq   = nm.zeros((nbin, 1))
  for q in range(nbin):
    lmin   = bin[q,0]
    lmax   = bin[q,1]
    rg     = nm.arange(lmin,lmax+1).astype(nm.int32)
    lsize = float(nm.size(rg))
    cq[q] = nm.sum(cl[rg])/lsize
  return nm.squeeze(cq)
        
def smicaBin(lmin,lmax,szbin):
  
  mbins = uniformbin(lmin,lmax,szbin)
  bins = nm.zeros((mbins.shape[0],lmax+1-lmin),dtype=nm.double)
  for l in range(lmax+1):
    clin = nm.zeros(lmax+1)
    clin[l] = 1
    bins[:,l-lmin] = bin_cl(clin,mbins)
  return bins

def checkhcl(hcl):
  hascl = nm.array([int(ss) for ss in hcl])
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
  return hascl,ncl
  
def generate_rq_hat(pars,rq,wq):
  if "meanfield" in pars and pars.int.meanfield:
    return rq*1.
  
  if "seed" in pars:
    ra.seed(pars.int.seed)
      
  m = rq.shape[1]
  rqh = rq*0.
  for iq in range(wq.shape[0]):
    xq = ra.standard_normal((m,wq[iq]))
    Sq = nm.dot(xq,xq.T) / wq[iq]
    rhalf = la.cholesky(rq[iq])
    rqh[iq] = nm.dot(rhalf,nm.dot(Sq,rhalf.T))
  return rqh

def build_rq_noise(cnq,np,hascl):
  ncnq = cnq
  if hascl[1] and hascl[2]:
    ncnq = nm.concatenate((cnq,cnq[-np:]))
  rq0 = nm.zeros((ncnq.shape[1],ncnq.shape[0],ncnq.shape[0]))
  for iq in range(ncnq.shape[1]):
    rq0[iq] = nm.diag(ncnq[:,iq])
  return rq0
  
def build_rq_cmb(cq,hascl,Acmb,nt,np):
  ne = hascl[1]*np
  nb = hascl[2]*np
  m = nt+ne+nb
  if len(Acmb)!=m:
    raise Exception("bad")
    
  if hascl[1:].sum()==0:
    # mode T pure
    rq = nm.zeros((len(cq),nt,nt))
    AAt = nm.outer(Acmb,Acmb)
    rq = AAt[nm.newaxis,:,:] * cq.flat[:][:,nm.newaxis,nm.newaxis]
    return rq
  At = Acmb[:nt]
  Ae = Acmb[nt:]
  Ab = Ae
  
  A = nm.zeros((m,nm.sum(hascl[:3])))
  rcq = nm.zeros((cq.shape[1],nm.sum(hascl[:3]),nm.sum(hascl[:3])))
  
  i=0
  # diagonal terms
  if hascl[0]:
    A[:nt,0] = At
    rcq[:,0,0] = cq[0]
  if hascl[1]:
    A[nt:,hascl[0]] = Ae
    rcq[:,hascl[0],hascl[0]] = cq[hascl[0]]
  if hascl[2]:
    A[nt:,hascl[:2].sum()] = Ab
    rcq[:,hascl[:2].sum(),hascl[:2].sum()] = cq[hascl[:2].sum()]

  if hascl[3]: #TE
    rcq[:,0,1] = cq[hascl[:3].sum()]
    rcq[:,1,0] = rcq[:,0,1]
  if hascl[4]: #TB
    rcq[:,0,hascl[:2].sum()] = cq[hascl[:4].sum()]
    rcq[:,hascl[:2].sum(),0] = rcq[:,0,hascl[:2].sum()]
  if hascl[5]: #EB
    rcq[:,hascl[:1].sum(),hascl[:2].sum()] = cq[hascl[:5].sum()]
    rcq[:,hascl[:2].sum(),hascl[:1].sum()] = rcq[:,hascl[:1].sum(),hascl[:2].sum()]  

  # do the tensor product
  rq = nm.tensordot(nm.dot(rcq,A.T),A,((1,),(1,)))
  
  return rq
  
def read_cl(pars,ncl,hascl):
  """returns cls as an array of shape (ncl,lmax+1) ell starts at 0"""
  
  #read cls
  # for now an ascii file
  print pars.cl
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

  lmin = 0
  lmax = cl.shape[1]-1

  #cut lmin/lmax if needed
  
  if hasattr(pars,"lmax"):
    lmax = pars.int.lmax
  cl = cl[:,:lmax+1]
  mcl = cl
  if hasattr(pars,"lmin"):
    lmin = pars.int.lmin
    cl = cl[:,lmin:]

  cl = cl.flat[:]
  mcl = mcl.flat[:]
  if "cl_save" in pars:
    f=open(pars.cl_save,"w")
    for ci in mcl:
      print >>f,ci
    f.close()
  return cl,mcl,lmin,lmax
  
def read_Acmb(pars):
  if hasattr(pars,"mixcol"):
    Acmb = pars.float_array.mixcol
  else:
    Acmb = read_somearray(pars.Acmb)
  return Acmb  
  
def read_bins(pars,cl,lmin,lmax,ncl,hascl):
  """read a binfile, returns an array of size (nbins*ncl)*((lmax-lmin+1)*ncl)"""
  hb=False
  if hasattr(pars,"bins"):
    hb = True
    bins = read_somearray(pars.bins)
    bins.shape=(-1,len(cl))
    cq = nm.dot(bins,cl)
    nbins = bins.shape[0]
  elif hasattr(pars,"binsz"):
    hb = True
    bins = smicaBin(lmin,lmax,pars.int.binsz)
    lins = nm.zeros((bins.shape[0]*ncl,bins.shape[1]*ncl))
    for i in range(hascl.sum()):
      lins[i*bins.shape[0]:(i+1)*bins.shape[0],i*bins.shape[1]:(i+1)*bins.shape[1]] = bins
    bins =lins
    bins.shape=(-1,len(cl))
    nbins = bins.shape[0]
    cq = nm.dot(bins,cl)
  else:
    bins = None
    cq = cl
    nbins = 0
    
  cq.shape = (ncl,-1)
  nq = cq.shape[1]
  if hb:
    oins = bins[:nq,:lmax+1-lmin]
  else:
    oins = None
  print bins.shape,nbins
  return hb,nbins,bins,oins,cq,nq
 
def read_somearray(somepath):
  # for now only ascii arrays
  return nm.loadtxt(somepath) 
  
def main(argv):
  pars = clik.miniparse(argv[1])

  #get hascl
  print "read cl"
  hascl,ncl = checkhcl(pars.int_array.has_cl)
  cl,mcl,lmin,lmax = read_cl(pars,ncl,hascl)
  
  print "build binning"
  hb,nbins,bins,oins,cq,nq = read_bins(pars,cl,lmin,lmax,ncl,hascl)

  #build rq_cmb
  print "build model"
  Acmb = read_Acmb(pars)
  rq_cmb = build_rq_cmb(cq,hascl,Acmb,pars.int.nT,pars.int.nP)

  rq = rq_cmb
  m = rq.shape[-1]
  rq0 = 0

  # add rq0 if needed
  h0=False
  if hasattr(pars,"Rq_0"):
    h0=True
    rq0n = read_somearray(pars.Rq_0)
    rq0n.shape = (nq,m,m)
    rq0 = rq0 + rq0n
  
  if hasattr(pars,"Nq") or hasattr(pars,"noise"):
    h0=True
    if hasattr(pars,"Nq"):
      cnq = read_somearray(pars.Nq)
    else :
      cnq = nm.array([do_nlbl(fwhm,ns,A,lmin,lmax,oins) for fwhm,ns,A in zip(pars.float_array.fwhm,pars.float_array.noise,Acmb)])    
  
    cnq.shape=(pars.int.nT+pars.int.nP,nq)
    rq0n = build_rq_noise(cnq,pars.int.nP,hascl)
    rq0 = rq0 + rq0n

  rq = rq+rq0

  # deal with wq
  print "build weights"
  if hasattr(pars,"fsky"):
    wl = 2*nm.arange(lmin,lmax+1)+1.
    if not hb:
      wq = wl
    else:
      ibins = (oins!=0)
      wq = nm.dot(ibins,wl)
  else:
    wq = read_somearray(pars.wq)

  #generate rqhat
  print "generate fake data"
  rqhat = generate_rq_hat(pars,rq,wq)


  #save all !
  root_grp,hf = php.baseCreateParobject(pars.res_object)
  if (not h0):
    rq0 = None

  lkl_grp = smh.base_smica(root_grp,hascl,lmin,lmax,pars.int.nT,pars.int.nP,wq,rqhat,Acmb,rq0,bins)
  
  
  hf.close()
  
    
  print "file created !"
  # reread the likelihood !
  if hasattr(clik,"clik"):
    res = php.add_selfcheck(pars.res_object,mcl)
    print "lkl for init cl %g"%res
  
import sys
if __name__=="__main__":
  main(sys.argv)
