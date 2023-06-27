#! PYTHONEXE

import numpy as nm
import healpy as hp

import sys
if __name__ == "__main__":
  sys.path = ["REPLACEPATH"]+sys.path

def mode(l,m,nside,field,ordering="ring",all=False):
  if m==0:
    value = 1
  else:
    value = (1+nm.sign(m)*1j)/2.
  nmode = ((l+1)*(l+2))/2
  Ind = hp.Alm.getidx(l,l,abs(m))
  dat = nm.zeros((nmode,),dtype=nm.complex)
  dat[Ind] = value
  if field.lower()=="t":
    res = hp.alm2map(dat,nside,l)
  else:
    datT = dat*0
    if field.lower()=="e":
      datE = dat
      datB = dat*0
    else:
      datB = dat
      datE = dat*0
    rr = hp.alm2map((datT,datE,datB),nside,l)
    res = nm.array(rr[1:]).flat[:]
  if ordering.lower() !="ring" :
    #fuck reorder !
    if field.lower()=="t":
      mres = res[hp.nest2ring(nside,list(range(12*nside**2)))]
      res = mres
    else:
      mres = res[nm.concatenate((hp.nest2ring(nside,list(range(12*nside**2))),hp.nest2ring(nside,list(range(12*nside**2)))+12*nside**2))]
      res = mres
      
  if all:
    return res,dat
    
  return res
  
def cYLM(nside=16,field="t",lmax=-1,ordering="ring",ls=None):
  if field.lower() not in ["t","e","b"]:
    raise ValueError("field must be one of T, E and B")
  
  field = field.lower()
  npix = 12*nside**2
  
  if field != "t":
    npix*=2
    
  if lmax==-1:
    lmax = 3*nside
  
  lmax=int(lmax)
  
  if ls == None:
    ls = nm.arange(lmax+1)
  
  nmode = nm.sum(2*ls+1)
  
  res = nm.zeros((npix,nmode),dtype=nm.float64)
  cnt=0
  for l in ls:
    for m in range(-l,l+1):
      rr = mode(l,m,nside,field,ordering)
      res[:,cnt]= rr
      cnt+=1
    
  return res
  
def applyMask(ylm,mask):
  if ylm.shape[0]==len(mask):
    return nm.compress(mask==1,ylm,axis=0)
  elif ylm.shape[0]==2*len(mask):
    return nm.compress(nm.tile(mask,2)==1,ylm,axis=0)
    
def applyls(ylm,ls):
  rls=[]
  cnt=0
  for l in range(max(ls)+1):
    if l in ls:
      rls+=(nm.arange(2*l+1)+cnt).tolist()
    cnt+=2*l+1
  return ylm[:,rls]

import numpy.linalg as la

def svdReduce(mylm,level=1,ls=None,nside=-1,all=False):
  if nside==-1:
    pside=nm.sqrt(mylm.shape[0]/12.)
    nside=2**nm.ceil(nm.log(pside)/nm.log(2))
  dlevel = (1-level)*nm.sqrt(3*nside**2/nm.pi)
  if ls!=None:
    mylm = applyls(mylm,ls)
  if level==1 and not all:
    q,r=la.qr(mylm)
    return q,r
  u,d,vt = la.svd(mylm,full_matrices=False)
  if level==1:
    complev=nm.ones(u.shape[1])
  else:
    complev=d>=dlevel
  u=nm.compress(complev,u,axis=1)
  vt=nm.compress(complev,vt,axis=0)
  df=nm.compress(complev,d,axis=0)
  if all:
    return u,df[:,nm.newaxis]*vt,d,df
  return u,df[:,nm.newaxis]*vt
  
def reducedYLM(ylm, mask, level):
  mylm=applyMask(ylm,mask)
  return svdReduce(mylm,level)
  
  
class YLM:
  head="YLM"
  def _normalize(self,nside=0,masks=None,ls=None):
    if masks == None:
      amask = nm.ones((2,nside**2*12),dtype=nm.uint8)
    elif len(masks)!=2:
      amask = nm.array((masks,masks)).astype(nm.uint8)
    else:
      amask = nm.array(masks).astype(nm.uint8)
    if nside == 0:
      nside = nm.sqrt(amask.shape[0]/12)
    if ls == None :
      ls = tuple(nm.arange(3*nside+1))
    return nside,amask,ls

  def _checksame(self,nside,masks,sigmask,ls,ordering,has_cl,ylm):
      same = (nside == 0 or ylm.nside==nside)
      same = same and (tuple(has_cl) == ylm.has_cl)
      same = same and (nside == 0 or ls ==None or tuple(ls) == ylm.ls)
      same = same and (ordering.lower() == ylm.ordering)
      same = same and (masks == None or sigmask == ylm.sigmask)
      return same
  
  def _store(self,nside,ls,ordering,sigmask,has_cl):
    self.nside = nside
    self.ls = tuple(ls)
    self.ordering = ordering.lower()
    self.sigmask = tuple(sigmask)
    self.has_cl = tuple(has_cl)

  def __init__(self,nside=0,masks=None,ls=None,has_cl=[1,0,0,0,0,0],ordering="RING",ylm=None,fromfile=""):

    nside,amasks,ls = self._normalize(nside,masks,ls)  
    sigmask = tuple(nm.sum(amasks,axis=1))    
    if fromfile:
      self.load(fromfile)
      if not self._checksame(nside,masks,sigmask,ls,ordering,has_cl,self):
        raise Exception("Not the same ylm")
      return

    
    if ylm!=None:
      if not self._checksame(nside,masks,sigmask,ls,ordering,has_cl,ylm):
        raise Exception("Not the same ylm")
      self.ylm=ylm.ylm
      self._store(ylm.nside,ylm.ls,ylm.ordering,ylm.sigmask,ylm.has_cl)
      return
    ylmT=nm.zeros((0,0))
    ylmE=nm.zeros((0,0))
    ylmB=nm.zeros((0,0))
    if has_cl[0] or has_cl[3] or has_cl[4]:
      ylmT=applyMask(cYLM(nside,"t",max(ls),ordering,ls),amasks[0])
    if has_cl[1] or has_cl[3] or has_cl[5]:
      ylmE=applyMask(cYLM(nside,"e",max(ls),ordering,ls),amasks[1])
      if ylmE.shape[0]==0:
        ylmE=nm.zeros((0,0))
    if has_cl[2] or has_cl[4] or has_cl[5]:
      ylmB=applyMask(cYLM(nside,"b",max(ls),ordering,ls),amasks[1])
      if ylmB.shape[0]==0:
        ylmB=nm.zeros((0,0))

    ylm=nm.zeros((ylmT.shape[0]+max(ylmE.shape[0],ylmB.shape[0]),ylmT.shape[1]+ylmE.shape[1]+ylmB.shape[1]))
    #print ylm.shape,ylmT.shape,ylmE.shape,ylmB.shape
    ylm[:ylmT.shape[0],:ylmT.shape[1]]=ylmT
    if  ylmE.shape[0]!=0:
      ylm[ylmT.shape[0]:,ylmT.shape[1]:ylmT.shape[1]+ylmE.shape[1]]=ylmE
    if  ylmB.shape[0]!=0:
      ylm[ylmT.shape[0]:,ylmT.shape[1]+ylmE.shape[1]:]=ylmB
    
    self.ylm=ylm
    self._store(nside,ls,ordering,sigmask,has_cl)
    
    return
  
  def _saveinfo(self,f):
    f.write(self.head)
    f.write(self.ordering[0])
    f.write(nm.array(self.nside,dtype=nm.int32).tostring())
    f.write(nm.array(len(self.ls),dtype=nm.int32).tostring())
    f.write(nm.array(self.ls,dtype=nm.int32).tostring())
    f.write(nm.array(self.sigmask,dtype=nm.int32).tostring())
    f.write(nm.array(self.has_cl,dtype=nm.int32).tostring())
    
  def _saveArray(self,f,arr):
    f.write(nm.array(arr.shape,dtype=nm.int32).tostring())
    f.write(nm.array(arr,dtype=nm.float64).tostring())
    
  def save(self,fi):
    f=open(fi,"w")
    self._saveinfo(f)
    self._saveArray(f,self.ylm)
    f.close()
    
  def _loadinfo(self,f):
    if f.read(len(self.head))!=self.head:
      raise Exception("argl 1")
    self.ordering={"r":"ring","n":"nested"}[f.read(1)]
    self.nside = nm.frombuffer(f.read(4),dtype=nm.int32)[0]
    nls = nm.frombuffer(f.read(4),dtype=nm.int32)[0]
    self.ls = tuple(nm.frombuffer(f.read(4*nls),dtype=nm.int32))
    self.sigmask = tuple(nm.frombuffer(f.read(4*2),dtype=nm.int32))
    self.has_cl = tuple(nm.frombuffer(f.read(4*6),dtype=nm.int32))

  def _loadArray(self,f):
    shape = tuple(nm.frombuffer(f.read(4*2),dtype=nm.int32))
    arr = (nm.frombuffer(f.read(8*shape[0]*shape[1]),dtype=nm.float64))
    arr.shape=shape
    return arr
    
  def load(self,fi):
    f=open(fi,"r")
    self._loadinfo(f)
    self.ylm = self._loadArray(f)
    f.close()
    
class UG(YLM):
  head="UG"
  def __init__(self,nside=0,masks=None,ls=None,has_cl=[1,0,0,0,0,0],ordering="RING",level=1,ylm=None,ylm_fromfile="",ug=None,fromfile=""):
    nside,amasks,ls = self._normalize(nside,masks,ls)  
    sigmask = tuple(nm.sum(amasks,axis=1))    
  
    if fromfile:
      self.load(fromfile)
      if not self._checksame(nside,masks,sigmask,ls,ordering,has_cl,level,self):
        raise Exception("Not the same UG")
      return

    if ug != None:
      if not self._checksame(nside,masks,sigmask,ls,ordering,has_cl,level,ug):
        raise Exception("Not the same UG")
      self.U=ug.U
      self.G=ug.G
      self._store(ug.nside,ug.ls,ug.ordering,ug.sigmask,ug.has_cl,ug.level)
      return
      
    ylm = YLM(nside,masks,ls,has_cl,ordering,ylm,ylm_fromfile)
    U,G = svdReduce(ylm.ylm,level)
    self.U = U
    self.G = G
    self._store(nside,ls,ordering,sigmask,has_cl,level)
    
  def _store(self,nside,ls,ordering,sigmask,has_cl,level):
    self.level=level
    YLM._store(self,nside,ls,ordering,sigmask,has_cl)
  def _checksame(self,nside,masks,sigmask,ls,ordering,has_cl,level,ug):
    same = YLM._checksame(self,nside,masks,sigmask,ls,ordering,has_cl,ug)
    same = same and (level == ug.level)
    return same 
  
  def _saveinfo(self,f):
    YLM._saveinfo(self,f)
    f.write(nm.array(self.level,dtype=nm.float64).tostring())
  
  def _loadinfo(self,f):
    YLM._loadinfo(self,f)
    self.level = nm.frombuffer(f.read(8),dtype=nm.float64)[0]
    
  def save(self,fi):
    f=open(fi,"w")
    self._saveinfo(f)
    self._saveArray(f,self.U)
    self._saveArray(f,self.G)
    f.close()
  def load(self,fi):
    f=open(fi,"r")
    self._loadinfo(f)
    self.U = self._loadArray(f)
    self.G = self._loadArray(f)
    f.close()

def prepUG(nside,mask,ls,has_cl,svdCut=None):
  ylms = YLM(nside,mask,ls,has_cl,"ring") 
  if svdCut!=None:
    UG = UG(nside,mask,ls,has_cl,"ring",svdCut,ylms) 
    return UG.U,UG.G
  return ylms.ylm,None

def prepdats(masks, maps,noise,ell,cl,has_cl,svdCut=None):
  import clik.lkl
  if len(maps)==3:
    nside = int(nm.sqrt(len(maps[0])/12))
  else:
    nside = int(nm.sqrt(len(maps)/12))

  U,G = prepUG(nside,masks,ell,has_cl,svdCut)
  
  if masks == None:
    masks = nm.ones(12*nside*nside,dtype=nm.uint8)
  return clik.lkl.powly_javel(masks, maps,noise,ell,cl,has_cl,U,G)
  
  
def hdf_lowl(hf,cl_fid,a_bar,H):
  hf.attrs["neff"] = len(a_bar)
  hf.create_dataset('cl_fid', data=cl_fid)  
  hf.create_dataset('a_bar', data=a_bar)  
  hf.create_dataset('H', data=H.flat[:])
  

def main(argv):
  import clik
  pars = clik.miniparse(argv[1])
  
  hascl = nm.array([int(ss) for ss in pars.str_array.has_cl])
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
  
  if nm.sum(hascl[1:])==0:
    rmaps = hp.read_map(pars.str.mapfile)
    print(rmaps[0])
  else:
    rmaps = [hp.read_map(pars.str.mapfile,i) for i in range(3)]
    print(rmaps[0][0],rmaps[1][0],rmaps[2][0])
    
  cl = nm.loadtxt(pars.str.clfile)
  wl = nm.loadtxt(pars.str.beamfile)
  cl.shape=[ncl,-1]
  
  allcl = nm.zeros((6,cl.shape[1]),dtype=nm.double)
  j=0
  for i in range(6):
    if hascl[i]:
      allcl[i] = cl[j]*wl
      j+=1
  clw = allcl.flat[:]
  print(clw)
  res = prepdats(None,rmaps,pars.float.noise,nm.arange(pars.int.lmin,pars.int.lmax+1),clw,hascl)
  
  return 0
if __name__ == "__main__":
  main(sys.argv)

