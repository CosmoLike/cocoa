import re
import numpy as nm
import os.path as osp

try:
  import pyfits as pf
except ImportError as e:
  # try astropy then
  from astropy.io import fits as pf
  

def scarray(li,scal=False):
  if len(li)==1 and scal:
    return li[0]
  else:
    if type(li[0])!=str:
      return nm.array(li)
    return li

def read_array(fname,dirname):
  fname = lookupfile(fname,dirname)
  
  try:
    pfits = pf.open(fname)
    ii=0
    while pfits[ii].data is None:
      ii+=1
    return pfits[ii].data
  except Exception:
    return nm.loadtxt(fname)
    
class transformme:
  def __init__(self,tfunc,pf,isar=False):
    self.tfunc = tfunc
    self.pf = pf
    self.df = None
    self.isar = isar
    self.scal=False
    
  def __getattr__(self,val):
    val = val.replace("_DOT_",".")
    val = val.replace("_dot_",".")
    if self.isar and val+".file" in self.pf.pf:
      vl = read_array(self.pf.pf[val+".file"],self.pf.localdir).flat[:]
      self.pf._access_list += [val+".file"]
    else:
      try:
        vl = self.pf.pf[val]
        self.pf._access_list += [val]
      except Exception as e:
        if self.df is None:
          raise e
        else:
          vl = self.df
          self.pf._access_list += [val+".file"]
    
    if self.isar:
      if isinstance(vl,str):
        vvl = vl.split()
      else:
        vvl = vl
      return scarray([self.tfunc(v) for v in vvl],self.scal)
    return self.tfunc(vl)
    
  def __call__(self,**kk):
    if "default" in kk:
      self.df = kk["default"]
    if kk.get("scalarize",False):
      self.scal = True
    return self

def getnextline(txtit):
  l = next(txtit)
  l = re.split("#|!",l)[0]
  if len(l.strip())==0:
    return getnextline(txtit)
  return l

def cleantxtit(txtit):
  while(1) :
    yield getnextline(txtit)

def parse_text(txtit):
  ctxtit = cleantxtit(txtit)
  cont = ""
  pf = {}
  for l in ctxtit:
    k = cont
    if not cont:
      k,v = re.findall("(?<!#)((?:\w|\.)+)\s*=\s*(.+)",l)[0]
      v = v.strip()
    else:
      v = l.strip()
      v = pf[k]+ " "+ v
    cont = "" 
    if v[-1]=="&": #continuation
        v = v[:-1]
        cont = k
    pf[k] = v
  return pf

def lookupfile(fi,dirs):
  dirs = list(dirs) + ["."]
  if osp.exists(fi):
    return osp.realpath(fi)
  for d in dirs:
    if osp.exists(osp.join(dd,fi)):
      return osp.realpath(osp.join(dd,fi))
  raise IOError("cannot find %s in any of the following directories %s"%(fi,dirs))


class miniparse(object):
  def _parse(self,txtit):
    cont = ""
    for l in txtit:
      
      # include 
      mtch = re.search("\s*#\s*include\s*[\"|\'](.+?)[\"|\']",l)
      if mtch:
        rfn = lookupfile(mtch.group(1),self.localdir)
        self.localdir += [osp.dirname(rfn)]
        self._parse(open(rfn))

      #undef 
      mtch = re.search("\s*#\s*undef\s*(.+)",l)
      if mtch:
        name = mtch.group(1).strip()
        if name in self.pf:
          del(self.pf[name])

      
      # deal with comments
      l = re.split("#|!",l)[0]
      if len(l.strip())==0:
        continue

      k = cont
      if not cont:
        k,v = re.findall("(?<!#)((?:\w|\.)+)\s*=\s*(.+)",l)[0]
        v = v.strip()
      else:
        v = l.strip()
        v = self.pf[k]+ " "+ v
      
      cont = "" 
      if v[-1]=="&": #continuation
        v = v[:-1]
        cont = k
      
      self.pf[k] = v

  def __init__(self, pfn,**kk):
    self.pf = {}

    self.localdir = [osp.dirname(osp.abspath(pfn or "."))]

    if pfn!=None:
      print("read parameter file %s"%pfn)
      self._parse(open(pfn))  
      #pff =open(pfn)
      #txt = "\n".join([to.split("#")[0] for to in pff])+"\n"
      #pf = dict(re.findall("(?<!#)((?:\w|\.)+)\s*=\s*(.+?)\n",txt))
      #pf = parse_text(pff)
      #self.pf.update(pf)
    
    self.pf.update(kk)
    self._access_list = []

  def keys(self,prefix=""):
    return [k for k in list(self.pf.keys()) if k[:len(prefix)]==prefix]
    
  def __repr__(self):
    rr = []
    print(self._access_list)
    for v in self._access_list:
      print(v) #ICICICICI
      rr += ["%s = %s"%(v,getattr(self,v))]
    return "\n".join(rr)
    
  def __contains__(self,val):
    res = val in self.pf or val+".file" in self.pf
    if res:
      self._access_list += [val]
    return res

  @property
  def bool(self):
    return transformme(lambda val:str(val).lower() in ("t","1","true"),self)

  @property
  def bool_array(self):
    return transformme(lambda val:str(val).lower() in ("t","1","true"),self,True)
    
  @property
  def int(self):
    return transformme(int,self)

  @property
  def int_array(self):
    return transformme(int,self,True)

  @property
  def float(self):
    return transformme(float,self)

  @property
  def float_array(self):
    return transformme(float,self,True)

  @property
  def str(self):
    return transformme(str,self)

  @property
  def str_array(self):
    return transformme(str,self,True)

  def __getattr__(self,val):
    res = getattr(self.str,val)
    return res
    
def fromargv():
  import sys
  argv = sys.argv

  if len(argv)!=2:
    print("usage: %s parfile\n"%(argv[0]))
    sys.exit(-1)

  pfn = argv[1]
  pf = miniparse(pfn)
  return pf