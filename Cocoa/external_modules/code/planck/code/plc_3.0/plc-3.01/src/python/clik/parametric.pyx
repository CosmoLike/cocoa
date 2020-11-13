cimport numpy as nm
import numpy as nm
nm.import_array()
cimport libc.stdlib as stdlib
cimport libc.stdio as stdio
import os.path as osp
import re

def loadcolumn(pth):
  reg = re.search("(.+)\[(\d+)\]",pth)
  if reg:
    pth = reg.group(1)
    col = int(reg.group(2))
    return nm.loadtxt(pth)[:,col]
  else:
    #print pth
    return nm.loadtxt(pth)

cdef extern from "errorlist.h":
  ctypedef struct error:
    pass
    
  void stringError(char* str, error *err)
  int getErrorValue(error* err)
  int isError(error *err) 
  void purgeError(error **err) 
  void printError(void* flog,error* err)

class CError(Exception):
  def __init__(self,val,str):
    self.val=val
    self.comment=str
  def __str__(self):
    return self.comment.strip()

cdef doError(error **err):
  cdef char estr[10000]
  if (isError(err[0])):
    stringError(estr,err[0])
    er=CError(getErrorValue(err[0]),estr)
    purgeError(err)
    
    return er
  return None

cdef extern from "clik_parametric.h":
  ctypedef struct c_parametric "parametric":
    int lmin,lmax,ndet,nfreq,nvar

  c_parametric *parametric_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
  void parametric_free(void** pegl)
  void parametric_compute(c_parametric *egl, double *pars, double* Rq, double *dRq, error **err)
  double parametric_get_default(c_parametric* egl,char *key, error **err)
  double parametric_get_value(c_parametric *egl, char *key, error **err)
  void parametric_dnofail(c_parametric* egl, int vl)
  void parametric_set_color(c_parametric *egl,double *color, error **err)
  void parametric_set_void(c_parametric *egl, int nvoid, int *voidlist,error **err)


  c_parametric *powerlaw_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
  c_parametric *powerlaw_free_emissivity_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)

cdef set_color(c_parametric *egl, object color,int ndet):
  cdef double *_color
  cdef error *_err
  cdef error **err

  _err = NULL
  err = &_err
  if color is not None:
    _color = <double*>stdlib.malloc(sizeof(double)*ndet*ndet)
    assert _color!=NULL,"Can't allocate color"
    if len(color)==ndet:
      for i in range(ndet):
        for j in range(ndet):
          _color[i*ndet+j] = color[i]*color[j]
    else:
      assert len(color)==ndet,"bad lenght for color vector"
      for i in range(ndet*ndet):
        _color[i] = color[i]
    parametric_set_color(egl,_color,err)
    er=doError(err)
    if er:
      raise er


cdef class parametric:
  
  def __cinit__(self):
    self.initfunc=NULL
    self.ismul=0
    self.component_name=""

  def set_name(self,name):
    self.component_name = name
  def get_name(self):
    if self.component_name =="":
      return self.__class__.__name__
    else:
      return "%s"%self.component_name
      
  def test_derivative(self,key,value=None,h = 1e-3):
    """test_derivative(key,value=None,h = 1e-3) -> exact, approx, diff, reldiff"""
    if value == None:
      value = self.parvalues
    
    value = nm.array(value)*1.
    ref,dref = self(value,True) 
  
    ik = self.varpar.index(key)
    vp = value*1.
    if vp[ik]==0:
      vp[ik] = h
      step = h
    else:
      vp[ik] = vp[ik]*(1+h)
      step = value[ik]*h
    vm = value*1.
    if vm[ik]==0:
      vm[ik] = -h
    else:
      vm[ik] = vm[ik]*(1-h)
    rp = self(vp)
    rm = self(vm)
    dap = (rp-rm)/(2*step)
    return dref[ik],dap,dref[ik]-dap,(dref[ik]-dap)/dref[ik]

  def __init__(self,detlist,vars,lmin,lmax,defs={},dnofail=False,color=None,voidmask=None,rename={}):
    """__init__(self,detlist,vars,lmin,lmax,defs={},dnofail=False,color=None,voidmask=None,rename={})"""
    cdef double p_detlist[2000]
    cdef char *defkey[2000]
    cdef char *defvalue[2000]
    cdef char *key[2000]
    cdef error *_err
    cdef error **err

    
    _err = NULL
    err = &_err
    
    ndef = len(defs)
    ndet = len(detlist)

    for i in range(ndet):
      p_detlist[i] = detlist[i]

    i = 0
    tmp_nk = []
    tmp_v =[]
    for k,v in defs.items():
      nk = rename.get(k,k)
      tmp_nk += [str.encode(nk)]
      defkey[i] = tmp_nk[-1]
      tmp_v += [str.encode(v)]
      defvalue[i] = tmp_v[-1]
      i+=1
    
    nvar = len(vars)
    tmp_key = []
    for i in range(nvar):
      tmp_key += [str.encode(rename.get(vars[i],vars[i]))]
      key[i] = tmp_key[-1]


    self.rename = rename
    self.emaner = dict([(rename[k],k) for k in rename])


    if self.initfunc==NULL:
      raise NotImplementedError("Must fill self.initfunc with a valid c function")
      #self.celf = parametric_init(ndet,p_detlist,ndef,defkey,defvalue,nvar,key,lmin,lmax,err)
    else:
      self.celf = (<simple_init>self.initfunc)(ndet,p_detlist,ndef,defkey,defvalue,nvar,key,lmin,lmax,err)
    er=doError(err)
    if er:
      raise er
    
    self._post_init(detlist,vars,lmin,lmax,defs,dnofail,color,voidmask,rename)

  def _post_init_1(self,color,ndet,voidmask,defs):
    cdef int voidlist[2000]
    cdef error *_err
    cdef error **err
    _err = NULL
    err = &_err
    
    self.defaults = defs

    set_color(self.celf,color,ndet)

    self.voidmask=None
    if voidmask:
      _voidlist = [i for i in range(len(voidmask)) if not bool(int(voidmask[i]))]
      nvoid = len(_voidlist)
      if nvoid!=0:
        for i in range(nvoid):
          voidlist[i]=_voidlist[i]
        parametric_set_void(self.celf,nvoid,voidlist,err)
        er=doError(err)
        if er:
          raise er
      self.voidmask = voidmask

  def _post_init_2(self,prs,lmin,lmax):
    dv = []
    for p in prs:
      if self.has_parameter(p):
        dv += [self.get_default_value(p)]
      else:
        dv +=[0]

    self.varpar = prs
    self.parvalues = dv
    #self.varpar = OrderedDict(zip(prs,dv))
    self.nell = lmax+1-lmin    
    self.name = self.__class__.__name__
    #print self.name
    #print dir(self)

  def _post_init(self,detlist,vars,lmin,lmax,defs,dnofail,color,voidmask,*other):
    cdef double _color[2000]
    cdef int i
    cdef error *_err
    cdef error **err
    cdef int voidlist[2000]

    _err = NULL
    err = &_err
    
    self._post_init_1(color,len(detlist),voidmask,defs)

    parametric_dnofail(self.celf,int(dnofail))
    prs = vars
    #print prs
    if not dnofail:
      prs = [p for p in vars if self.has_parameter(p)]
      #print prs
      if len(prs)!= len(vars):
        parametric_free(<void**>&(self.celf))
        self.__init__(detlist,prs,lmin,lmax,defs,False,color,voidmask,*other)

    self._post_init_2(prs,lmin,lmax)
    
  def get_default_value(self,key):
    cdef error *_err
    cdef error **err
    _err = NULL
    err = &_err
    
    nk = str.encode(self.rename.get(key,key))
    res = parametric_get_default(self.celf,nk, err)
    er=doError(err)
    if er:
      raise er
    return res


  def has_parameter(self,key):
    try:
      #print key
      self.get_default_value(key)
      return True
    except Exception,e:
      return False

  def __call__(self,pars,derivatives=False):
    cdef error *_err
    cdef error **err
    cdef double *_drq
    cdef double *_rq
      
    
    if len(pars)!=self.celf.nvar:
      raise Exception("Bad shape (expecting (%d) got (%d))"%(self.celf.nvar,len(pars)))
    rq = nm.zeros((self.nell,self.celf.ndet,self.celf.ndet),dtype=nm.double)
    _rq = <double*> nm.PyArray_DATA(rq)
    if derivatives:
      drq = nm.zeros((self.celf.nvar,self.nell,self.celf.ndet,self.celf.ndet),dtype=nm.double)
      _drq = <double*> nm.PyArray_DATA(drq)
      
    else:
      _drq = NULL
    pars_proxy=nm.PyArray_ContiguousFromAny(pars,nm.NPY_DOUBLE,1,1)
    _err = NULL
    err = &_err
    parametric_compute(self.celf,  <double*> nm.PyArray_DATA(pars_proxy), _rq,_drq, err)
    er=doError(err)
    if er:
      raise er
    if derivatives:
      return rq,drq
    return rq
    
  def __dealloc__(self):
    if self.celf!=NULL:
      parametric_free(<void**>&(self.celf))

  def get_template(self,data_dir="",data_path="",data_file="",data=None):
    return None


cdef class powerlaw(parametric):
  def __cinit__(self):
    self.initfunc = <void*>powerlaw_init

cdef class powerlaw_free_emissivity(parametric):
  def __cinit__(self):
    self.initfunc = <void*>powerlaw_free_emissivity_init

def get_data_path(plg=None):
  import os
  if "CLIK_DATA" in os.environ:
    res=os.environ["CLIK_DATA"]
    if plg:
      return osp.join(res,plg)
  return ""

cdef class parametric_mul(parametric):
  def _post_init_2(self,prs,lmin,lmax):
    parametric._post_init_2(self,prs,lmin,lmax)
    self.ismul=1

cdef class parametric_template(parametric):
  def __cinit__(self):
    self.initfunc = NULL
    self.template_name = ""
    self.plugin_name = ""

  def __init__(self,detlist,vars,lmin,lmax,defs={},dnofail=False,color=None,voidmask=None,rename={},data_dir="",data_path="",data_file="",data=None):
    """__init__(self,detlist,vars,lmin,lmax,defs={},dnofail=False,color=None,voidmask=None,rename={},data_dir="",data_path="",data_file="",data=None)"""
    cdef double p_detlist[2000]
    cdef char *defkey[2000]
    cdef char *defvalue[2000]
    cdef char *key[2000]
    cdef double *template
    cdef error *_err
    cdef error **err

    _err = NULL
    err = &_err
    
    ndef = len(defs)
    ndet = len(detlist)
    
    for i in range(ndet):
      p_detlist[i] = detlist[i]

    i = 0
    tmp_nk = []
    tmp_v =[]
    for k,v in defs.items():
      nk = rename.get(k,k)
      tmp_nk += [str.encode(nk)]
      defkey[i] = tmp_nk[-1]
      tmp_v += [str.encode(v)]
      defvalue[i] = tmp_v[-1]
      i+=1
        
    nvar = len(vars)
    tmp_key = []
    for i in range(nvar):
      tmp_key += [str.encode(rename.get(vars[i],vars[i]))]
      key[i] = tmp_key[-1]


    self.rename = rename
    self.emaner = dict([(rename[k],k) for k in rename])
    
    self._template=None
    tmp = self.get_template(data_dir,data_path,data_file,data)
    self._template = tmp*1.
    if len(self._template.flat[:])==0:
      template=NULL
    else:
      template = <double*> nm.PyArray_DATA(self._template)
    
    if self.initfunc==NULL:
      raise NotImplementedError("Must fill self.initfunc with a valid c function")
      #self.celf = parametric_init(ndet,p_detlist,ndef,defkey,defvalue,nvar,key,lmin,lmax,err)
    else:
      self.celf = (<template_init>self.initfunc)(ndet,p_detlist,ndef,defkey,defvalue,nvar,key,lmin,lmax,template,err)
        
    er=doError(err)
    if er:
      raise er

    self._post_init(detlist,vars,lmin,lmax,defs,dnofail,color,voidmask,rename,data_dir,data_path,data_file,data)

  
  def get_template(self,data_dir="",data_path="",data_file="",data=None):
    if data_dir=="" and data_path=="" and data_file=="" and data is None and self._template is not None:
      return self._template
    if data is None:
      if data_path:
        pth = data_path
      else:
        bpth = get_data_path(self.plugin_name)
        if data_dir:
          bpth = data_dir
        fpth = self.template_name
        if data_file:
          fpth = data_file  
        if isinstance(fpth,str):
          fpth = [fpth]
        pth = [osp.join(bpth,fpthI) for fpthI in fpth]
      if isinstance(pth,str):
        pth = [pth]
      
      tmp = nm.concatenate([loadcolumn(pp) for pp in pth])
    else:
      
      tmp = nm.array(data)
    return tmp

cdef class parametric_template_mul(parametric_template):
  def _post_init_2(self,prs,lmin,lmax):
    parametric_template._post_init_2(self,prs,lmin,lmax)
    self.ismul=1
    
cdef class parametric_pol(parametric):
  def __init__(self,detlist_T,detlist_P,has_TEB,vars,lmin,lmax,defs={},dnofail=False,color=None,voidmask=None,rename={}):
    """__init__(self,detlist,vars,lmin,lmax,defs={},dnofail=False,color=None,voidmask=None,rename={})"""
    cdef double p_detlist[2000]
    cdef char *defkey[2000],
    cdef char *defvalue[2000]
    cdef char *key[2000]
    cdef error *_err
    cdef error **err
    cdef int p_has_TEB[3]
    
    _err = NULL
    err = &_err
  
    ndef = len(defs)
    ndet_T = len(detlist_T)
    ndet_P = len(detlist_P)

    for i in range(ndet_T):
      p_detlist[i] = detlist_T[i]
    for i in range(ndet_P):
      p_detlist[i+ndet_T] = detlist_P[i]

    i = 0
    tmp_nk = []
    tmp_v =[]
    for k,v in defs.items():
      nk = rename.get(k,k)
      tmp_nk += [str.encode(nk)]
      defkey[i] = tmp_nk[-1]
      tmp_v += [str.encode(v)]
      defvalue[i] = tmp_v[-1]
      i+=1
    
    
    p_has_TEB[0] = int(has_TEB[0])
    p_has_TEB[1] = int(has_TEB[1])
    p_has_TEB[2] = int(has_TEB[2])

    nvar = len(vars)
    tmp_key = []
    for i in range(nvar):
      tmp_key += [str.encode(rename.get(vars[i],vars[i]))]
      key[i] = tmp_key[-1]

    
    self.rename = rename
    self.emaner = dict([(rename[k],k) for k in rename])

    
    if self.initfunc==NULL:
      raise NotImplementedError("Must fill self.initfunc with a valid c function")
      #self.celf = parametric_init(ndet,p_detlist,ndef,defkey,defvalue,nvar,key,lmin,lmax,err)
    else:
      self.celf = (<pol_init>self.initfunc)(ndet_T,ndet_P,p_has_TEB,p_detlist,ndef,defkey,defvalue,nvar,key,lmin,lmax,err)
    er=doError(err)
    if er:
      raise er
    
    

    self._post_init(detlist_T,detlist_P,has_TEB,vars,lmin,lmax,defs,dnofail,color,voidmask,rename)

  def _post_init(self,detlist_T,detlist_P,has_TEB,vars,lmin,lmax,defs,dnofail,color,voidmask,*other):
    cdef double _color[2000]
    cdef int i
    cdef error *_err
    cdef error **err
    cdef int voidlist[2000]

    _err = NULL
    err = &_err

    self._post_init_1(color,len(detlist_T)+len(detlist_P),voidmask,defs)
    
    parametric_dnofail(self.celf,int(dnofail))
    prs = vars
    if not dnofail:
      prs = [p for p in vars if self.has_parameter(p)]
      if len(prs)!= len(vars):
        parametric_free(<void**>&(self.celf))
        self.__init__(detlist_T,detlist_P,has_TEB,prs,lmin,lmax,defs,False,color,voidmask,*other)

    self._post_init_2(prs,lmin,lmax)

cdef class parametric_pol_template(parametric_pol):
  def __cinit__(self):
    self.initfunc = NULL
    self.template_name = ""
    self.plugin_name = ""

  def __init__(self,detlist_T,detlist_P,has_TEB,vars,lmin,lmax,defs={},dnofail=False,color=None,voidmask=None,rename={},data_dir="",data_path="",data_file="",data=None):
    """__init__(self,detlist,vars,lmin,lmax,defs={},dnofail=False,color=None,voidmask=None,rename={},data_dir="",data_path="",data_file="",data=None)"""
    cdef double p_detlist[2000]
    cdef char *defkey[2000]
    cdef char *defvalue[2000]
    cdef char *key[2000]
    cdef double *template
    cdef error *_err
    cdef error **err
    cdef int p_has_TEB[3]
    
    _err = NULL
    err = &_err
    
    ndef = len(defs)
    ndet_T = len(detlist_T)
    ndet_P = len(detlist_P)

    for i in range(ndet_T):
      p_detlist[i] = detlist_T[i]
    for i in range(ndet_P):
      p_detlist[i+ndet_T] = detlist_P[i]

    i = 0
    tmp_nk = []
    tmp_v =[]
    for k,v in defs.items():
      nk = rename.get(k,k)
      tmp_nk += [str.encode(nk)]
      defkey[i] = tmp_nk[-1]
      tmp_v += [str.encode(v)]
      defvalue[i] = tmp_v[-1]
      i+=1
    
    p_has_TEB[0] = int(has_TEB[0])
    p_has_TEB[1] = int(has_TEB[1])
    p_has_TEB[2] = int(has_TEB[2])

    nvar = len(vars)
    tmp_key = []
    for i in range(nvar):
      tmp_key += [str.encode(rename.get(vars[i],vars[i]))]
      key[i] = tmp_key[-1]


    self.rename = rename
    self.emaner = dict([(rename[k],k) for k in rename])

    
    self._template=None
    tmp = self.get_template(data_dir,data_path,data_file,data)
    self._template = tmp*1.
    if len(self._template.flat[:])==0:
      template=NULL
    else:
      template = <double*> nm.PyArray_DATA(self._template)
    
    if self.initfunc==NULL:
      raise NotImplementedError("Must fill self.initfunc with a valid c function")
      #self.celf = parametric_init(ndet,p_detlist,ndef,defkey,defvalue,nvar,key,lmin,lmax,err)
    else:
      self.celf = (<pol_template_init>self.initfunc)(ndet_T,ndet_P,p_has_TEB,p_detlist,ndef,defkey,defvalue,nvar,key,lmin,lmax,template,err)
        
    er=doError(err)
    if er:
      raise er

    self._post_init(detlist_T,detlist_P,has_TEB,vars,lmin,lmax,defs,dnofail,color,voidmask,rename,data_dir,data_path,data_file,data)

  
  def get_template(self,data_dir="",data_path="",data_file="",data=None):
    if data_dir=="" and data_path=="" and data_file=="" and data is None and self._template is not None:
      return self._template
    if data is None:
      if data_path:
        pth = data_path
      else:
        bpth = get_data_path(self.plugin_name)
        if data_dir:
          bpth = data_dir
        fpth = self.template_name
        if data_file:
          fpth = data_file  
        if isinstance(fpth,str):
          fpth = [fpth]
        pth = [osp.join(bpth,fpthI) for fpthI in fpth]
      if isinstance(pth,str):
        pth = [pth]
      
      tmp = nm.concatenate([loadcolumn(pp) for pp in pth])
    else:
      
      tmp = nm.array(data)
    return tmp

component_list = ["powerlaw","powerlaw_free_emissivity"]
simple_parametric_list = component_list

def register_plugin(plg,gl,verb):
  import sys
  if verb:
    print(plg)
  #print("import %s"%plg)
  mlg =__import__("clik."+plg,fromlist=[plg])
  #print("done %s"%plg)
  global component_list
  component_list += mlg.component_list
  for cp in mlg.component_list:
    #print("register %s"%cp)
    setattr(gl,cp,getattr(mlg,cp))
    if verb:
      print("add %s"%cp)

import sys
def register_all(gl=sys.modules[__name__],verb=False):
  #print gl
  import os  
  plgs = [plg.strip() for plg in os.environ.get("CLIK_PLUGIN","").split(",") if plg.strip()]

  for plg in plgs:
    try:
      register_plugin(plg,gl,verb)
    except Exception as e:
      print ("cannot register %s (%s)"%(plg,e))
      #print e
      pass


def rename_machine(component, bdefs, rename_func=None,data_dir="",data_path="",data_file="",data=None):
  import types
  renf = rename_func
  if rename_func is None:
    renf = norename
  if isinstance(rename_func,dict):
    renf = rename_dict(rename_func)
  def rename_update(defs,vars,rename):
    rename = rename.copy()
    rups = {}
    bdef = bdefs.copy()
    bdef.update(defs)
    vv = tuple(vars)+tuple(bdef.keys())+tuple(rename.values())
    for v in vv:
      renf(v,rups)
    for k in rename:
      if rename[k] in rups:
        oo = rename[k]
        rename[k] = oo
        del(rups[k])
    rename.update(rups)
    return rename,bdef

  rmch = type(component.__name__,(component,),{})

  rdata_dir = data_dir
  rdata_path = data_path
  rdata_file = data_file
  rdata = data

  if issubclass(component,parametric_pol):  
    if issubclass(component,parametric_pol_template):
      def __init__(self,detlist_T,detlist_P,has_TEB,vars,lmin,lmax,defs={},dnofail=False,color=None,voidmask=None,rename={},data_dir="",data_path="",data_file="",data=None):
        rename,bdef = rename_update(defs,vars,rename)
        component.__init__(self,detlist_T,detlist_P,has_TEB,vars,lmin,lmax,bdef,dnofail,color,voidmask,rename, rdata_dir, rdata_path, rdata_file, rdata)
    else:      
      def __init__(self,detlist_T,detlist_P,has_TEB,vars,lmin,lmax,defs={},dnofail=False,color=None,voidmask=None,rename={}):
        rename,bdef = rename_update(defs,vars,rename)
        component.__init__(self,detlist_T,detlist_P,has_TEB,vars,lmin,lmax,bdef,dnofail,color,voidmask,rename)    
  else:
    if issubclass(component,parametric_template):
      def __init__(self,detlist,vars,lmin,lmax,defs={},dnofail=False,color=None,voidmask=None,rename={},data_dir="",data_path="",data_file="",data=None):
        rename,bdef = rename_update(defs,vars,rename)
        component.__init__(self,detlist,vars,lmin,lmax,bdef,dnofail,color,voidmask,rename, rdata_dir, rdata_path, rdata_file, rdata)
    else:
      def __init__(self,detlist,vars,lmin,lmax,defs={},dnofail=False,color=None,voidmask=None,rename={}):
        rename,bdef = rename_update(defs,vars,rename)
        component.__init__(self,detlist,vars,lmin,lmax,bdef,dnofail,color,voidmask,rename)
  if sys.version_info[0]==2:
    rmch.__init__ = types.UnboundMethodType(__init__,None,rmch)
  else:
    rmch.__init__ = __init__
  
  return rmch


def norename(v,rups):
  return

def rename_replace(before,after):
  def rename(v,rups):
    if after in v:
      rv = v.replace(after,before)
      rups[v]=rv
  return rename

def rename_dict(renint):
  def rename(v,rups):
    if v in renint: 
      rups[v] = renint[v]
    
  return rename

register_all()