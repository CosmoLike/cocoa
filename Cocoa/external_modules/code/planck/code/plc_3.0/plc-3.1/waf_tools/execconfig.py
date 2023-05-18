from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import waflib.TaskGen
import waflib.Task as Task
from waflib import Utils
import waflib
import os.path as osp
import os

def uniqify(lst):
  rlst = []
  for v in lst:
    #print v, rlst
    if v in rlst:
      #print "caught"
      continue
    rlst.append(v)
  return rlst
def ptrquote(st):
  res = ""
  for v in st:
    if v=='"':
      res +='\\"'
    else:
      res+=v
  return res

@waflib.TaskGen.feature("build_pkgconfig")
def build_pkgconfig(self):
  from waflib.Tools.ccroot import USELIB_VARS
  if self.flavor=='c':  
    USELIB_VARS['build_pkgconfig']   = set(['INCLUDES', 'DEFINES', 'CPPFLAGS', 'CFLAGS']+['LIB', 'STLIB', 'LIBPATH', 'STLIBPATH', 'LINKFLAGS', 'RPATH', 'LINKDEPS'])
    cf = ['CPPFLAGS', 'CFLAGS']
    addlib = ["clik"]
  else:
    USELIB_VARS['build_pkgconfig']   =set(['FCFLAGS','DEFINES','INCLUDES']+['LIB','STLIB','LIBPATH','STLIBPATH','LINKFLAGS','RPATH','LINKDEPS'])
    cf = ['FCFLAGS']
    addlib = ["clik","clik_f90"]

  #USELIB_VARS['cprogram']
  self.process_use()
  self.propagate_uselib_vars()
  vrs = dict([(v,list((self.env[v]))) for v in USELIB_VARS['build_pkgconfig']])
  includepath = ptrquote(" ".join([self.env.CPPPATH_ST%v for v in uniqify(vrs["INCLUDES"])]))
  libpath = ptrquote(" ".join([self.env.LIBPATH_ST%v for v in uniqify(vrs["LIBPATH"])]))
  rlibpath = ptrquote(" ".join([self.env.RPATH_ST%v for v in uniqify(vrs["RPATH"])]))
  stlibpath = ptrquote(" ".join([self.env.LIBPATH_ST%v for v in uniqify(vrs["STLIBPATH"])]))
  libs = ptrquote(" ".join([self.env.LIB_ST%v for v in uniqify(vrs["LIB"]+addlib)]))
  stlibs = ptrquote(" ".join([self.env.STLIB_ST%v for v in uniqify(vrs["STLIB"])]))
  defines = ptrquote(" ".join([self.env.DEFINES_ST%v for v in uniqify(vrs["DEFINES"])]))
  cfs = []
  #print cf
  for tt in cf+["LINKFLAGS"]:
    #print tt,vrs[tt]
    cfs += vrs[tt]
  #print cfs
  cflags = ptrquote(" ".join(uniqify(cfs)))

  #print "YEAH:"
  #print includepath
  #print libpath
  #print rlibpath
  #print stlibpath
  #print libs
  #print stlibs
  #print cflags
  #print defines
  
  alibs = ""
  if libs:
    alibs += (self.env.SHLIB_MARKER or "") +" ".join([rlibpath,libpath,libs])
  if stlibs:
    alibs += (self.env.STLIB_MARKER or "") +" ".join([srlibpath,stlibs])

  f=open(osp.join(self.env.BINDIR,self.target),"w")
  print(config_tpl%(" ".join((includepath,defines,cflags)),alibs), file=f)
  f.close()  
  os.chmod(osp.join(self.env.BINDIR,self.target),Utils.O755)

config_tpl = """#! /usr/bin/env python
# don't do much for now
from optparse import OptionParser
parser = OptionParser()

parser.add_option("--cflags", action="store_true",
                  help="only the cflags")
parser.add_option("--libs", action="store_true",
                  help="only libflags")

(options, args) = parser.parse_args()


res={}
cflags = "%s"
libs = "%s"

if (not options.libs) and (not options.cflags):
  options.libs=True
  options.cflags=True

if options.cflags:
  print cflags,
if options.libs:
  print libs,
print

"""

@waflib.TaskGen.feature("*")
@waflib.TaskGen.before_method('process_source')
def process_execrule(self):
  
  if not getattr(self,'execrule',None):
    return
  self.meths.remove('process_source')
  name=str(getattr(self,'name',None)or self.target or self.execrule)
  cls=Task.task_factory(name,self.execrule,getattr(self,'vars',[]),shell=getattr(self,'shell',True),color=getattr(self,'color','BLUE'))
  tsk=self.create_task(name)
  if getattr(self,'target',None):
    if isinstance(self.target,str):
      self.target=self.target.split()
    if not isinstance(self.target,list):
      self.target=[self.target]
    for x in self.target:
      if isinstance(x,str):
        tsk.outputs.append(self.path.find_or_declare(x))
      else:
        x.parent.mkdir()
        tsk.outputs.append(x)
    if getattr(self,'install_path',None):
      self.bld.install_files(self.install_path,tsk.outputs,chmod=Utils.O755)
  if getattr(self,'source',None):
    tsk.inputs=self.to_nodes(self.source)
    self.source=[]
  if getattr(self,'scan',None):
    cls.scan=self.scan
  if getattr(self,'cwd',None):
    tsk.cwd=self.cwd
  if getattr(self,'update_outputs',None)or getattr(self,'on_results',None):
    Task.update_outputs(cls)
  if getattr(self,'always',None):
    Task.always_run(cls)
  for x in['after','before','ext_in','ext_out']:
    setattr(cls,x,getattr(self,x,[]))
