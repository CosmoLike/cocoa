from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from waflib import Logs, Utils

import sys
import os.path as osp

SNIP_FUNCTION = '''
int main(int argc, char **argv) {
void (*p)();
(void)argc; (void)argv;
p=(void(*)())(%s);
return !p;
}'''
def to_header(**dct):
  if'header_name'in dct:
    dct=Utils.to_list(dct['header_name'])
    return''.join(['#include <%s>\n'%x for x in dct])
  return''

def add_lib_option(libname,opt,default="/usr/local/",install=True):
  import optparse
  grp=optparse.OptionGroup(opt.parser,libname+" options")
  grp.add_option("--%s_islocal"%libname,action="store_true",default=False,help="%s has been installed previsouly using --%s_install or --%s_forceinstall"%(libname,libname,libname))
  grp.add_option("--%s_prefix"%libname,action="store",default="",help="%s include/lib path prefix"%libname)
  grp.add_option("--%s_include"%libname,action="store",default="",help="%s include path"%libname)
  grp.add_option("--%s_lib"%libname,action="store",default="",help="%s lib path"%libname)
  grp.add_option("--%s_link"%libname,action="store",default="",help="%s link line"%libname)
  if install:
    grp.add_option("--%s_installifneeded"%libname,action="store_true",default=False,help="if %s is not found, install it try to install"%libname,dest="%s_install"%libname)
    grp.add_option("--%s_install"%libname,action="store_true",default=False,help="install %s"%libname,dest="%s_forceinstall"%libname)
  opt.add_option_group(grp)
    

def noemptylist(li):
  return [ll for ll in li if ll]
  
def opt_to_libpaths(ctx,name):
  prefix = [getattr(ctx.options,name+"_prefix","usr/local/")]
  include = [v for v in getattr(ctx.options,name+"_include","").split(":") if v]
  lib = [v for v in getattr(ctx.options,name+"_lib","").split(":") if v]
  link = libsfromlinkline(getattr(ctx.options,name+"_link"))
  if getattr(ctx.options,name+"_islocal",None):
    prefix = [ctx.env.localpref] + prefix
    
  for prf in prefix[::-1]:
    if prf:
      include = [osp.join(prf,"include")] + include 
      lib = [osp.join(prf,"lib")] + lib 
  
  return prefix,include,lib,link
      
def libsfromlinkline(ll):
  return [lb.strip() for lb in ll.split("-l") if lb.strip()]
  
def magic_join(pth,adx):
  if isinstance(pth,str):
    pth = set(pth.split(":"))
  return ":".join([osp.join(pp,adx) for pp in set(pth)])
      
def add_lib(conf,prefix,include,libpath,libname, funcname="",headername="",libs = [], uselib=[],defines=[],frameworkpath=[],framework=[],flagline=""):
  #print conf,prefix,include,libpath,libname, funcname,headername,libs , uselib,defines
  
  for inc in include:
    if inc:
      conf.env.append_value("INCLUDES_%s"%(libname),inc)
  if libs == []:
    libs = [libname]
  if type(uselib)==type(""):
    uselib = [uselib]
  if type(funcname)==type(""):
    funcname=[funcname]
  if type(defines)==type(""):
    defines=[defines]
    
  conf.parse_flags(flagline,libname)

  conf.check_cc(lib=libs, libpath = noemptylist(libpath),rpath=noemptylist(libpath) ,uselib_store=libname,mandatory=1,uselib=uselib+[libname],defines=defines,frameworkpath=frameworkpath,framework=framework)
  for fnc in funcname:
    conf.check_cc(
      errmsg="failed",msg="checking for function '%s'"%fnc,
      fragment=to_header(header_name=headername)+SNIP_FUNCTION%fnc,header_name=headername,uselib=" ".join([libname]+uselib),mandatory=1,frameworkpath=frameworkpath,framework=framework)
    conf.undefine("HAVE_"+fnc.upper())

def add_lib_f90(conf,prefix,include,libpath,libname, funcname="",headername="",libs = [], uselib=[],defines=[],frameworkpath=[],framework=[],flagline=""):
  #print conf,prefix,include,libpath,libname, funcname,headername,libs , uselib,defines
  
  for inc in include:
    if inc:
      conf.env.append_value("INCLUDES_%s"%(libname),inc)
  if libs == []:
    libs = [libname]
  if type(uselib)==type(""):
    uselib = [uselib]
  if type(funcname)==type(""):
    funcname=[funcname]
  if type(defines)==type(""):
    defines=[defines]
    
  conf.parse_flags(flagline,libname)
  # do nothing for now...
  for inc in libpath:
    if inc:
      conf.env.append_value("LIBPATH_%s"%(libname),inc)
      conf.env.append_value("RPATH_%s"%(libname),inc)
  for inc in libs:
    if inc:
      conf.env.append_value("LIB_%s"%(libname),inc)
  
  #print dir(conf)
  #conf.check_fortran(lib=libs, libpath = noemptylist(libpath),rpath=noemptylist(libpath) ,uselib_store=libname,mandatory=1,uselib=uselib+[libname],defines=defines,frameworkpath=frameworkpath,framework=framework)
  for fnc in funcname:
    #self.check_cc(fragment=FC_FRAGMENT,compile_filename='test.f',features='fc fcprogram',msg='Compiling a simple fortran app')
    #print " ".join([libname]+uselib)
    conf.check_cc(
      errmsg="failed (check whether lib is compiled in 32 or 64bits)",msg='checking for module %s'%fnc,
      uselib=" ".join([libname]+uselib),mandatory=1,frameworkpath=frameworkpath,framework=framework, fragment = "program test\n  use %s\n end program test\n"%fnc,compile_filename='test.f90',features='fc fcprogram')
    conf.undefine("HAVE_"+fnc.upper())




add_lib_dict = {
  "c" : add_lib,
  "f90" : add_lib_f90
}
def shouldIinstall_all(ctx,name):
  if ctx.options.install_all_deps==False:
    return False
  rr = [(vv,getattr(ctx.options,vv,"")) for vv in dir(ctx.options) if name in vv and vv!=name+"_forceinstall" ]
  for vv in rr:
    if vv[1]:
      return False
  return True
def upgrade(ctx,name):
  rr = [(vv,getattr(ctx.options,vv,"")) for vv in dir(ctx.options) if name in vv and vv not in (name+"_forceinstall",name+"_install") ]
  for vv in rr:
    if vv[1]:
      return False
  return getattr(ctx.options,name+"_install",False) or ctx.options.upgrade_all_deps


def conf_lib(ctx,name,_libs,testfunc=[],testinclude=[],add_inc_path=[],defines=[],frameworkpath=[],framework=[],install=False,msg="",uselib=[],flagline="",opt_name="",add_lib_code="c",forceinstall=False):
  if not opt_name:
    opt_name=name
  # do install if needed
  #print install and (getattr(ctx.options,opt_name+"_install",False) or getattr(ctx.options,opt_name+"_forceinstall",False))
  #print install ,(getattr(ctx.options,opt_name+"_install",False) ,getattr(ctx.options,opt_name+"_forceinstall",False))
  iall =shouldIinstall_all(ctx,name)
  if install and (upgrade(ctx,name) or forceinstall or getattr(ctx.options,opt_name+"_forceinstall",False) or iall):
    # first try without install !
    setattr(ctx.env,"has_"+name,True)
    setattr(ctx.options,"%s_islocal"%opt_name,True)
    try:
      #print forceinstall==False and getattr(ctx.options,opt_name+"_forceinstall",False)==False
      #print forceinstall==False,getattr(ctx.options,opt_name+"_forceinstall",False)==False
      #print "LALALA",getattr(ctx.options,opt_name+"_forceinstall","MUMUMU")
      assert forceinstall==False and getattr(ctx.options,opt_name+"_forceinstall",False)==False and iall==False
      conf_lib(ctx,name,_libs,testfunc,testinclude,add_inc_path,defines,frameworkpath,framework,False,msg,uselib,flagline,opt_name,add_lib_code)
      return
    except Exception as e:
      if forceinstall==False and getattr(ctx.options,opt_name+"_forceinstall",False)==False and iall==False:
        Logs.pprint("RED","%s not found, try to install it"%name)
      install(ctx)
  # compute paths
  prefix,include,lib,link = opt_to_libpaths(ctx,opt_name)
  
  # libs to test for
  libs=_libs
  # if a link option is given, it overides the libs
  if link:
    libs = link
    
  # extra includes ?
  if add_inc_path:
    extinc = []
    for adi in add_inc_path:
      extinc += [osp.join(inc,adi) for inc in include]
    include += extinc
  
  try:
    add_lib_dict[add_lib_code](ctx,prefix,include,lib,name,
          testfunc,testinclude,
          libs=libs,uselib=uselib,defines=defines,frameworkpath=frameworkpath,framework=framework,flagline=flagline)

    setattr(ctx.env,"use_%s"%name,name)
    setattr(ctx.env,"has_%s"%name,name)

  except Exception as e:
    ctx.env["INCLUDES_%s"%name]=[]
    ctx.env["DEFINES_%s"%name]=[]
    if not getattr(ctx.env,"has_"+name,False):
      if not getattr(ctx.env,"silent_"+name,False):
        Logs.pprint("BLUE","Optional %s not found"%name)
        Logs.pprint("BLUE","Compilation will continue without it")
    else:
      Logs.pprint("RED","%s not found"%name)
      Logs.pprint("PINK", "check that %s_prefix or %s_lib and %s_include command line options point toward your %s install"%(opt_name,opt_name,opt_name,opt_name))
      Logs.pprint("PINK", "or check that %s is compiled in %d bit"%(name,{True:64}.get(ctx.options.m64,32)))
      if msg:
        Logs.pprint("PINK", msg)      
      if install:
        Logs.pprint("PINK", "or install automatically using cmdline option --%s_install"%(opt_name))      
      raise e

def getfromurl(fromurl,tofile):
  try:
    import urllib.request, urllib.error, urllib.parse
    urlopen = urllib.request.urlopen
  except ImportError:
    # are we on py2 ?
    import urllib2
    urlopen = urllib2.urlopen
  
  try:
    luaf = urlopen(fromurl)
    aa = luaf.read()
    luaf.close()
  except Exception as e:
    import subprocess
    print("curl %s %s"%(fromurl,tofile))
    aa = subprocess.check_output(["curl",fromurl])

  f=open(tofile,"wb")
  f.write(aa)  
  f.close()
  
def installsmthg_pre(ctx,where,what,whereto="build/"):

  from waflib import Utils,Errors
  import re
  import os.path as osp
  import tarfile
  import shutil
  import os
  
  #ctx.env = Environment.Environment(filename="build/c4che/default.cache.py")
  Logs.pprint("PINK","Install '%s'"%what)
  if osp.exists(osp.join(whereto,what)):
    Logs.pprint("PINK","%s already downloaded"%what)
  else:
    Logs.pprint("PINK","download from "+where)
    getfromurl(where,osp.join(whereto,what))
    
  try:
    tf = tarfile.open(osp.join(whereto,what))
  except tarfile.ReadError:
    os.remove(osp.join(whereto,what))
    return installsmthg_pre(ctx,where,what,whereto)
  #Logs.pprint("RED","LALALALA")
  ss = set([v.name.split("/")[0] for v in tf.getmembers()])
  assert len(ss)==1,"cannot understand the directory structure of %s"%osp.join(whereto,what)
  for ff in [ff.name for ff in tf.getmembers()]:
    if osp.exists(osp.join(whereto,ff)):
      if osp.isdir(osp.join(whereto,ff)):
        shutil.rmtree(osp.join(whereto,ff))
      else:
        os.remove(osp.join(whereto,ff))
  tf.close()
  Logs.pprint("PINK","untar "+what)
  if ctx.exec_command("cd %s/;tar -zxf %s"%(whereto,what))!=0:
    raise Errors.WafError("Cannot untar "+what)
  return ss.pop()

def installsmthg_post(ctx,where,what,extra_config=""):
  from waflib import Utils,Errors
  CCMACRO = "\"%s %s\""%(ctx.env.CC[0]," ".join(ctx.env.mopt))
  CCMACRO = "CC=%s CXX=%s "%(CCMACRO,CCMACRO)
  CPPMACRO = "CPP=\"%s -E\" CXXCPP=\"g++ -E\" "%(ctx.env.CC[0])
  cmdline = "cd build/%s; ./configure --prefix=%s %s  %s %s; make clean;make -j %d ;make install"%(where,ctx.env.mprefix,extra_config,CCMACRO, CPPMACRO,ctx.options.jobs)
  Logs.pprint("PINK",cmdline)
  if ctx.exec_command(cmdline)!=0:
    raise Errors.WafError("Cannot build %s"%what)
  #Logs.pprint("GREEN","You can now run ./waf configure, adding the following option '--%s_islocal'"%what)

def check_python_module(ctx,name,extracmd=""):
  import sys
  import imp
  import site
  site.addsitedir(ctx.env.PYTHONDIR)
  if ctx.env.PYTHONDIR not in sys.path:
    sys.path=[ctx.env.PYTHONDIR]+sys.path
  try:
    ctx.start_msg("Checking python module '%s'"%name)
    imp.find_module(name)
    if extracmd:
      #print extracmd
      exec(extracmd)
    ctx.end_msg(True)
  except Exception as e: 
    ctx.end_msg(False)
    raise e

def add_python_option(ctx,name):
  import optparse
  grp = ctx.parser.get_option_group("--nopyo")
  if grp==None:
    grp=optparse.OptionGroup(ctx.parser,"python module options")
  
  grp.add_option("--%s_install"%name,action="store_true",default=False,help="install %s"%name,dest="%s_forceinstall"%name)
  grp.add_option("--%s_installifneeded"%name,action="store_true",default=False,help="if %s is not found, install it try to install"%name,dest="%s_install"%name)
  ctx.add_option_group(grp)  

def configure_python_module(ctx,name,url,packtgz,pack,cmdline=None,extracmd="",forceinstall=False,postinstall=None):
  import waflib.Logs
  import os
  from waflib import Errors
  import os.path as osp
  import autoinstall_lib as atl
  ctx.load("python")
  doit = False
  import sys

  iall = shouldIinstall_all(ctx,name)
  try:
    assert forceinstall==False and getattr(ctx.options,name+"_forceinstall")==False and iall==False
    check_python_module(ctx,name,extracmd)
  except Exception as e: 
    if upgrade(ctx,name) or getattr(ctx.options,name+"_forceinstall",False) or iall:
      waflib.Logs.pprint("PINK","Install python module '%s'"%name)
      cmdline = "PYTHONUSERBASE=%s %s -m pip install --user --upgrade %s"%(ctx.env["PREFIX"],ctx.env["PYTHON"][0],name)
      ret = ctx.exec_command(cmdline)
      waflib.Logs.pprint("PINK",cmdline)
      ret = ctx.exec_command(cmdline)
      if ret!=0:
        # default to old way
        configure_python_module_old(ctx,name,url,packtgz,pack,cmdline,extracmd,forceinstall,postinstall)
      check_python_module(ctx,name,extracmd)
      if postinstall:
        postinstall()

def configure_python_module_old(ctx,name,url,packtgz,pack,cmdline=None,extracmd="",forceinstall=False,postinstall=None):
  import waflib.Logs
  import os
  from waflib import Errors
  import os.path as osp
  import autoinstall_lib as atl
  ctx.load("python")
  doit = False
  import sys

  iall = shouldIinstall_all(ctx,name)
  try:
    assert forceinstall==False and getattr(ctx.options,name+"_forceinstall")==False and iall==False
    check_python_module(ctx,name,extracmd)
  except Exception as e: 
    if upgrade(ctx,name) or getattr(ctx.options,name+"_forceinstall",False) or iall:
      waflib.Logs.pprint("PINK","Install python module '%s'"%name)
      packname = atl.installsmthg_pre(ctx,url,packtgz)
      if not osp.exists(ctx.env.PYTHONDIR):
        os.makedirs(ctx.env.PYTHONDIR)
      if cmdline==None:
        cmdline =  "cd build/%s; PYTHONPATH=%s:$PYTHONPATH %s setup.py build_ext -L=%s ;PYTHONPATH=%s:$PYTHONPATH %s setup.py install --install-lib=%s --install-scripts=%s"%(packname,ctx.env.PYTHONDIR,ctx.env.PYTHON[0],ctx.env.LIBPATH_PYEMBED[0],ctx.env.PYTHONDIR,ctx.env.PYTHON[0],ctx.env.PYTHONDIR,ctx.env.BINDIR)
      waflib.Logs.pprint("PINK",cmdline)
      ret = ctx.exec_command(cmdline)
      if ret!=0:
        raise Errors.ConfigurationError("Cannot build %s"%name)
      # deal with eggs...
      if (not osp.exists(osp.join(ctx.env.PYTHONDIR,name))) and (not osp.exists(osp.join(ctx.env.PYTHONDIR,name+".py"))):
        eggdir = [v for v in os.listdir(ctx.env.PYTHONDIR) if name in v and osp.isdir(osp.join(ctx.env.PYTHONDIR,v))]
        if eggdir:
          eggdir = eggdir[0]
        else:
          eggdir = name
        if eggdir!=name and eggdir!=name+".py":
          mdir = [v for v in os.listdir(osp.join(ctx.env.PYTHONDIR,eggdir)) if name in v][0]
          import os
          os.symlink(osp.join(ctx.env.PYTHONDIR,eggdir,mdir),osp.join(ctx.env.PYTHONDIR,name))
      check_python_module(ctx,name,extracmd)
      if postinstall:
        postinstall()
    else:
      raise e

def get_lib_url(ctx,packname,default=["",""]):
  if packname+"_url" in ctx.env:
    #print("-->>",ctx.env[packname+"_url"],ctx.env[packname+"_tar"])
    return ctx.env[packname+"_url"],ctx.env[packname+"_tar"]
  return default

