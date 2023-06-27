import autoinstall_lib as atl
from waflib import Logs

version = "lua-5.1.4"
tool = "lua"

print "-> loading %s autoinstall (using version %s)"%(tool,version)

def options(opt):
  atl.add_lib_option(tool,opt,install=True)
  
def configure(ctx):
  atl.conf_lib(ctx,tool,"lua","lua_newstate","lua.h",defines=["HAS_LUA"],install=installlua)
    
def installlua(ctx):
  atl.installsmthg_pre(ctx,"http://www.lua.org/ftp/"+version+".tar.gz",version+".tar.gz")
  import Options, Environment
  import re
  import sys
  plat = "linux"
  if sys.platform.lower()=="darwin":
    plat = "macosx"

  basepath = "build/%s/src/"%version
  Logs.pprint("PINK","patching %sMakefile"%basepath)
  f=open("%sMakefile"%basepath)
  mk=f.read()
  f=open("%s/Makefile"%basepath,"w")
  mk = re.sub("PLAT= none","PLAT= "+plat,mk)
  mk = re.sub("CFLAGS= -O2 -Wall ","CFLAGS= -O2 -Wall %s -fPIC "%(" ".join(ctx.env.mopt)),mk)
  mk = re.sub("CC= gcc","CC= gcc %s "%(" ".join(ctx.env.mopt)),mk)
  print >>f,mk
  f.close()
  
  basepath = "build/%s/"%version
  Logs.pprint("PINK","patching %sMakefile"%basepath)
  f=open("%sMakefile"%basepath)
  mk=f.read()
  f=open("%sMakefile"%basepath,"w")
  mk = re.sub("PLAT= none","PLAT= "+plat,mk)
  mk = re.sub("/usr/local","%s"%ctx.env.mprefix,mk)
  print >>f,mk
  f.close()
  
  basepath = "build/%s/src/"%version
  Logs.pprint("PINK","patching %sluaconf.h"%basepath)
  f=open("%s/luaconf.h"%basepath)
  mk=f.read()
  f=open("%s/luaconf.h"%basepath,"w")
  mk = re.sub("/usr/local","%s"%ctx.env.mprefix,mk)
  print >>f,mk
  f.close()

  atl.installsmthg_post(ctx,version,tool)

