import autoinstall_lib as atl
from waflib import Logs
from waflib import Utils,Errors
import waflib

import os.path as osp
    
def options(ctx):
  atl.add_lib_option("healpix",ctx,install=True)
  
twice = False
def configure(ctx):
  #atl.conf_lib(ctx,"cfitsio",["cfitsio"],"fits_init_cfitsio","fitsio.h",msg="",opt_name="healpix",uselib=["cshlib"],install=install_cfitsio)
  ctx.env.th = False
  atl.conf_lib(ctx,"chealpix",["chealpix","m","cfitsio"],"pix2vec_ring","chealpix.h",msg="or check that the path also point toward your cfitsio install",opt_name="healpix",uselib=["cpic","cfistio"],install=install_healpix)
  atl.conf_lib(ctx,"healpix_f90",["healpix","cfitsio"],"HEALPIX_TYPES",msg="or check that the path also point toward your cfitsio install",opt_name="healpix",add_lib_code="f90",uselib=["fpic","cfistio"],install=install_healpix)

#def configure_old(ctx):
#  iall = atl.shouldIinstall_all(ctx,"healpix")
#  if ctx.options.healpix_install or ctx.options.healpix_forceinstall or iall:
#    #print "do install"
#    ctx.options.healpix_islocal=True
#    ctx.options.healpix_forceinstall=True
#  #try:
#  #  atl.conf_lib(ctx,"chealpix",["chealpix","m","cfitsio"],"pix2vec_ring","chealpix.h",msg="or check that the path also point toward your cfitsio install",opt_name="healpix",uselib=["cshlib"])
#  #  atl.conf_lib(ctx,"healpix_f90",["healpix","cfitsio"],"HEALPIX_TYPES",msg="or check that the path also point toward your cfitsio install",opt_name="healpix",add_lib_code="f90",uselib=["fcshlib"])
#  #except:
#    #if not ctx.options.healpix_install:
#    #  raise
#    #else:
#    #Logs.pprint("PINK","healpix not found. Try to install it")
#    Logs.pprint("PINK","Try to install healpix it")
#    ctx.options.healpix_islocal=True
#    #atl.conf_lib(ctx,"cfitsio",["cfitsio"],"fits_init_cfitsio","fitsio.h",msg="cfitsio will be installed",opt_name="healpix",uselib=["cshlib"])
#    #if not bool(ctx.env.has_cfitsio):
#    install_cfitsio(ctx)
#    atl.conf_lib(ctx,"cfitsio",["cfitsio"],"fits_init_cfitsio","fitsio.h",msg="cfitsio will be installed",opt_name="healpix",uselib=["cshlib"])
#    if not bool(ctx.env.has_cfitsio):
#      raise Errors.WafError("Cannot build %s"%"cfitsio")
#    install_healpix(ctx)
#  atl.conf_lib(ctx,"cfitsio",["cfitsio"],"fits_init_cfitsio","fitsio.h",msg="",opt_name="healpix",uselib=["cshlib"])
#  atl.conf_lib(ctx,"chealpix",["chealpix","m","cfitsio"],"pix2vec_ring","chealpix.h",msg="or check that the path also point toward your cfitsio install",opt_name="healpix",uselib=["cpic"])
#  atl.conf_lib(ctx,"healpix_f90",["healpix","cfitsio"],"HEALPIX_TYPES",msg="or check that the path also point toward your cfitsio install",opt_name="healpix",add_lib_code="f90",uselib=["fpic"])
#      
#def install_cfitsio(ctx):
#  atl.installsmthg_pre(ctx,"ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/cfitsio3280.tar.gz","cfitsio3280.tar.gz")
#  CCMACRO = "\"%s %s\""%(ctx.env.CC[0],ctx.env.mopt)
#  CCMACRO = "CC=%s CXX=%s "%(CCMACRO,CCMACRO)
#  CPPMACRO = "CPP=\"%s -E\" CXXCPP=\"g++ -E\" "%(ctx.env.CC[0])
#  cmdline = "cd build/%s; ./configure --prefix=%s %s  %s %s; make clean;make -j %d ;make -j %d shared;make install"%("cfitsio",ctx.env.mprefix,"",CCMACRO, CPPMACRO,ctx.options.jobs,ctx.options.jobs)
#  Logs.pprint("PINK",cmdline)
#  if ctx.exec_command(cmdline)!=0:
#    raise Errors.WafError("Cannot build %s"%"cfitsio")
    
def install_healpix(ctx):
  if ctx.env.th:
    return 
  import os
  hpdir = "Healpix_2.20a"
  atl.installsmthg_pre(ctx,"http://sourceforge.net/projects/healpix/files/Healpix_2.20a/Healpix_2.20a_2011Feb09.tar.gz/download","Healpix_2.20a_2011Feb09.tar.gz")
  fpic_c = [vv for vv in ctx.env.CFLAGS_cshlib if "-fpic" in vv.lower()]
  fpic_f90 = [vv for vv in ctx.env.CFLAGS_cshlib if "-fpic" in vv.lower()]
  
  dii={"CC":ctx.env.CC[0],"CFLAGS":" ".join(ctx.env.CCFLAGS+fpic_c+[" -DDONT_USE_PSHT"]),"LIBDIR":ctx.env.LIBDIR,"INCDIR":ctx.env.PREFIX+"/include","FC":ctx.env.FC+" "+" ".join(ctx.env.FCFLAGS+fpic_f90+[" -DDONT_USE_PSHT"]),"FFLAGS":" "}
  # if I am here, I found cfitsio
  # could it be somewhere else ?
  cfitsiopath=""

  for pth in ctx.env.LIBPATH_cfitsio+["/usr/local/lib","/usr/lib","/lib"]+os.environ.get("LD_LIBRARY_PATH","").split(":")+os.environ.get("DYLD_LIBRARY_PATH","").split(":"):
    #print pth    
    #print osp.join(pth,ctx.env.cshlib_PATTERN%"cfitsio")
    if osp.exists(osp.join(pth,ctx.env.cshlib_PATTERN%"cfitsio")):
      cfitsiopath = pth
      break
  if not bool(cfitsiopath):
    raise Exception("cannot find cfitsio !")
  dii["CFITSIOPATH"]=cfitsiopath
  dii["CFITSIOPATHINC"]=osp.realpath(cfitsiopath+"/../include")
  #dii["CFITSIOPATH"]=ctx.env.LIBPATH_cfitsio[0]
  #dii["CFITSIOPATHINC"]=ctx.env.INCLUDES_cfitsio[0]
  #print dii

  f=open(osp.join("build",hpdir,"conf_1.cmd"),"w")
  print >>f,cnf_tmpl_1%dii
  f.close()
  # prepare a few things
  try:
    os.mkdir(osp.join("build",hpdir,"lib"))
  except Exception,e:
    #print e
    pass
  try:
    os.mkdir(osp.join("build",hpdir,"include"))
  except Exception,e:
    #print e
    pass

  import re
  # first part
  cmdline = "cd build/%s; ./configure <conf_1.cmd"%hpdir
  Logs.pprint("PINK",cmdline)
  try:
    rr = ctx.cmd_and_log(cmdline, output=waflib.Context.BOTH)
  except Exception,e:
    Logs.pprint("PINK","first pass failed. Keep going...")
    rr=("",)
  addN = ""
  if re.findall("The following line should be inserted into your home",rr[0]) :
    addN = "n\n"

  dii["addN"]=addN
  f=open(osp.join("build",hpdir,"conf.cmd"),"w")
  print >>f,cnf_tmpl%dii
  f.close()
  
  
  cmdline = "cd build/%s; ./configure <conf.cmd"%hpdir
  Logs.pprint("PINK",cmdline)
  try:
    rr = ctx.cmd_and_log(cmdline, output=waflib.Context.BOTH)
  except Exception,e:
    raise Errors.WafError("Healpix configure failed. Cannot build healpix")
  if re.findall("Something\s+went\s+wrong",rr[0]+rr[1]) :
    raise Errors.WafError("Healpix configure failed. Cannot build healpix")
  
  cmdline = "cd build/%s; make c-all f90-modules "%hpdir
  if ctx.exec_command(cmdline)!=0:
    raise Errors.WafError("Cannot build healpix. Look in build/config.log to understand why")
    
  import shutil,os
  for fi in os.listdir(osp.join("build",hpdir,"lib")):
    #print "copy",osp.join("build",hpdir,"lib",fi),osp.join(ctx.env.LIBDIR,fi)
    shutil.copyfile(osp.join("build",hpdir,"lib",fi),osp.join(ctx.env.LIBDIR,fi))
  for fi in os.listdir(osp.join("build",hpdir,"include")):
    #print "copy",osp.join("build",hpdir,"include",fi),osp.join(ctx.env.PREFIX,"include",fi)
    shutil.copyfile(osp.join("build",hpdir,"include",fi),osp.join(ctx.env.PREFIX,"include",fi))
  
  ctx.env.th=True    
cnf_tmpl_1="""6
n
0
"""
cnf_tmpl="""6
%(addN)s2
%(CC)s
-O2 -Wall %(CFLAGS)s


%(CFITSIOPATH)s
%(CFITSIOPATHINC)s
y
3
%(FC)s




%(CC)s
-O %(CFLAGS)s


%(CFITSIOPATH)s

0
0
"""
