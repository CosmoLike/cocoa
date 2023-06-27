from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import autoinstall_lib as atl
from waflib import Logs
from waflib import Utils,Errors
import waflib

import os.path as osp
    
def options(ctx):
  atl.add_lib_option("cfitsio",ctx,install=True)
  
twice = False
def configure(ctx):
  atl.conf_lib(ctx,"cfitsio",["cfitsio"],"fits_init_cfitsio","fitsio.h",msg="",opt_name="cfitsio",uselib=["cshlib"],install=install_cfitsio)
  
def install_cfitsio(ctx):
  #url,tar = atl.get_lib_url(ctx,"cfitsio",("ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/cfitsio3280.tar.gz","cfitsio3280.tar.gz"))
  url,tar = atl.get_lib_url(ctx,"cfitsio",("https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio_latest.tar.gz","cfitsio_latest.tar.gz"))
  
  dirname = atl.installsmthg_pre(ctx,url,tar)
  CCMACRO = "\"%s %s\""%(ctx.env.CC[0]," ".join(ctx.env.mopt))
  CCMACRO = "CC=%s CXX=%s "%(CCMACRO,CCMACRO)
  CPPMACRO = "CPP=\"%s -E\" CXXCPP=\"g++ -E\" "%(ctx.env.CC[0])
  cmdline = "cd build/%s; ./configure --prefix=%s %s  %s %s; make clean;make -j %d ;make -j %d shared;make install"%(dirname,ctx.env.mprefix,"",CCMACRO, CPPMACRO,ctx.options.jobs,ctx.options.jobs)
  Logs.pprint("PINK",cmdline)
  if ctx.exec_command(cmdline)!=0:
    raise Errors.WafError("Cannot build %s"%"cfitsio")
    
