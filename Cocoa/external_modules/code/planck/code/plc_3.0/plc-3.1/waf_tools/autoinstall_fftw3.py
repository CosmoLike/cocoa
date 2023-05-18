from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import autoinstall_lib as atl

version = "fftw-3.2.2"
tool = "fftw3"

print("-> loading %s autoinstall (using version %s)"%(tool,version))

def options(opt):
  atl.add_lib_option(tool,opt,install=True)
  
def configure(ctx):
  atl.conf_lib(ctx,tool,["fftw3"],"fftw_execute","fftw3.h",defines="HAS_FFTW3",install=installfftw3)
    
def installfftw3(ctx):
  filen = version+".tar.gz"
  atl.installsmthg_pre(ctx,"http://www.fftw.org/"+filen,filen)
  atl.installsmthg_post(ctx,filen,"fftw","--enable-shared")
