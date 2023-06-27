def options(ctx):
  ctx.load("python")
  ctx.add_option("--h5py_install",action="store_true",default="",help="try to install h5py")
  
def configure(ctx):
  import waflib.Logs
  from waflib import Errors
  ctx.load("python")
  doit = False
  import sys
  sys.path+=[ctx.env.PYTHONDIR]
  
  if ctx.options.h5py_install or ctx.options.upgrade_all_deps or ctx.options.install_all_deps:
    try:
      import h5py
      raise Exception()
    except Exception,e: 
      doit=True
    if doit:
      import os.path as osp
      import autoinstall_lib as atl
      atl.installsmthg_pre(ctx,"http://h5py.googlecode.com/files/h5py-1.3.1.tar.gz","h5py-1.3.1.tar.gz")
      if ctx.env.INCLUDES_hdf5:
        HDF5_DIR=osp.split(ctx.env.INCLUDES_hdf5[0])[0]
      else:
        fi = ctx.find_file("hdf5.h",ctx.env.INCLUDES_pmc)
        print fi
        HDF5_DIR=osp.split(osp.split(fi)[0])[0]
      HDF5_API="18"
      print HDF5_DIR
      cmdline =  "cd build/%s; HDF5_DIR=%s HDF5_API=%s python setup.py install --install-purelib=%s"%("h5py-1.3.1",HDF5_DIR,HDF5_API,ctx.env.PYTHONDIR)
      waflib.Logs.pprint("PINK",cmdline)
      if ctx.cmd_and_log(cmdline)!=0:
        raise Errors.ConfigurationError("Cannot build h5py")
  try:
    import h5py
  except Exception,e: 
    if not doit:
      waflib.Logs.pprint("PINK", "You can install automatically h5py using cmdline option --h5py_install")      
    else:
      waflib.Logs.pprint("RED", "Autoinstall h5py has failed !")      
    raise e
    