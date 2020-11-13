from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from waflib import Logs
from waflib.Logs import colors
import sys
import os.path as osp
import re

clik_version = "XXX"
plc_version = "XXX"

sys.path+=["waf_tools"]
import autoinstall_lib as atl

from waflib.Configure import conf


def get_version(ctx):
  res = ctx.cmd_and_log("hg identify --id", output=waflib.Context.STDOUT, quiet=waflib.Context.BOTH)
  svnversion = res
  f=open("svnversion","w")
  print(svnversion, file=f)
  f.close()
  
def get_tag(ctx):
  res = ctx.cmd_and_log("hg tags", quiet=waflib.Context.BOTH)
  clik_v = None
  plc_v = None
  for r in res.split("\n"):
    if not r:
      continue
    v = (r.split()[0]).strip()
    if v.startswith("clik_") and not clik_v:
      clik_v = v[len("clik_"):]
    elif v.startswith("plc_") and not plc_v:
      plc_v = v[len("plc_"):]
  
  global clik_version,plc_version
  if clik_v:
    clik_version = clik_v
  if plc_v:
    plc_version = plc_v


def options(ctx):
  ctx.add_option('--forceinstall_all_deps',action='store_true',default=False,help='Install all dependencies',dest="install_all_deps")
  ctx.add_option('--install_all_deps',action='store_true',default=False,help='Install all dependencies (if they have not already been installed)',dest="upgrade_all_deps")
  
  ctx.load("local_install","waf_tools")
  ctx.load("try_icc","waf_tools")
  ctx.load("try_ifort","waf_tools")
  ctx.load("mbits","waf_tools")
  ctx.load("any_lapack","waf_tools")
  ctx.load("cfitsio","waf_tools")

  ctx.load("pmclib","waf_tools")
  ctx.load("python")

  try:
    import waflib.Configure
    if not osp.exists(osp.join(osp.split(waflib.__file__)[0],"extras/cython.py")):
      waflib.Configure.download_tool("cython",ctx=ctx)  
    ctx.load("cython",dowload=True)
  except Exception as e:
    pass
  
  import optparse
  grp=optparse.OptionGroup(ctx.parser,"Python options")
  grp.add_option('--nopyc',action='store_false',default=1,help='Do not install bytecode compiled .pyc files (configuration) [Default:install]',dest='pyc')
  grp.add_option('--nopyo',action='store_false',default=1,help='Do not install optimised compiled .pyo files (configuration) [Default:install]',dest='pyo')

  grp.add_option("--no_pytools",action="store_true",default=False,help="do not build the python tools")
  ctx.add_option_group(grp)
  options_numpy(ctx)
  options_pyfits(ctx)
  options_cython(ctx)
  
  grp=optparse.OptionGroup(ctx.parser,"Plugins options")
  #grp.add_option("--no_bopix",action="store_true",default=True,help="do not build bopix")
  #grp.add_option("--bopix",action="store_true",default=False,help="build bopix")
  grp.add_option("--no_lowlike",action="store_true",default=False,help="do not build lowlike")
  grp.add_option("--wmap_src",action="store",default="",help="location of wmap likelihood sources")
  grp.add_option("--wmap_7_install",action="store_true",default=False,help="download wmap 7 likelihood for me")
  grp.add_option("--wmap_9_install",action="store_true",default=False,help="download wmap 9 likelihood for me")
  grp.add_option("--wmap_install",action="store_true",default=False,help="download latest wmap likelihood for me")
  #grp.add_option("--wmap_dh_install",action="store_true",default=False,help="download D&H modified wmap likelihood for me")
  grp.add_option("--no_lenslike",action="store_true",default=False,help="do not build the lensing likelihood")
  
  ctx.add_option_group(grp)
  
  ctx.add_option("--extra_libpath",action="store",default="",help="libpath for extra lib to be linked")
  ctx.add_option("--extra_lib",action="store",default="",help="extra lib to be linked")
  ctx.add_option("--debug",action="store_true",default=False,help="if present compile with debug on")

def configure(ctx):
  import os
  import os.path as osp
  allgood = True

  try:
    # read extra env file
    f=open("clik_extra_env")
    for li in f:
      li = li.split("#")[0]
      if (not li.strip()):
        continue
      left,right = li.split("=")
      left = left.strip()
      right = right.strip()
      ctx.env[left] = right
  except Exception as e:
    pass

  try:
    ctx.load("try_icc","waf_tools")
  except Exception as e:
    Logs.pprint("RED","No suitable c compiler found (cause: '%s')"%e)
    ctx.fatal('The configuration failed') 
  ctx.load("mbits","waf_tools")
  ctx.load("osx_shlib","waf_tools")
  #ctx.load("c_openmp","waf_tools")
  #ctx.check_openmp_cflags(mandatory=False)

  try:
    ctx.load("try_ifort","waf_tools")
    ctx.env.has_f90 = True
  except Exception as e:
    Logs.pprint("RED","No suitable fortran compiler found (cause: '%s')"%e)
    ctx.fatal('The configuration failed') 
    ctx.env.has_f90 = False
    allgood = False
  ctx.load("local_install","waf_tools")
  if ctx.env.GCC_VERSION and ctx.env.IFORT_VERSION:
    ctx.start_msg("check compatibility of ifort and gcc versions")
    gcc_v = [int(v) for v in ctx.env.GCC_VERSION.split(".")]
    ifort_v = [int(v) for v in ctx.env.IFORT_VERSION.split(".")]
    if gcc_v[0]==4 and gcc_v[1]==9:
      if ifort_v[0]<14 or (ifort_v[0]==14 and ifort_v[1]==0 and ifort_v[2]<4):
        ctx.end_msg("%s can only be linked with ifort >14.0.4 (you have %s)"%( ctx.env.GCC_VERSION, ctx.env.IFORT_VERSION),color="RED")
        ctx.fatal('The configuration failed') 
    ctx.end_msg("ok")


  if not ctx.options.no_pytools:
    import sysconfig
    ctx.env.append_value("LIB_PYEMBED",['m','dl','util','pthread'])
    ctx.env.append_value("LIB_PYEXT",['m','dl','util','pthread'])
    ctx.env.append_value("LIBPATH_PYEMBED",[sysconfig.get_config_var("LIBDIR")])
    ctx.env.append_value("LIBPATH_PYEXT",[sysconfig.get_config_var("LIBDIR")])
    ctx.load("python")
    
    if ctx.env.PYTHON[0]!=sys.executable:
      from waflib.Logs import warn
      warn("reverting to current executable")
      ctx.env.PYTHON[0]=sys.executable
      os.environ["PATH"]=":".join(set(os.environ["PATH"].split(":")+[osp.dirname(sys.executable)]))
    try:
      # remove unwanted flags for darwin
      ctx.env.FRAMEWORK_ST=""
      ctx.check_python_version()
      #ctx.env["PYTHONDIR"]=ctx.get_python_variables(["get_python_lib(standard_lib=0, prefix=%r) or ''"%ctx.env['PREFIX']])[0]
      ctx.env["PYTHONDIR"]="%s/lib/python/site-packages"%ctx.env['PREFIX']
      #print(ctx.env["PYTHONDIR"])
      ctx.check_python_headers("pyembed")
      ctx.env.INCLUDES_PYEXT = ctx.env.INCLUDES_PYEMBED
      # remove unwanted flags for darwin
      _remove_arch(ctx,"CFLAGS_PYEXT")
      _remove_arch(ctx,"LINKFLAGS_PYEXT")
      from distutils.sysconfig import get_config_var
      for v in ctx.env["DEFINES"]:
        if "PYTHON" in v:
          ctx.undefine(v)
    except Exception as e:
      #ctx.options.no_pytools = True
      Logs.pprint("BLUE","No suitable python distribution found")
      Logs.pprint("BLUE","Cause : '%s'"%e)
      Logs.pprint("BLUE","Compilation will continue without it (but I strongly advise that you install it)")
      allgood = False
    
    
  # dl
  ctx.check_cc(lib="dl",mandatory=1,uselib_store="dl")
  ctx.check_cc(lib="dl",mandatory=0,defines=["HAS_RTLD_DEFAULT"],fragment="#include <dlfcn.h> \nint main() {void* tt = RTLD_DEFAULT;}",msg="checking for RTLD_DEFAULT in dl",uselib_store="dl")
  
  # rpath
  ctx.env.append_value("RPATH",ctx.env.PREFIX+"/lib")

  # debug
  ctx.start_msg("compile with debug option")
  if ctx.options.debug:
    ctx.env.FCFLAGS_fpic.append("-g")
    ctx.env.FCFLAGS.append("-g")
    ctx.env.CFLAGS_fpic.append("-g")
    ctx.env.CFLAGS.append("-g")
    ctx.end_msg("ON")
  else:
    ctx.end_msg("OFF")
  
  #configure pmc
  ctx.env.has_pmc = False
  ctx.env.silent_pmc = True
  ctx.load("pmclib","waf_tools")

  
  if (not ctx.env.has_pmc) or ("HAS_LAPACK" not in ctx.env.DEFINES_pmc):
    #configure lapack
    ctx.env.has_lapack = True
    ctx.load("any_lapack","waf_tools")

  
  #configure cfitsio
  ctx.env.has_cfitsio = True
  ctx.load("cfitsio","waf_tools")

  #lenslike
  ctx.env.has_lenslike = not (ctx.options.no_lenslike or not osp.exists("src/lenslike"))
  
  #bopix
  #ctx.env.has_bopix = not ((not ctx.env.bopix) or ctx.options.no_bopix or not osp.exists("src/bopix"))
  
  #lowlike
  ctx.env.has_lowlike =   not (ctx.options.no_lowlike or not osp.exists("src/lowlike"))

  
  #camspec
  ctx.env.has_camspec = osp.exists("src/camspec")
  ctx.env.has_camspec_v3 = osp.exists("src/camspec/temp_like_v3.f90")

  #cmbonly
  ctx.env.has_cmbonly = osp.exists("src/cmbonly")

  #actspt
  ctx.env.has_actspt = osp.exists("src/actspt")
  #ctx.env.has_actspt = False

  #gibbs
  ctx.env.has_gibbs = osp.exists("src/gibbs")

  #egfs 
  #ctx.env.has_egfs = osp.exists("src/egfs")
  #ctx.env.has_egfs = False

  #bflike
  ctx.env.has_bflike = osp.exists("src/bflike")

  #bicep
  #ctx.env.has_bicep = osp.exists("src/bicep") and False

  #mspec
  ctx.env.has_mspec = osp.exists("src/mspec")

  #momento
  #ctx.env.has_momento = osp.exists("src/momento")

  #lollipop
  #ctx.env.has_lollipop = osp.exists("src/lollipop")

  # wmap
  if (ctx.options.wmap_7_install or ctx.options.wmap_9_install or ctx.options.wmap_install) and not ctx.options.wmap_src :
    #if (ctx.options.wmap_7_install or ctx.options.wmap_9_install or ctx.options.wmap_install) and not ctx.options.wmap_src :
    if ctx.options.wmap_7_install:
      atl.installsmthg_pre(ctx,"http://lambda.gsfc.nasa.gov/data/map/dr4/dcp/wmap_likelihood_sw_v4p1.tar.gz","wmap_likelihood_sw_v4p1.tar.gz","src/")
      ctx.options.wmap_src = "likelihood_v4p1" 
    else:
      atl.installsmthg_pre(ctx,"http://lambda.gsfc.nasa.gov/data/map/dr5/dcp/likelihood/wmap_likelihood_sw_v5.tar.gz","wmap_likelihood_sw_v5.tar.gz","src/")
      ctx.options.wmap_src = "wmap_likelihood_v5"   
    ctx.env.wmap_src = ctx.options.wmap_src

  if ctx.options.wmap_src:
    ld = os.listdir(osp.join("src",ctx.options.wmap_src))
    for v in ld:
      if "7yr" in v:
        ctx.env.wmap_version = 7
        break
      elif "9yr" in v:
        ctx.env.wmap_version = 9
        break
  
  if not ctx.options.no_pytools:
    # check_python version
    ctx.start_msg("Python version allows for automatic install")
    vr = sys.version_info[0:3]
    ver = ".".join(str(v) for v in sys.version_info[0:3])
    if vr[0]==2 and ((vr[1]<7) or (vr[1]==7 and vr[2]<9)):
      ctx.end_msg("no (%s<2.7.9, the extra python packages will need to be manually installed)"%ver,"YELLOW")
    else:
      ctx.end_msg("yes (%s)"%ver)
    try:
      configure_numpy(ctx)
      configure_pyfits(ctx)
      configure_cython(ctx)
      
    except Exception as e:
      #ctx.options.no_pytools = True
      Logs.pprint("BLUE","No suitable python distribution found")
      Logs.pprint("BLUE","Cause : '%s'"%e)
      Logs.pprint("BLUE","Compilation will continue without it (but I strongly advise that you install it)")
      allgood = False
  
  ctx.env.has_plik = False
  if osp.exists("src/plik"):
    ctx.env.has_plik = True
  ctx.env.plik_release = False
  if not osp.exists("src/plik/smica_ext.c"):
    ctx.env.plik_release=True
  # go through the component plugins
  if osp.exists("src/plik/component_plugin"):
    for plg in os.listdir("src/plik/component_plugin"):
      if plg[0]==".":
        continue
      if plg=="rel2015" and ctx.env.plik_release==False:
        continue
      pth = osp.join("src/plik/component_plugin",plg)
      if not osp.isdir(pth):
        continue
      ctx.start_msg("Add plugin : '%s'"%plg)
      try:
        f = open(osp.join(pth,"plugin.txt"))
        txt = f.read()
        decr = {"source":"","data":"","python":""}
        decr.update(dict(re.findall("(.+?)\s+(?::|=)\s+(.+)",txt)))
        for k in decr:
          decr[k] = [l.strip() for l in re.split("(?:\s|,|;)+",decr[k]) if l.strip()]
          for f in decr[k]:
            if not osp.exists(osp.join(pth,f)):
              raise Exception((osp.join(pth,f))+" does not exists")
      except Exception as e:
        ctx.end_msg("Ignored (%s)"%e,"YELLOW")
        continue
      ctx.end_msg("ok")
      ctx.env.append_unique("PLG",[plg])
      ctx.env.append_unique("PLG_%s_SRC"%plg,decr["source"])
      ctx.env.append_unique("PLG_%s_DATA"%plg,decr["data"])
      ctx.env.append_unique("PLG_%s_PYTHON"%plg,decr["python"])


  try:
    f=open("svnversion")
  except :
    get_version(ctx)
    f=open("svnversion")
  ctx.env.svnversion = f.read()
  f.close()
  
  # extra libs
  if ctx.options.extra_lib:
    libpath = ctx.options.extra_libpath.split(":")
    lib = ctx.options.extra_lib.split(",")
    ctx.env.LIBPATH_extra = libpath
    ctx.env.LIB_extra = lib
    
  if allgood==False:
    print("\n%s*----------------------------------------------------*"%(colors.PINK+colors.BOLD))
    print("%s|%s                                                    %s|"%(colors.PINK+colors.BOLD,colors.NORMAL,colors.PINK+colors.BOLD))
    print("%s|%s   Configure partial.                               %s|"%(colors.PINK+colors.BOLD,colors.NORMAL,colors.PINK+colors.BOLD))
    print("%s|%s   You can build now but you will be lacking some   %s|"%(colors.PINK+colors.BOLD,colors.NORMAL,colors.PINK+colors.BOLD))
    print("%s|%s   features.                                        %s|"%(colors.PINK+colors.BOLD,colors.NORMAL,colors.PINK+colors.BOLD))
    print("%s|%s   run './waf install' now !                        %s|"%(colors.PINK+colors.BOLD,colors.NORMAL,colors.PINK+colors.BOLD))
    print("%s|%s                                                    %s|"%(colors.PINK+colors.BOLD,colors.NORMAL,colors.PINK+colors.BOLD))
    print("%s*----------------------------------------------------*%s\n"%(colors.PINK+colors.BOLD,colors.NORMAL))
    
  else:
    print("\n%s*----------------------------------------------------*"%(colors.PINK+colors.BOLD))
    print("%s|%s                                                    %s|"%(colors.PINK+colors.BOLD,colors.NORMAL,colors.PINK+colors.BOLD))
    print("%s|%s   Configure step ok.                               %s|"%(colors.PINK+colors.BOLD,colors.NORMAL,colors.PINK+colors.BOLD))
    print("%s|%s                                                    %s|"%(colors.PINK+colors.BOLD,colors.NORMAL,colors.PINK+colors.BOLD))
    print("%s|%s   run './waf install' now !                        %s|"%(colors.PINK+colors.BOLD,colors.NORMAL,colors.PINK+colors.BOLD))
    print("%s|%s                                                    %s|"%(colors.PINK+colors.BOLD,colors.NORMAL,colors.PINK+colors.BOLD))
    print("%s*----------------------------------------------------*%s\n"%(colors.PINK+colors.BOLD,colors.NORMAL))
    
  
def build(ctx):
  import os
    
  ctx.recurse("src")
  if not ctx.options.no_pytools:
    ctx.recurse("src/python")
  
  if not osp.exists("%s/share"%ctx.env["PREFIX"]):
    os.mkdir("%s/share"%ctx.env["PREFIX"])
  if not osp.exists("%s/share/clik"%ctx.env["PREFIX"]):
    os.mkdir("%s/share/clik"%ctx.env["PREFIX"])

  # install data
  #if ctx.env.has_egfs:
  #  ctx.install_files('${PREFIX}/share/clik/egfs', 
  #                  'src/egfs/egfs_data/clustered_150.dat src/egfs/egfs_data/clustered_flat.dat src/egfs/egfs_data/ksz_ov.dat src/egfs/egfs_data/ksz_patchy.dat src/egfs/egfs_data/tsz.dat src/egfs/egfs_data/clustered_1108.4614.dat')

  for plg in ctx.env.PLG:
    data = getattr(ctx.env,"PLG_%s_DATA"%plg)
    ppp = "%s/share/clik/%s"%(ctx.env["PREFIX"],plg)
    if not osp.exists(ppp):
      os.mkdir("%s/share/clik/%s"%(ctx.env["PREFIX"],plg))
    if data:
      ctx.install_files('${PREFIX}/share/clik/%s'%plg,[osp.join('src/plik/component_plugin/%s'%plg,dt) for dt in data])

  ctx.add_post_fun(post)

def _remove_arch(ctx,evn):
  if sys.platform.lower()=="darwin":
    cflags_pyext = getattr(ctx.env,evn)
    cflags_pyext_new = []
    inarch = 0
    for cf in cflags_pyext:
      #print cf,inarch
      if inarch == 1:
        inarch = 0
        continue
      if cf == "-arch":
        inarch = 1
        continue
      cflags_pyext_new +=[cf]
    setattr(ctx.env,evn,cflags_pyext_new)


def dist(ctx):
  print("private")
  import re
  get_tag(ctx)
  ctx.base_name = 'clik-'+clik_version
  get_version(ctx)
  dist_list =  "svnversion waf wscript **/wscript src/minipmc/* src/cldf/* waf_tools/*.py waf_tools/*.txt readme.md"
  dist_list += "src/* src/python/clik/*.py src/python/clik/*.pxd src/python/clik/*.pyx "
  dist_list += "src/python/tools/*.py "
  dist_list += "examples/*.par examples/*.dat "
  dist_list += "src/plik/component_plugin/** src/plik/* "
  dist_list += "src/* src/camspec/* "
  dist_list += "src/actspt/* "
  dist_list += "src/lowlike/* "
  dist_list += "src/gibbs/* "
  dist_list += "src/mspec/* "
  #dist_list += "src/lollipop/* "
  dist_list += "src/bflike/* "
  dist_list += "src/cmbonly/* "
  #dist_list += "src/ftau/* "
  #dist_list += "src/simbal/* "
  #dist_list += "src/momento/* "
  dist_list += "src/simall/* "
  #dist_list += "src/bicep/* "
  dist_list += "src/lenslike/plenslike/*.c src/lenslike/plenslike/*.h "
  
  exclude_list = ["src/plik/component_plugin/rel2015/* ","src/plik/*_rl.*"]
  excl_list = ctx.path.ant_glob(exclude_list)
  
  #print ctx.path.ant_glob(dist_list)
  excl_list = ctx.path.ant_glob(exclude_list)
  files = ctx.path.ant_glob(dist_list)
  giles = []
  print(excl_list)
  for f in files:
    if f in excl_list:
      continue
    giles +=[f]
  ctx.files = giles

import waflib
class Dist_public(waflib.Scripting.Dist):
  cmd = 'dist_public'
  fun = 'dist_public'
  
def dist_public(ctx):
  print("public")
  import re
  get_tag(ctx)
  ctx.base_name = 'plc-'+plc_version
  get_version(ctx)
  dist_list =  "Makefile setup.py svnversion readme.md waf wscript **/wscript src/minipmc/* src/cldf/* waf_tools/*.py  waf_tools/*.txt "
  dist_list += "src/python/clik/*.py src/python/clik/*.pxd src/python/clik/*.pyx "
  #dist_list += "src/python/tools/*.py "
  dist_list += "src/* src/camspec/* "
  dist_list += "src/actspt/* "
  dist_list += "src/lowlike/* "
  dist_list += "src/gibbs/* "
  dist_list += "src/cmbonly/* "
  dist_list += "src/lowlike/* "
  dist_list += "src/bflike/* "
  dist_list += "src/simall/* "
  
  dist_list += " src/plik/component_plugin/rel2015/* src/plik/* "
  dist_list += "src/lenslike/plenslike/*.c src/lenslike/plenslike/*.h "
  dist_list+=" src/python/tools/".join(["clik_cldf_dump","clik_cldf_ls","clik_explore_1d.py",
              "clik_get_selfcheck.py",
              "clik_example_py.py",
              "clik_print.py",
              "prepare_wmap.py",
              "clik_change_lrange.py", "clik_change_lrange_bychannel_pol.py"])
  
  exclude_list = []
  exclude_list += ["src/actspt/test.f90","src/actspt/test_actspt.f90"]
  exclude_list += ["src/bflike/test_bflike.f90","src/bflike/test_bflike_smw.f90"]
  exclude_list += ["src/camspec/CAMtst.f90"]
  exclude_list += ["src/cldf/test_cldf.c"]
  exclude_list += ["src/cmbonly/plik_cmbonly_test.f90"]
  exclude_list += ["src/gibbs/test_comm.c","src/gibbs/validate_comm_lowl.c",]
  exclude_list += ["src/lowlike/test.F90"]
  exclude_list += ["src/plik/*_ext.*"]
  
  
  excl_list = ctx.path.ant_glob(exclude_list)
  files = ctx.path.ant_glob(dist_list)
  giles = []
  for f in files:
    #print f.abspath(),
    if f in excl_list:
      #print "X"
      continue
    #print ""
    giles +=[f]
  ctx.files = giles
  
class Dist_public_plc2py3(waflib.Scripting.Dist):
  cmd = 'dist_public_plc2py3'
  fun = 'dist_public_plc2py3'
  
def dist_public_plc2py3(ctx):
  print("public")
  import re
  get_tag(ctx)
  ctx.base_name = 'plc-'+plc_version
  get_version(ctx)
  dist_list =  "Makefile setup.py svnversion readme.md waf wscript **/wscript src/minipmc/* src/cldf/* waf_tools/*.py  waf_tools/*.txt "
  dist_list += "src/python/clik/*.py src/python/clik/*.pxd src/python/clik/*.pyx "
  #dist_list += "src/python/tools/*.py "
  dist_list += "src/* src/camspec/* "
  dist_list += "src/actspt/* "
  dist_list += "src/lowlike/* "
  dist_list += "src/gibbs/* "
  dist_list += "src/cmbonly/* "
  dist_list += "src/lowlike/* "
  dist_list += "src/bflike/* "
  #dist_list += "src/simall/* "
  
  dist_list += " src/plik/component_plugin/rel2015/* src/plik/* "
  dist_list += "src/lenslike/plenslike/*.c src/lenslike/plenslike/*.h "
  dist_list+=" src/python/tools/".join(["clik_cldf_dump","clik_cldf_ls","clik_explore_1d.py",
              "clik_get_selfcheck.py",
              "clik_example_py.py",
              "clik_print.py",
              "prepare_wmap.py",
              "clik_change_lrange.py"])
  
  exclude_list = []
  exclude_list += ["src/actspt/test.f90","src/actspt/test_actspt.f90"]
  exclude_list += ["src/bflike/test_bflike.f90","src/bflike/test_bflike_smw.f90"]
  exclude_list += ["src/camspec/CAMtst.f90","temp_like.f90", "temp_like_v3.f90.c"]
  exclude_list += ["src/cldf/test_cldf.c"]
  exclude_list += ["src/cmbonly/plik_cmbonly_test.f90"]
  exclude_list += ["src/gibbs/test_comm.c","src/gibbs/validate_comm_lowl.c","comm_gauss_br_mod_v3.f90"]
  exclude_list += ["src/lowlike/test.F90"]
  exclude_list += ["src/plik/*_ext.*"]
  
  
  excl_list = ctx.path.ant_glob(exclude_list)
  files = ctx.path.ant_glob(dist_list)
  giles = []
  for f in files:
    #print f.abspath(),
    if f in excl_list:
      #print "X"
      continue
    #print ""
    giles +=[f]
  ctx.files = giles
def post(ctx):
  import shutil
  from waflib import Utils
  import os
  if ctx.cmd == 'install':
    # install the module file. This is a cheap trick... grml
    shutil.copy('build/clik.mod',ctx.env.INCDIR)
    shutil.copy('build/clik_plik.mod',ctx.env.INCDIR)
    # go around a waf bug which set the wrong chmod to fortran exec
    os.chmod("%s/bin/clik_example_f90"%ctx.env["PREFIX"],Utils.O755)
    build_env_files(ctx)
    
def build_env_files(ctx):
  import os
  import os.path as ops
  #full_libpath = set(ctx.env.LIBPATH_fc_runtime + ctx.env.LIBPATH_lapack)
  full_libpath = set(ctx.env.LIBPATH_fc_runtime  + ctx.env.LIBPATH_lapack + ctx.env.LIBPATH_cfitsio + [ctx.env.PREFIX+"/lib"])
  #print full_libpath
  #tcsh and co
  shell = "csh"
  name = "clik_profile.csh"
  shebang = "#! /bin/tcsh"
  extra=""
  single_tmpl = "setenv %(VAR)s %(PATH)s\n"
  block_tmpl = """
if !($?%(VAR)s) then
setenv %(VAR)s %(PATH)s
else
setenv %(VAR)s %(PATH)s:${%(VAR)s}
endif
"""
  multi_tmpl = """if !($?%(VAR)s) then
setenv %(VAR)s %(PATH)s
else
set newvar=$%(VAR)s
set newvar=`echo ${newvar} | sed s@:%(PATH)s:@:@g`
set newvar=`echo ${newvar} | sed s@:%(PATH)s\$@@` 
set newvar=`echo ${newvar} | sed s@^%(PATH)s:@@`  
set newvar=%(PATH)s:${newvar}                     
setenv %(VAR)s %(PATH)s:${newvar} 
endif"""
  
  __dofile(ctx,name,shell,extra,multi_tmpl,single_tmpl,full_libpath)

  #bash and co
  shell = "sh"
  name = "clik_profile.sh"
  shebang = "#! /bin/sh"
  extra= """function addvar () {
local tmp="${!1}" ;
tmp="${tmp//:${2}:/:}" ; tmp="${tmp/#${2}:/}" ; tmp="${tmp/%:${2}/}" ;
export $1="${2}:${tmp}" ;
}"""
  single_tmpl = "%(VAR)s=%(PATH)s\nexport %(VAR)s\n"
  block_tmpl = """
if [ -z "${%(VAR)s}" ]; then 
%(VAR)s=%(PATH)s
else
%(VAR)s=%(PATH)s:${%(VAR)s}
fi
export %(VAR)s
"""
  multi_tmpl = """if [ -z "${%(VAR)s}" ]; then 
%(VAR)s=%(PATH)s
export %(VAR)s
else
addvar %(VAR)s %(PATH)s
fi"""
    
  __dofile(ctx,name,shell,extra,multi_tmpl,single_tmpl,full_libpath)
  
  print("\n%s*----------------------------------------------------*"%(colors.PINK+colors.BOLD))
  print("%s|%s                                                    %s|"%(colors.PINK+colors.BOLD,colors.NORMAL,colors.PINK+colors.BOLD))
  print("%s|%s   Source clik_profile.sh (or clik_profile.csh)     %s|"%(colors.PINK+colors.BOLD,colors.NORMAL,colors.PINK+colors.BOLD))
  print("%s|%s   to set the environment variables needed by clik  %s|"%(colors.PINK+colors.BOLD,colors.NORMAL,colors.PINK+colors.BOLD))
  print("%s|%s                                                    %s|"%(colors.PINK+colors.BOLD,colors.NORMAL,colors.PINK+colors.BOLD))
  print("%s*----------------------------------------------------*%s\n"%(colors.PINK+colors.BOLD,colors.NORMAL))
  
def __dofile(ctx,name,shell,extra,multi_tmpl,single_tmpl,full_libpath):
  import sys
  if sys.platform.lower()=="darwin":
    LD_LIB = "DYLD_LIBRARY_PATH"
  else:
    LD_LIB = "LD_LIBRARY_PATH"
  f = open(osp.join(ctx.env.BINDIR,name),"w")
  print("# this code cannot be run directly", file=f)
  print("# do 'source %s' from your %s shell or put it in your profile\n"%(osp.join(ctx.env.BINDIR,name),shell), file=f)
  print(extra,"\n", file=f)
  print(multi_tmpl%{"PATH":ctx.env.BINDIR,"VAR":"PATH"}, file=f)
  print(multi_tmpl%{"PATH":ctx.env.PYTHONDIR,"VAR":"PYTHONPATH"}, file=f)
  for pt in full_libpath:
    print(multi_tmpl%{"PATH":pt,"VAR":LD_LIB}, file=f)
  print(single_tmpl%{"PATH":ctx.env.PREFIX,"VAR":"CLIK_PATH"}, file=f)
  print(single_tmpl%{"PATH":osp.join(ctx.env.PREFIX,"share/clik"),"VAR":"CLIK_DATA"}, file=f)
  print(single_tmpl%{"PATH":",".join(ctx.env.PLG),"VAR":"CLIK_PLUGIN"}, file=f)
  f.close()

def options_pyfits(ctx):
  import autoinstall_lib as atl
  atl.add_python_option(ctx,"pyfits")
def options_numpy(ctx):
  atl.add_python_option(ctx,"numpy")
def options_cython(ctx):
  atl.add_python_option(ctx,"cython")


def configure_numpy(ctx):
  import autoinstall_lib as atl
  url,tar = atl.get_lib_url(ctx,"numpy",["http://sourceforge.net/projects/numpy/files/NumPy/1.6.0/numpy-1.6.0.tar.gz/download","numpy-1.6.0.tar.gz"])
  atl.configure_python_module(ctx,"numpy",url,tar,"numpy")
  import imp
  numpy = imp.load_module("numpy",*imp.find_module("numpy"))
  ctx.env.append_value("INCLUDES_PYEXT",numpy.get_include())

def configure_pyfits(ctx):
  import autoinstall_lib as atl
  try:
    atl.check_python_module(ctx,"astropy")
  except ImportError as e:
    Logs.pprint("PINK","Astropy not found, try pyfits instead")
    url,tar = atl.get_lib_url(ctx,"pyfits",["http://pypi.python.org/packages/source/p/pyfits/pyfits-3.2.2.tar.gz","pyfits-3.2.2.tar.gz"])
    atl.configure_python_module(ctx,"pyfits",url,tar,"pyfits")
  
def configure_cython(ctx):
  import autoinstall_lib as atl
  from waflib import Utils
  import os.path as osp
  import os
  
  os.environ["PATH"] = ctx.env["PREFIX"][0]+"/bin:"+os.environ["PATH"]+":"+osp.dirname(osp.realpath(" ".join(ctx.env.PYTHON)))
  vv=False
  def postinstallcython():
    ctx.env.CYTHON=[osp.join(ctx.env.BINDIR,"cython")]
    f=open(osp.join(ctx.env.BINDIR,"cython"))
    cytxt = f.readlines()
    f.close()
    cytxt[1:1] = ["import sys\n","sys.path+=['%s']\n"%str(ctx.env.PYTHONDIR)]
    f=open(osp.join(ctx.env.BINDIR,"cython"),"w")
    f.write("".join(cytxt))
    f.close()
    os.chmod(osp.join(ctx.env.BINDIR,"cython"),Utils.O755)

  url,tar = atl.get_lib_url(ctx,"cython",["https://pypi.python.org/packages/b7/67/7e2a817f9e9c773ee3995c1e15204f5d01c8da71882016cac10342ef031b/Cython-0.25.2.tar.gz#md5=642c81285e1bb833b14ab3f439964086","Cython-0.25.2.tar.gz"])
  atl.configure_python_module(ctx,"cython",url,tar,"Cython",postinstall=postinstallcython)

  try:
    # check for cython
    atl.check_python_module(ctx,"cython")
    vv=True
    version_str = "unknown"
    mincyver = 20
    if int(ctx.env["PYTHON_VERSION"].split(".")[0])>2 and int(ctx.env["PYTHON_VERSION"].split(".")[1])>6:
      mincyver = 28
    ctx.start_msg("Checking cython version (>0.%s)"%mincyver)
    import Cython.Compiler.Version
    version_str = Cython.Compiler.Version.version
    version = [int(v) for v in version_str.split(".")]
    #print version
    assert version[0]>1 or version[1]>=mincyver
    ctx.end_msg(version_str)
  except Exception as e:
    if vv:
      ctx.end_msg("no (%s)"%version_str,'YELLOW')
    # no cython, install it !
    url,tar = atl.get_lib_url(ctx,"cython",["https://pypi.python.org/packages/b7/67/7e2a817f9e9c773ee3995c1e15204f5d01c8da71882016cac10342ef031b/Cython-0.25.2.tar.gz#md5=642c81285e1bb833b14ab3f439964086","Cython-0.25.2.tar.gz"])
    atl.configure_python_module(ctx,"cython",url,tar,"Cython",postinstall=postinstallcython)

  try:
    ctx.load("cython")
    
  except:
    ctx.env.CYTHON=[osp.join(ctx.env.BINDIR,"cython")]
    f=open(osp.join(ctx.env.BINDIR,"cython"))
    cytxt = f.readlines()
    f.close()
    cytxt[1:1] = ["import sys\n","sys.path+=['%s']\n"%str(ctx.env.PYTHONDIR)]
    f=open(osp.join(ctx.env.BINDIR,"cython"),"w")
    f.write("".join(cytxt))
    f.close()
    os.chmod(osp.join(ctx.env.BINDIR,"cython"),Utils.O755)
  
