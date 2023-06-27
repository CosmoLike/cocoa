from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from waflib import Options
from waflib import Logs
from waflib import  Context
from waflib import Errors
import re  
import waflib
import os

def options(ctx):
  import optparse
  grp = ctx.parser.get_option_group("--gfortran")
  if grp==None:
    grp=optparse.OptionGroup(ctx.parser,"compiler options")
  grp.add_option("--gcc",action="store_true",default=False,help="Do not test for icc and only use gcc")
  grp.add_option("--icc",action="store_true",default=False,help="Do not test for gcc and only use icc")
  grp.add_option("--clang",action="store_true",default=False,help="Do not test for gcc and only use clang")
  ctx.add_option_group(grp)  

def show_linkline(ctx):
  ctx.start_msg("icc link line")
  ctx.end_msg(" ".join(["-L%s"%vv for vv in ctx.env.LIBPATH_cc_runtime])+" "+" ".join(["-l%s"%vv for vv in ctx.env.LIB_cc_runtime]))

def do_icc(ctx):
  ctx.env.CC=[]  
  ctx.env.LINK_CC=[]
  ctx.load('icc')
  ctx.check_cc(
      errmsg="failed",msg="Compile a test code with icc",
      mandatory=1,fragment = "#include <stdio.h>\nmain() {fprintf(stderr,\"hello world\");}\n",compile_filename='test.c',features='c cprogram')
  ctx.start_msg("retrieve icc link line")
  ctx.env["CCFLAGS_cc_omp"]=[]
  if int(ctx.env["CC_VERSION"][0])>15:
    ctx.env.append_value("CCFLAGS_cc_omp","-qopenmp")  
  else:
    ctx.env.append_value("CCFLAGS_cc_omp","-openmp")
  try:
    #print "%s %s -dryrun -dynamiclib -shared-intel -no-cxxlib dummy.f90"%(ctx.env.FC," ".join(ctx.env.FCFLAGS))
    llgo,llge = ctx.cmd_and_log("%s %s -dryrun -dynamiclib -shared-intel -no-cxxlib dummy.c"%(" ".join(ctx.env.CC)," ".join(ctx.env.CCFLAGS+ctx.env.CCFLAGS_cc_omp)), output=waflib.Context.BOTH)
    #print "RET",llgo,llge
    L = set([ll.strip() for ll in re.findall("-L(.+)\s*\\\\", llge.split("ld ")[1]) if ("icc" in ll.lower()) or ("intel" in ll.lower())])
    l = set([ll.strip() for ll in re.findall("-l(.+)\s*\\\\", llge.split("ld ")[1])])
    rL = set()
    rl = set()
    for Li in L:
      oli = os.listdir(Li)
      for li in l:
        if ctx.env.cshlib_PATTERN%li in oli:
          rl.add(li)
          rL.add(Li)
  except:
    ctx.end_msg(False)
    raise
  for pth in list(rL) + ["/lib","/lib64"]:
    ctx.env.append_value("LIBPATH",pth)
    ctx.env.append_value("RPATH",pth)
  
  ctx.env.append_value("LIB",list(rl)+["pthread"])
  ctx.end_msg(True)
  show_linkline(ctx)
  ctx.env.SHLIB_MARKER = ""
  ctx.env.STLIB_MARKER = ""

  ctx.env.has_icc = True

def do_clang(ctx):
  ctx.env.CC=[]
  ctx.env.LINK_CC=[]
  ctx.load('clang')
  ctx.start_msg("Check clang version") 
  v90 = ctx.cmd_and_log(" ".join(ctx.env.CC)+" --version",quiet=Context.STDOUT).split("\n")[0].strip()
  #version90 = re.findall("(4\.[0-9]\.[0-9])",v90)
  #if len(version90)<1:
    #Logs.pprint("PINK","Can't get gfortran version... Let's hope for the best")
  #  ctx.end_msg("not found, let's hope for the best...",color="PINK")
  #else:
  #  version90 = version90[0]
  #  vmid = int(version90.split(".")[1])
  #  if vmid<2:
  #    ctx.end_msg(v90,color="YELLOW")
  #    raise Errors.WafError("gcc version need to be above 4.2 got %s"%version90)
  ctx.end_msg(v90)
  ctx.check_cc(
    errmsg="failed",msg="Compile a test code with clang",
    mandatory=1,fragment = "#include <stdio.h>\nmain() {fprintf(stderr,\"hello world\");}\n",compile_filename='test.c',features='c cprogram')
  ctx.env["CCFLAGS_cc_omp"]=[]
  # try if omp is installed somewhere !
  try:
    ctx.env["CCFLAGS_cc_omp"]=['-Xpreprocessor','-fopenmp','-lomp']
    ctx.check_cc(lib="omp", libpath = "/usr/local/lib",rpath="/usr/local/lib" ,uselib_store="cc_omp",mandatory=1,uselib="cc_omp")
  except Exception as e:
    ctx.env["CCFLAGS_cc_omp"]=[]
    
  #ctx.env.append_value("CCFLAGS_cc_omp","-fopenmp")

def do_gcc(ctx):
  ctx.env.CC=[]
  ctx.env.LINK_CC=[]
  ctx.load('gcc')
  ctx.start_msg("Check gcc version") 
  version = [int(v) for v in ctx.env.CC_VERSION]
  if version[0]<4 or (version[0]==4 and version[1]<2):
    ctx.end_msg(".".join(ctx.env.CC_VERSION),color="YELLOW")
    raise Errors.WafError("gcc version need to be above 4.2 got %s"%".".join(ctx.env.CC_VERSION))
  ctx.end_msg(".".join(ctx.env.CC_VERSION))
  ctx.env.GCC_VERSION = ".".join(ctx.env.CC_VERSION)
  ctx.check_cc(
    errmsg="failed",msg="Compile a test code with gcc",
    mandatory=1,fragment = "#include <stdio.h>\nmain() {fprintf(stderr,\"hello world\");}\n",compile_filename='test.c',features='c cprogram')
  ctx.env["CCFLAGS_cc_omp"]=[]
  ctx.env.append_value("CCFLAGS_cc_omp","-fopenmp")

def configure_iccfirst(ctx):
  import os
  if "CC" in os.environ:
    ctx.start_msg("Using C compiler path from 'CC' environment variable")
    ctx.end_msg(os.environ["CC"])

  from waflib import Logs
  ctx.env.has_icc = False
  if not Options.options.gcc:
    try:
      do_icc(ctx)    
      return
    except:
      if Options.options.icc:
        raise
      Logs.pprint("PINK", "icc not found, defaulting to gcc")
  do_gcc(ctx)


def configure_gccfirst(ctx):
  from waflib import Logs
  from waflib import  Context
  from waflib import Errors

  ctx.env.has_icc = False
  import re  
  if not Options.options.icc and not Options.options.clang :
    try:
      do_gcc(ctx)  
      return
    except Exception as e:
      if Options.options.gcc:
        raise
      Logs.pprint("PINK", "gcc not found(cause : %s)"%e)
  if not Options.options.clang :
    try:
      do_icc(ctx)
      return 
    except Exception as e:
      if Options.options.icc:
        raise
      Logs.pprint("PINK", "icc not found (cause : %s)"%e)
  do_clang(ctx)

def configure(ctx):
  configure_gccfirst(ctx)
  ctx.env.append_value("CFLAGS_cshlib",ctx.env.LINKFLAGS_cshlib)  
  ctx.env["CFLAGS_cpic"]=[]
  ctx.env.append_value("CFLAGS_cpic",[flg for flg in ctx.env.CFLAGS_cshlib if "-fpic" in flg.lower()])
  
