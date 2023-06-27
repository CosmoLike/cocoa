from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import re
import sys
from waflib import Options
import os.path as osp
from waflib import Logs
from waflib import  Context
from waflib import Errors

### common tools
def show_linkline(ctx,flavour,wafname):
  ctx.start_msg("%s link line"%flavour)
  ctx.end_msg(" ".join(["-L%s"%vv for vv in ctx.env["LIBPATH_"+wafname]])+" "+" ".join(["-l%s"%vv for vv in ctx.env["LIB_"+wafname]]))

def retrieve_intel_linkline(ctx,flavour,wafname,execnm,flags,rl0=[],optionname=""):
  import waflib
  import os
  
  if optionname and not getattr(ctx.options,optionname,""):
    ctx.start_msg("retrieve %s link line"%flavour)
  
    try:
        # do a dry run to retrieve all of the command
      magic_cmd = "%s %s -dryrun -dynamiclib -shared-intel -no-cxxlib dummy.f90"%(execnm,flags)
      llgo,llge = ctx.cmd_and_log(magic_cmd, output=waflib.Context.BOTH)
      
      # now try to parse the hell out of it !
      # first, look for the location of the ld command  
      # assume it is spelled ld
      parts = re.split("^\s*ld\s*\\\\",llge,flags=re.M)
      if len(parts)!=2:
        # it is not called ld, or there is a specific path!!!
        # is it set by a -Qlocation command ?
        gotcha = re.findall("-Qlocation,ld,.+?\s",llge)
        if gotcha:
          # I am picking the last one...
          nloc = gotcha[-1].strip()
          # remove a posible trailling \
          if nloc[-1]=="\\":
            nloc[-1]=" "
          nloc = nloc.strip()
          # location should be after the last comma
          ldpath = nloc.split(",")[-1]
          ldname = re.escape(ldpath+"/ld")
          
          parts = re.split("^\s*%s\s*\\\\"%ldname,llge,flags=re.M)
          if len(parts)!=2:
            # ah ! That did not work ! I am getting a bit desperate here...
            # at this point I will grab anything that starts from the begining of the line with something different from a space and ends with  /ld some spaces and a trailing \ in it you know !
            parts = re.split("^\S.+?/ld\s+\\\\",llge,flags=re.M)
            if len(parts)!=2:
              # ok, I call it quit...
              raise
              
      # now parts[1] contains the link command ! Yeah !
      
      # grab the paths (assuming that they contains either intel or ifort in them)
      L = set([ll.strip() for ll in re.findall("^\s+-L(.+)\s*\\\\", parts[1],flags=re.M)  if ("ifort" in ll.lower()) or ("intel" in ll.lower())])

      # grab the libs
      l = set([ll.strip() for ll in re.findall("^\s+-l(.+)\s*\\\\", parts[1],flags=re.M)])
      
      # reduce to what's really relevant
      rL = set()
      rl = set(rl0)
      for Li in L:
        if osp.exists(Li):
          oli = os.listdir(Li)
          for li in l:
            if ctx.env.cshlib_PATTERN%li in oli:
              rl.add(li)
              rL.add(Li)
    except Exception as e:
      #print(e)
      ctx.end_msg(False)
      Logs.pprint("RED","Cannot retrieve the %s link line"%flavour)
      if optionname:
        Logs.pprint("RED","Please run manually")
        Logs.pprint("RED","  '%s'"%magic_cmd)
        Logs.pprint("RED","build the link line and pass it using the option")
        Logs.pprint("RED","  '%s'"%optionname)                
      ctx.fatal('The configuration failed') 
      raise
    for pth in list(rL) + ["/lib","/lib64"]:
      if osp.exists(pth):
        ctx.env.append_value("LIBPATH_%s"%wafname,pth)
        ctx.env.append_value("RPATH_%s"%wafname,pth)
    
    ctx.env.append_value("LIB_%s"%wafname,list(rl)+["pthread"])
    ctx.end_msg(True)
  show_linkline(ctx,flavour,wafname)


def options(ctx):
  import optparse
  grp = ctx.parser.get_option_group("--gcc")
  if grp==None:
    grp=optparse.OptionGroup(ctx.parser,"compiler options")
  grp.add_option("--gfortran",action="store_true",default=False,help="Do not test for ifort and only use gfortran")
  grp.add_option("--ifort",action="store_true",default=False,help="Do not test for gfortran and only use ifort")
  grp.add_option("--fortran_flagline",action="store",default="",help="flagline to link fortran object to c using ld")
  ctx.add_option_group(grp)  
  
def configure_(ctx):
  if ctx.options.fortran_flagline:
    conf.parse_flags(ctx.options.fortran_flagline,uselib="fc_runtime")
  if sys.platform.lower()=="darwin":
    ctx.env.fcshlib_PATTERN = 'lib%s.dylib'
  
  
  ctx.env.has_ifort = False
  if not Options.options.gfortran:
    try:
      ifort_conf(ctx)
      return
    except Exception as e:
      if Options.options.ifort:
        raise
      Logs.pprint("PINK", "ifort not found, defaulting to gfortran (cause: '%s')"%e)
  gfortran_conf(ctx)

def configure(ctx): 
  import os
  if "FC" in os.environ:
    ctx.start_msg("Using fortran compiler path from 'FC' environment variable")
    ctx.end_msg(os.environ["FC"])
  configure_(ctx) 
  ctx.env.append_value("FCFLAGS_fcshlib",ctx.env.LINKFLAGS_fcshlib)  
  ctx.env["FCFLAGS_fpic"]=[]
  ctx.env.append_value("FCFLAGS_fpic",[flg for flg in ctx.env.FCFLAGS_fcshlib if "-fpic" in flg.lower()])
  #ctx.env.append_value("FCFLAGS_fpic","-fpe0")

def ifort_conf(ctx):
  import waflib
  import os
  ctx.env.FC=[]
  ctx.load('ifort')
  if sys.platform.lower()=="darwin":
    ctx.env.LINKFLAGS_fcshlib = ['-dynamiclib']
  ctx.env.append_value('FCFLAGS',ctx.env.mopt)
  ctx.env["FCFLAGS_fc_omp"]=[]
  
  ctx.env.FCSHLIB_MARKER = [""]
  ctx.env.FCSTLIB_MARKER = [""]
  ctx.start_msg("Check ifort version") 
  v90 = ctx.cmd_and_log(" ".join(ctx.env.FC)+" --version",quiet=Context.STDOUT).split("\n")[0].strip()
  v90 = v90.split("\n")[0].strip().split(" ")[2]
  ctx.end_msg(v90)
  ctx.env.IFORT_VERSION = v90
  majver = int(v90.split(".")[0])
  rl0 = []
  if majver>13:
    rl0 = ["irc"]
  if majver>15:
    ctx.env.append_value("FCFLAGS_fc_omp","-qopenmp")  
  else:
    ctx.env.append_value("FCFLAGS_fc_omp","-openmp")
  
  ctx.check_cc(
    errmsg="failed",msg='Compile a test code with ifort',
    mandatory=1,fragment = "program test\n  WRITE(*,*) 'hello world'\n end program test\n",compile_filename='test.f90',features='fc fcprogram')

  retrieve_intel_linkline(ctx,"ifort","fc_runtime"," ".join(ctx.env.FC)," ".join(ctx.env.FCFLAGS+ctx.env.FCFLAGS_fc_omp),rl0,"fortran_flagline")
  #if not ctx.options.fortran_flagline:
  #  ctx.start_msg("retrieve ifort link line")
  #  try:
  #    #print "%s %s -dryrun -dynamiclib -shared-intel -no-cxxlib dummy.f90"%(ctx.env.FC," ".join(ctx.env.FCFLAGS))
  #    llgo,llge = ctx.cmd_and_log("%s %s -dryrun -dynamiclib -shared-intel -no-cxxlib dummy.f90"%(" ".join(ctx.env.FC)," ".join(ctx.env.FCFLAGS+ctx.env.FCFLAGS_fc_omp)), output=waflib.Context.BOTH)
  #    #print "RET",llgo,llge
  #    L = set([ll.strip() for ll in re.findall("^\s+-L(.+)\s*\\\\", re.split("^\s*ld\s*\\\\",llge,flags=re.M)[1],flags=re.M)  if ("ifort" in ll.lower()) or ("intel" in ll.lower())])
  #    l = set([ll.strip() for ll in re.findall("^\s+-l(.+)\s*\\\\", re.split("^\s*ld\s*\\\\",llge,flags=re.M)[1],flags=re.M)])
  #    rL = set()
  #    rl = set(rl0)
  #    for Li in L:
  #      if osp.exists(Li):
  #        oli = os.listdir(Li)
  #        for li in l:
  #          if ctx.env.cshlib_PATTERN%li in oli:
  #            rl.add(li)
  #            rL.add(Li)
  #  except:
  #    ctx.end_msg(False)
  #    raise
  #  for pth in list(rL) + ["/lib","/lib64"]:
  #    if osp.exists(pth):
  #      ctx.env.append_value("LIBPATH_fc_runtime",pth)
  #      ctx.env.append_value("RPATH_fc_runtime",pth)
  #  
  #  ctx.env.append_value("LIB_fc_runtime",list(rl)+["pthread"])
  #  ctx.end_msg(True)
  #show_linkline(ctx)
  ctx.env.has_ifort = True

  
def gfortran_conf(ctx):
  ctx.env.FC=[]
  ctx.env.FCFLAGS = []
  ctx.load('gfortran')
  ctx.env["FCFLAGS_fc_omp"]=[]
  ctx.env.append_value("FCFLAGS_fc_omp","-fopenmp")
  ctx.env.append_value("FCFLAGS","-DGFORTRAN")
  ctx.env.append_value("FCFLAGS","-ffixed-line-length-0")
  ctx.env.append_value("FCFLAGS","-ffree-line-length-0")
  mopt = ctx.env.mopt
  if sys.platform.lower()=="darwin":
    if "i386" in ctx.env.mopt:
      ctx.env.append_value('FCFLAGS','-m32')
      mopt = ["-m32"]
    else:
      ctx.env.append_value('FCFLAGS','-m64')
      mopt = ["-m64"]
  else:
    ctx.env.append_value('FCFLAGS',ctx.env.mopt)
  ctx.start_msg("Check gfortran version") 
  v90 = ctx.cmd_and_log(" ".join(ctx.env.FC)+" --version",quiet=Context.STDOUT).split("\n")[0].strip()
  version90 = re.findall("(4\.[0-9]\.[0-9])",v90)
  if len(version90)<1:
    #Logs.pprint("PINK","Can't get gfortran version... Let's hope for the best")
    ctx.end_msg("not found, let's hope for the best...",color="PINK")
  else:
    version90 = version90[0]
    vmid = int(version90.split(".")[1])
    if vmid<3:
      ctx.end_msg(v90,color="YELLOW")
      raise Errors.WafError("gfortran version need to be above 4.3 got %s"%version90)
    ctx.end_msg(v90)
  
  # kludge !
  ctx.env.FCSHLIB_MARKER = [""]
  ctx.env.FCSTLIB_MARKER = mopt
  ctx.check_cc(
      errmsg="failed",msg='Compile a test code with gfortran',
      mandatory=1,fragment = "program test\n  WRITE(*,*) 'hello world'\n end program test\n",compile_filename='test.f90',features='fc fcprogram')

  ctx.start_msg("retrieve gfortran link line")
  lgfpath = ctx.cmd_and_log(" ".join(ctx.env.FC)+" %s -print-file-name=libgfortran.dylib"%(" ".join(mopt)),quiet=Context.STDOUT)    
  lpath = [osp.dirname(osp.realpath(lgfpath))]
  lgfpath = ctx.cmd_and_log(" ".join(ctx.env.FC)+" %s -print-file-name=libgomp.dylib"%(" ".join(mopt)),quiet=Context.STDOUT)    
  lpath += [osp.dirname(osp.realpath(lgfpath))]
  lpath = set(lpath)

  ctx.env.append_value("LIB_fc_runtime",["gfortran","gomp"])
  ctx.env.append_value("LIBPATH_fc_runtime",list(lpath))
  ctx.env.append_value("RPATH_fc_runtime",list(lpath))
  ctx.end_msg(True)
  
  show_linkline(ctx,"gfortran","fc_runtime")




def NOTUSEDANYMORE_OLD_ifort_conf_(ctx):
  ctx.env.FC=[]
  ctx.load('ifort')
  if sys.platform.lower()=="darwin":
    ctx.env.LINKFLAGS_fcshlib = ['-dynamiclib']
  ctx.env.append_value('FCFLAGS',ctx.env.mopt)
  ctx.env.append_value("FCFLAGS_fc_omp","-openmp")
  ctx.env.FCSHLIB_MARKER = [""]
  ctx.env.FCSTLIB_MARKER = [""]
  
  ctx.check_cc(
    errmsg="failed",msg='Compile a test code with ifort',
    mandatory=1,fragment = "program test\n  WRITE(*,*) 'hello world'\n end program test\n",compile_filename='test.f90',features='fc fcprogram')
  if not ctx.options.fortran_flagline:
    ctx.start_msg("retrieve ifort link line")
    if "/" not in ctx.env.FC[0]:
      ctx.env.FC = ctx.cmd_and_log("which %s"%ctx.env.FC[0]).strip()
      #print ctx.env.FC
    ifort_path = osp.dirname(osp.realpath(ctx.env.FC[0]))
    
    #print ifort_path
    if ctx.options.m32:
      try:
        f=open(osp.join(ifort_path,'ifortvars_ia32.sh'))
      except:
        ctx.end_msg(False)
        raise Errors.WafError("Can't locate ifort configuration file")
    else:
      try:
        f=open(osp.join(ifort_path,'ifortvars_intel64.sh'))
      except:
        ctx.end_msg(False)
        raise Errors.WafError("Can't locate ifort configuration file")

    txt = f.read()
    f.close()
    #print txt
    if sys.platform.lower()=="darwin":
      sp = "DYLD_LIBRARY_PATH"
    else:
      sp = "LD_LIBRARY_PATH"
    res = re.findall("\s"+sp+"\s*=\s*\"(.+)\"",txt)[0]
    for pth in res.split(":"):
      ctx.env.append_value("LIBPATH_fc_runtime",pth)
      ctx.env.append_value("RPATH_fc_runtime",pth)
    ctx.env.append_value("LIB_fc_runtime",["ifcore","intlc","ifport","imf","irc","svml","iomp5","pthread"])
    ctx.end_msg(True)
  show_linkline(ctx)  
