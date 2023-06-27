from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#try to support many flavours of lapack
import autoinstall_lib as atl
from waflib import Logs
import os.path as osp
    
def options(ctx):
  atl.add_lib_option("pmc",ctx,install=False)
  
def configure(ctx):
  ctx.env.has_pmc = False
  #pmc_config_path = ctx.find_program("pmc-config",path_list=[ctx.options.pmc_prefix+"/bin"])[0]
  try:
    pmc_config_path = ctx.find_program("pmc-config",path_list=[ctx.options.pmc_prefix+"/bin"])
    pmcflagline = ctx.cmd_and_log(pmc_config_path)
  except Exception as e:
    pmcflagline='' 
  atl.conf_lib(ctx,"pmc","pmc","init_distribution","pmc.h",["pmclib","pmctools"],defines=["HAS_PMC"],flagline=pmcflagline)
    