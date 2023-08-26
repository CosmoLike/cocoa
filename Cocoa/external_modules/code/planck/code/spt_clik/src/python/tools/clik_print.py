#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import numpy as nm
import clik
import clik.hpy as hpy

def main(argv):
  if len(sys.argv)!=2:
    print("clik version "+clik.version())
    print("usage : %s lkl_file"%(sys.argv[0]))
    sys.exit(1)

  if clik.try_lensing(argv[1]):
    main_lensing(argv)
    return
  main_CMB(argv)

def main_CMB(argv):
  clikl = clik.clik(sys.argv[1])

  extn = clikl.extra_parameter_names
  
  lkl = hpy.File(sys.argv[1],"r")["clik"]
  
  print("clik lkl file =  %s"%sys.argv[1])
  print("  number of likelihoods = %d"%lkl.attrs["n_lkl_object"])
  print("  lmax ( "+ " ".join([nl+" = %d"%ll for nl,ll in zip(("TT","EE","BB","TE","TB","EB"),lkl.attrs["lmax"]) if ll >-1])+" )")
  print("  number of varying extra parameters %d"%(len(extn)))
  for n in extn:
    print("    %s"%n)
  if "prior" in lkl:
    names = lkl["prior"].attrs["name"]
    names = [names[i*256:(i+1)*256].strip() for i in range(len(names)/256)]
    print("  gaussian priors on %s"%", ".join(names))
    loc = lkl["prior/loc"][:]
    print("  at \n    %s"%" ".join([str(l) for l in loc]))
    var = lkl["prior/var"][:]
    if len(var)==len(loc):
      var = nm.diag(var)
    var.shape=(len(loc),-1)
    print("  with variance")
    print("\n".join(["    "+" ".join([str(v) for v in vl]) for vl in var]))
  if "default" in lkl:
    loc = lkl["default/loc"][:] 
    print("  number of fixed parameters = %d"%len(loc))
    nms = lkl["default"].attrs["name"]
    nms = [nms[i*256:i*256+256].strip() for i in range(len(loc))]
    for n,l in zip(nms,loc):
      print("    %s = %g"%(n,l))
  options = clikl.get_options()
  if len(options)!=0:
    print("  number of initialization options %d"%(len(options)))
    for n in options:
      print("    %s"%n)
  
  ilkl = 0
  for lkli_n in ("lkl_%d"%v for v in range(lkl.attrs["n_lkl_object"])):
    lkli = lkl[lkli_n]
    print("\n  %s"%lkli_n)
    print("    lkl_type = %s"%lkli.attrs["lkl_type"])
    print("    unit = %g"%lkli.attrs["unit"])
    
    if "lmax" in lkli.attrs:
      lmax = lkli.attrs["lmax"]
      lmin = 0
      if "lmin" in lkli.attrs:
        lmin = lkli.attrs["lmin"]
      ellh = False
    else:
      ell = lkli.attrs["ell"]
      lmax = nm.max(ell)      
      lmin = nm.min(ell)
      ellh = not nm.alltrue((ell[1:]-ell[:-1]) == 1)
    
    print("    "+" ".join([nl+" = [%d , %d]"%(lmin,lmax) for nl,hl in zip(("TT","EE","BB","TE","TB","EB"),lkli.attrs["has_cl"]) if hl ])+" (discontinous)"*ellh)
    
    if "wl" in lkli.attrs:
      print("    has window function")
    if "nbins" in lkli.attrs:
      print("    nbins = %d"%lkli.attrs["nbins"])
    if lkli.attrs["lkl_type"]=="smica":
      print("    component 0 : CMB")
      for nc in range(1,lkli.attrs["n_component"]):
        if "component_name" in lkli["component_%d"%nc]:
          print("    component %d : %s"%(nc,lkli["component_%d"%nc].attrs["component_name"]))
        else:
          print("    component %d : %s"%(nc,lkli["component_%d"%nc].attrs["component_type"]))

    extn = clikl.get_extra_parameter_names_by_lkl(ilkl)
    print("    number of extra parameters = %d %s"%(len(extn),extn))
    ilkl +=1

def main_lensing(argv):
  lkl = clik.clik_lensing(sys.argv[1])
  print("clik lensing file = %s"%sys.argv[1])
  if hpy.cldf.is_cldf(argv[1]):
    ff = hpy.File(argv[1])
    ty = ff["clik_lensing/itype"]
    if ty==0:
      renorm = 1
      ren1 = 0
    else:
      renorm = ff["clik_lensing/renorm"]
      ren1 = 0
      if ty in [3,4]:
        ren1 = ff["clik_lensing/ren1"]
  else:
    rr = file(argv[1]).read()
    if ("# format: mono") in rr:
      renorm = 1
      ren1 = 0
    elif ("# format: qecl") in rr:
      renorm = 1
      ren1 = 0
    else:
      renorm = 1
      ren1 = 1
    
  print("  Renormalization: %s\n  N1 computation: %s"%(["Off","On"][renorm],["Off","On"][ren1]))
  print("  lmax = %s (PP TT EE BB TE TB EB)"%lkl.lmax)
  print("  number of extra parameters %d"%len(lkl.extra_parameter_names))
  for nn in lkl.extra_parameter_names:
    print("     %s"%nn)

if __name__=="__main__":
  main(sys.argv)