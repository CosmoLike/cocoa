#! PYTHONEXE
import sys
sys.path = ["REPLACEPATH"]+sys.path

import numpy as nm
import clik.parobject as php
import clik
import re
import os.path as osp
import os

class forfile:
  def __init__(self,fi):
    
    if isinstance(fi,file):
      self.fi = fi
    else :
      self.fi=open(fi)
    self.bf=""
  def read(self,fmt=''):
    if self.bf=='':
      sz = nm.fromstring(self.fi.read(4),dtype=nm.int32)[0]
      #print "want %d bytes"%sz
      self.bf = self.fi.read(sz)
      #print self.bf
      sz2 =nm.fromstring(self.fi.read(4),dtype=nm.int32)[0]
      #print sz2 
      assert sz==sz

    if fmt=='':
      self.bf=''
      return
    res = [self.cvrt(ff) for ff in fmt.strip().split()]
    if len(res)==1:
      return res[0]
    return tuple(res)
  
  def cvrt(self,fmt):
    cmd = re.findall("([0-9]*)([i|f])([0-9]+)",fmt)[0]
    dtype = nm.dtype({"f":"float","i":"int"}[cmd[1]]+cmd[2])
    itm = nm.array(1,dtype=dtype).itemsize
    nelem=1
    if cmd[0]: 
      nelem = int(cmd[0])
    res = nm.fromstring(self.bf[:itm*nelem],dtype=dtype)
    self.bf=self.bf[itm*nelem:]
    if nelem==1:
      return res[0]
    return res
  
  def close(self):
    self.bf=''
    self.fi.close()
  
    
def main(argv):
  pars = clik.miniparse(argv[1])
  if pars.int(default=1).camspec_version == 1:
    return main_v1(argv)
  return main_v2(argv)

def main_v1(argv):
  pars = clik.miniparse(argv[1])
  fcs = forfile(pars.input_likefile)

  Nspec,nX = fcs.read("i32 i32")

  lminX = []
  lmaxX = []
  np = []
  npt = []
  for i in range(Nspec):
    r = fcs.read("i32 "*4)
    lminX += [r[0]]
    lmaxX += [r[1]]
    np += [r[2]]
    npt += [r[3]]

  X = fcs.read("%df64"%nX)
  fcs.read()
  c_inv = fcs.read("%df64"%(nX)**2)
  

  tsz = nm.loadtxt(pars.input_tszfile)
  ksz = nm.loadtxt(pars.input_kszfile)
  tszXcib = nm.loadtxt(pars.input_tszXcibfile)

  tsz_temp = nm.zeros(5001)
  ksz_temp = nm.zeros(5001)
  tszXcib_temp = nm.zeros(5001)

  tsz_temp[tsz[0,0]:]=tsz[:5001-tsz[0,0],1]
  ksz_temp[ksz[0,0]:]=ksz[:5001-ksz[0,0],1]
  tszXcib_temp[tszXcib[0,0]:]=tszXcib[:5001-tszXcib[0,0],1]

  beam_Nspec = 0
  if "input_beamfile" in pars:
    fcs = forfile(pars.input_beamfile)
    beam_Nspec,num_modes_per_beam,beam_lmax = fcs.read("i32 "*3)
    cov_dim = beam_Nspec * num_modes_per_beam
    beam_modes = fcs.read("%df64"%(beam_Nspec*(beam_lmax+1)*num_modes_per_beam))
    fcs.read()
    beam_cov_inv = fcs.read("%df64"%(cov_dim*cov_dim))

  lmin = nm.min(lminX)
  lmax = nm.max(lmaxX)
  
  hascl = [0]*6
  hascl[0] = 1
  hascl = nm.array(hascl,dtype=nm.int)
  
  root_grp,hf = php.baseCreateParobject(pars.res_object)
  lkl_grp = php.add_lkl_generic(root_grp,"CAMspec",1,hascl,lmax,lmin)

  lkl_grp.attrs["Nspec"] = Nspec
  
  lkl_grp.attrs["lminX"] = lminX
  lkl_grp.attrs["lmaxX"] = lmaxX
  
  lkl_grp.attrs["np"] = np
  lkl_grp.attrs["npt"] = npt

  lkl_grp.attrs["nX"]          = nX
  
  lkl_grp.attrs["lmax_sz"]     = 5000

  lkl_grp.create_dataset('X', data=X)
  lkl_grp.create_dataset('c_inv', data=c_inv)
  lkl_grp.create_dataset('tsz', data=tsz_temp)
  lkl_grp.create_dataset('ksz', data=ksz_temp)
  lkl_grp.create_dataset('tszXcib', data=tszXcib_temp)

  lkl_grp.attrs["beam_Nspec"] = beam_Nspec
  if beam_Nspec!=0:
    lkl_grp.attrs["beam_lmax"] = beam_lmax
    lkl_grp.attrs["num_modes_per_beam"] = num_modes_per_beam
    lkl_grp.attrs["cov_dim"] = cov_dim
    lkl_grp.create_dataset('beam_cov_inv', data=beam_cov_inv)
    # beware beam_modes wants to be transposed...
    beam_modes.shape = (num_modes_per_beam,beam_lmax+1,beam_Nspec)
    bm = beam_modes.T*1.
    lkl_grp.create_dataset('beam_modes', data=bm.flat[:])

  php.add_pid(lkl_grp,pars.str(default="").pid)

  lkl_grp.attrs["has_dust"] = pars.int(default=0).has_dust
  lkl_grp.attrs["has_calib_prior"] = pars.int(default=1).has_calib_prior

  if "beam_mode" in pars:
    beam_mode = nm.loadtxt(pars.str.beam_mode)
    beam_flag = [int(v in ("1","t")) for v in pars.str.free_beam_flag.ljust(cov_dim,"0")]
    marge_flag = [1-v for v in beam_flag]
    print(marge_flag)
    lkl_grp.attrs["marge_mode"] = beam_mode.flat[:]
    lkl_grp.attrs["marge_flag"] = nm.array(marge_flag).flat[:]

  if "bs_factor" in pars:
    lkl_grp.attrs["bs_factor"] = pars.float.bs_factor

  hf.close()
  
def main_v2(argv):
  pars = clik.miniparse(argv[1])

  version = pars.int(default=2).camspec_version

  root_grp,hf = php.baseCreateParobject(pars.res_object)
  hascl = nm.array(pars.bool_array.hascl,dtype=nm.int)
  lmin = pars.int.lmin
  lmax = pars.int.lmax
  print(hascl,hascl.dtype)
  lkl_grp = php.add_lkl_generic(root_grp,"CAMspec",1,hascl,lmax,lmin)

  lkl_grp.attrs["camspec_version"] = version

  if "pre_marged" in pars:
    print(pars.int.pre_marged)
    lkl_grp.attrs["pre_marged"] = pars.int.pre_marged

  if "spec_flag" in pars:
    spec_flag = pars.bool_array.spec_flag
    assert len(spec_flag)==6,"argh len spec_flag !=6"
    lkl_grp.attrs["spec_flag"] = [int(v) for v in spec_flag]

    lmins = pars.int_array.lmins
    assert len(lmins)==6,"argh len lmins !=6"
    lkl_grp.attrs["spec_lmins"] = lmins
    lmaxs = pars.int_array.lmaxs
    assert len(lmaxs)==6,"argh len lmaxs !=6"
    lkl_grp.attrs["spec_lmaxs"] = lmaxs

  lkl_grp.attrs["camspec_beam_mcmc_num"] = pars.int(default=1).camspec_beam_mcmc_num

  lkl_grp.attrs["bs_factor"] = pars.float(default=1).bs_factor

  lkl_grp.attrs["sz_prior"] = pars.int(default=1).sz_prior
  
  nuisance_pars_v2 = ["aps100",  "aps143",  "aps217",  "acib143",  "acib217",  
                   "asz143",  "psr",  "cibr",  "ncib143",  "ncib",  "cibrun",  
                   "xi",  "aksz",  "wig1_143",  "wig1_217",  "wig1_r",  
                   "wig1_L",  "wig1_sigma",  "wig2_143",  "wig2_217",  "wig2_r",
                   "wig2_L",  "wig2_sigma",  "cal0",  "cal1",  "cal2",  "calTE",
                   "calEE",  "bm_1_1",  "bm_1_2",  "bm_1_3",  "bm_1_4",  
                   "bm_1_5",  "bm_2_1",  "bm_2_2",  "bm_2_3",  "bm_2_4",  
                   "bm_2_5",  "bm_3_1",  "bm_3_2",  "bm_3_3",  "bm_3_4",  
                   "bm_3_5",  "bm_4_1",  "bm_4_2",  "bm_4_3",  "bm_4_4",  
                   "bm_4_5"]
  nuisance_pars_v3 = ["aps100",  "aps143",  "aps217",  "acib143",  "acib217",  
                      "asz143",  "psr",  "cibr",  "ncib143",  "ncib",  "cibrun",  
                      "xi",  "aksz",  "dust100",  "dust143",  "dust217",  
                      "dust143x217",  "wig1_143",  "wig1_217",  "wig1_r",  
                      "wig1_L",  "wig1_sigma",  "wig2_143",  "wig2_217",  "wig2_r",  
                      "wig2_L",  "wig2_sigma",  "calPlanck",  "cal0",  "cal1",  
                      "cal2",  "calTE",  "calEE",  "bm_1_1",  "bm_1_2",  "bm_1_3",
                        "bm_1_4",  "bm_1_5",  "bm_2_1",  "bm_2_2",  "bm_2_3",  
                        "bm_2_4",  "bm_2_5",  "bm_3_1",  "bm_3_2",  "bm_3_3",  
                        "bm_3_4",  "bm_3_5",  "bm_4_1",  "bm_4_2",  "bm_4_3",  
                        "bm_4_4",  "bm_4_5"]

  if version == 2:
    nuisance_pars = nuisance_pars_v2
  else:
    nuisance_pars = nuisance_pars_v3

  lkl_grp.attrs["n_nuisance"] = len(nuisance_pars)
  lkl_grp.attrs["nuisance"] = php.pack256(*nuisance_pars)

  import tempfile
  #dr = tempfile.mkdtemp()
  dr = "tmp_camspec"
  os.mkdir(dr)
  print(dr)
  import shutil
  #like_file[20], sz143_file[20], tszxcib_file[20], ksz_file[20], beam_file[20],data_vector[20]

  shutil.copy(pars.like_file.strip(),dr+"/like_file")
  shutil.copy(pars.sz143_file.strip(),dr+"/sz143_file")
  shutil.copy(pars.tszxcib_file.strip(),dr+"/tszxcib_file")
  shutil.copy(pars.ksz_file.strip(),dr+"/ksz_file")
  shutil.copy(pars.beam_file.strip(),dr+"/beam_file")
  if "data_vector" in pars:
    lkl_grp.attrs["data_vector_flag"] = 1
    shutil.copy(pars.data_vector.strip(),dr+"/data_vector")
  if "camspec_fiducial_cl" in pars:
    ext=""
    lkl_grp.attrs["minimum.theory_cl"] = 0
    if "minimum.theory_cl" in pars.camspec_fiducial_cl.strip():
      ext = ".minimum.theory_cl"
      lkl_grp.attrs["minimum.theory_cl"] = 1
    shutil.copy(pars.camspec_fiducial_cl.strip(),dr+"/camspec_fiducial_cl"+ext)
    shutil.copy(pars.camspec_fiducial_foregrounds.strip(),dr+"/camspec_fiducial_foregrounds")
  if version == 3:
    if "cib217_file" in pars:  
      shutil.copy(pars.cib217_file.strip(),dr+"/cib_file")
      lkl_grp.attrs["cib_consistency_flag"] = 1
    if "dust217_file" in pars:
      lkl_grp.attrs["dust_flag"] = 1
      shutil.copy(pars.dust217_file.strip(),dr+"/dust217_file")
      shutil.copy(pars.dust143_file.strip(),dr+"/dust143_file")
      shutil.copy(pars.dust143x217_file.strip(),dr+"/dust143x217_file")
      shutil.copy(pars.dust100_file.strip(),dr+"/dust100_file")
      
  php.add_external_data(dr,lkl_grp,tar=True)

  if "default_name" in pars:
    name = pars.str_array.default_name
    loc = pars.float_array.default_val
    assert len(name)==len(loc),"name and loc have different sizes"
    php.add_default(hf["clik"],name,loc,nuisance_pars)
  

  hf.close()

  shutil.rmtree(dr)

import sys
if __name__=="__main__":
  main(sys.argv)
