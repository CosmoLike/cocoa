from clik.parametric cimport c_parametric, error, doError, parametric, parametric_template, parametric_pol, parametric_pol_template,parametric_mul,parametric_template_mul
from clik.parametric import norename,rename_machine,rename_replace

cdef extern c_parametric *bleak_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)
cdef extern c_parametric *cnoise_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)
cdef extern c_parametric *cnoise_freepow_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *dip_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *feature_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *mul0_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *beamnl_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err)
cdef extern c_parametric *nslb_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err)
cdef extern c_parametric *trans_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err)
cdef extern c_parametric *trans2_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err)
cdef extern c_parametric *transG_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* template, error **err)

cdef class bleak(parametric_pol_template):
  def __cinit__(self):
    self.initfunc = <void*> bleak_init;
    self.template_name = "sky_template_M40404000_F100_143_217_353.dat"
    self.plugin_name = "systematics"

cdef class cnoise(parametric_pol_template):
  def __cinit__(self):
    self.initfunc = <void*> cnoise_init;
    self.template_name = "cnoise_F100_143_217_353.dat"
    self.plugin_name = "systematics"

cdef class cnoise_freepow(parametric_pol):
  def __cinit__(self):
    self.initfunc = <void*> cnoise_freepow_init;
    
cdef class dip(parametric):
  def __cinit__(self):
    self.initfunc = <void*> dip_init;

cdef class mul0(parametric_mul):
  def __cinit__(self):
    self.initfunc = <void*> mul0_init;
    
cdef class feature(parametric):
  def __cinit__(self):
    self.initfunc = <void*> feature_init;

cdef class beamnl(parametric_template_mul):
  def __cinit__(self):
    self.initfunc = <void*> beamnl_init;
    self.template_name = "beamnl_DX11_HM.dat"
    self.plugin_name = "systematics"

cdef class nslb(parametric_pol_template):
  def __cinit__(self):
    self.initfunc = <void*> nslb_init;
    self.template_name = "sky_and_beam_template_v15_F100_143_217_353.dat"
    self.plugin_name = "systematics"

cdef class trans(parametric_pol_template):
  def __cinit__(self):
    self.initfunc = <void*> trans_init;
    self.template_name = "sky_and_beam_template_v15_F100_143_217_353.dat"
    self.plugin_name = "systematics"

cdef class trans2(parametric_pol_template):
  def __cinit__(self):
    self.initfunc = <void*> trans2_init;
    self.template_name = "sky_and_beam_template_v15_F100_143_217_353.dat"
    self.plugin_name = "systematics"

cdef class transG(parametric_pol_template):
  def __cinit__(self):
    self.initfunc = <void*> transG_init;
    self.template_name = "sky_and_beam_template_v15_F100_143_217_353.dat"
    self.plugin_name = "systematics"

cnoise_EEBB = rename_machine(cnoise,{},rename_replace("cnoise","cnEEBB"),data_file="base_cnoise_EE_BB_rd12_rc3plik.txt")
cr_143E_HM = rename_machine(cnoise,{"A_cr_143E_143_143_EE":"1","A_cr_143E_100_143_EE":"1","A_cr_143E_143_217_EE":"1"},rename_replace("cnoise","cr_143E"),data_file="cr_hdr_HM_2017.dat")
cr_143E_OE = rename_machine(cnoise,{"A_cr_143E_143_143_EE":"1","A_cr_143E_100_143_EE":"1","A_cr_143E_143_217_EE":"1"},rename_replace("cnoise","cr_143E"),data_file="cr_hdr_OE_2017.dat")
cnoise_e2e = rename_machine(cnoise,{"A_cnoise_e2e_143_143_EE":"1","A_cnoise_e2e_100_100_EE":"1","A_cnoise_e2e_217_217_EE":"1"},rename_replace("cnoise","cnoise_e2e"),data_file="cnoise_e2e_v2.dat")


cnoise_gpe = rename_machine(cnoise,{},norename,data_file="cnoise_GPE_F100_143_217_353.dat")
cnoise_t2 = rename_machine(cnoise,{},norename,data_file="cnoise_F100_143_217_353_t2.dat")
cnoise_t3 = rename_machine(cnoise,{},norename,data_file="cnoise_F100_143_217_353_t3.dat")
cnoise_v17 = rename_machine(cnoise,{},norename,data_file="cnoise_F100_143_217_353_v17.dat")
bleak_v15 = rename_machine(bleak,{},norename,data_file="sky_template_v15_F100_143_217_353.dat")
cleak_defs = {}
for f1 in (100,143,217):
  for f2 in (100,143,217):
    for s1 in "TE":
      for s2 in "TE":
        cleak_defs["A_cleak_%d_%d_%s%s"%(f1,f2,s1,s2)] = "1"
cleak_v1 = rename_machine(cnoise,cleak_defs,rename_replace("cnoise","cleak"),data_file="cleak_eh_dx11_full_naive_v1.dat")
cleak_v2 = rename_machine(cnoise,cleak_defs,rename_replace("cnoise","cleak"),data_file="cleak_eh_dx11_full_wght_v2.dat")
cleak_rd12rc3_v1 = rename_machine(cnoise,cleak_defs,rename_replace("cnoise","cleak"),data_file="cleak_eh_rd12rc3_v1.dat")
cleak_rd12rc3_v1_oe = rename_machine(cnoise,cleak_defs,rename_replace("cnoise","cleak"),data_file="cleak_eh_rd12rc3_v1_oe.dat")

sbpx_defs = {}
for f1 in (100,143,217):
  for f2 in (100,143,217):
    for s1 in "TE":
      for s2 in "TE":
        sbpx_defs["A_sbpx_%d_%d_%s%s"%(f1,f2,s1,s2)] = "1"
subpix_rd12rc3_v4_hm = rename_machine(cnoise,sbpx_defs,rename_replace("cnoise","sbpx"),data_file="sbpx_tmpl_v4_hm.dat")

beamnl_DX11_defs = {"beamnl_nfreq_template":"3",
                    "beamnl_lmax_template":"3000",
                    "beamnl_freq_0":"100",
                    "beamnl_freq_1":"143",
                    "beamnl_freq_2":"217"}

beamnl_DX11 = rename_machine(beamnl,beamnl_DX11_defs,norename)

component_list = ["bleak","cnoise","dip","cnoise_gpe","cnoise_t2","cnoise_t3","cnoise_v17","bleak_v15","feature","cnoise_freepow",
                  "cleak_v1","cleak_v2","mul0","beamnl","beamnl_DX11","cleak_rd12rc3_v1","cnoise_EEBB","nslb","cleak_rd12rc3_v1_oe","trans","trans2","cr_143E_OE","cr_143E_HM","subpix_rd12rc3_v4_hm","cnoise_e2e","transG"]
 