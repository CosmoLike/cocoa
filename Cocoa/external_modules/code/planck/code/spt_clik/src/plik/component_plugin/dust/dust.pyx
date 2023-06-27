from clik.parametric cimport c_parametric, error, doError, parametric, parametric_template, parametric_pol, parametric_pol_template
from clik.parametric import powerlaw_free_emissivity,rename_machine,rename_replace
from clik.systematics import cnoise

cdef extern c_parametric *galactic_component_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *hgal_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *kgal_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *gpe_dust_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *gpegal_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *gal545_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *gal545_80pc_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern double c_dust_spectrum "dust_spectrum" (double nu, double T_dust, double beta_dust, double nu0)
cdef extern double c_non_thermal_spectrum "non_thermal_spectrum" (double nu, double alpha_non_thermal, double nu0)
cdef extern double c_dBdT "dBdT" (double nu, double nu0)
cdef extern c_parametric *gal_TE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *gal_EE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *t1gal_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)


def dust_spectrum(nu,T_dust=18.0,beta_dust=1.8,nu0=143.0):
  return c_dust_spectrum(<double>nu,<double>T_dust,<double>beta_dust,<double>nu0)

def non_thermal_spectrum(nu,alpha_non_thermal=-1.0,nu0=143.0):
  return c_non_thermal_spectrum(<double>nu,<double>alpha_non_thermal,<double>nu0)

def dBdT(nu,nu0=143.0):
  return c_dBdT(<double>nu,<double>nu0)


cdef class galametric(parametric):
  def __cinit__(self):
    self.initfunc = <void*>galactic_component_init

cdef class gal_TE(parametric_pol):
  def __cinit__(self):
    self.initfunc = <void*>gal_TE_init

cdef class gal_EE(parametric_pol):
  def __cinit__(self):
    self.initfunc = <void*>gal_EE_init

cdef class gpe_dust(parametric):
  def __cinit__(self):
    self.initfunc = <void*> gpe_dust_init;

cdef class hgal(parametric):
  def __cinit__(self):
    self.initfunc = <void*> hgal_init;

cdef class kgal(parametric):
  def __cinit__(self):
    self.initfunc = <void*> kgal_init;

cdef class gpegal(parametric):
  def __cinit__(self):
    self.initfunc = <void*> gpegal_init;

cdef class gal545(parametric):
  def __cinit__(self):
    self.initfunc = <void*> gal545_init;

cdef class gal545_80pc(parametric):
  def __cinit__(self):
    self.initfunc = <void*> gal545_80pc_init;

cdef class t1gal(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> t1gal_init;
    self.template_name = "y1y2_temp2_50pc.dat[1]"
    self.plugin_name = "dust"

def galf_rename_func(v,rups):
  if v.startswith("galf"):
    rv = v.replace("galf","pwfe")
    rvl = rv.split("_")
    if rvl[-1].isdigit() and not rvl[-2].isdigit():
      rvl+=[rvl[-1]]
    rups[v]="_".join(rvl)

galf = rename_machine(powerlaw_free_emissivity,{"galf_A_143":"0","galf_A_100":"0","galf_A_100_143":"0","galf_A_100_217":"0","galf_A_100_353":"0","galf_A_353":"0","galf_A_143_353":"0","galf_A_217_353":"0"},galf_rename_func)

def ngpegal_rename_func(v,rups):
  if v.startswith("t1gal"):
    rv = v.replace("t1gal","gal545")
    rvl = rv.split("_")
    if rvl[-1].isdigit() and not rvl[-2].isdigit():
      rvl+=[rvl[-1]]
    rups[v]="_".join(rvl)

ngpegal = rename_machine(t1gal,{"gal545_l_pivot":"200"},rename_replace("t1gal","gal545"),data_file="gpe_dust.dat")

gal545_2017_defs = {"gal545_abs" : "1",
                    "gal545_l_pivot" : "200",
                    "gal545_A_100" : "7",
                    "gal545_A_143" : "9",
                    "gal545_A_143_217" : "20",
                    "gal545_A_217" : "83",
                  }

gal545_2017_rename = {"gal545_A_217" : "A_cnoise_217_217_TT",
                      "gal545_A_100" : "A_cnoise_100_100_TT",
                      "gal545_A_143" : "A_cnoise_143_143_TT",
                      "gal545_A_143_217" : "A_cnoise_143_217_TT",
                      "gal545_abs" : "cnoise_abs",
                      "gal545_l_pivot" : "cnoise_l_pivot"}
              
import os
import os.path as osp
gal545_2017 = rename_machine(cnoise,gal545_2017_defs,gal545_2017_rename,data_path=osp.join(os.environ.get("CLIK_DATA"),"dust","dust_2017_SG_v1.dat"))
  

component_list = ["galametric","gpe_dust","gal_EE","gal_TE","galf","hgal","kgal","t1gal","gpegal","gal545","gal545_80pc","ngpegal","gal545_2017"]

