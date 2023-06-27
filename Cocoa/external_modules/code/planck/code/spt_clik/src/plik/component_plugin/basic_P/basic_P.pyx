from clik.parametric cimport c_parametric, error, doError, parametric, parametric_template, parametric_pol, parametric_pol_template
from clik.parametric import powerlaw_free_emissivity,rename_machine

cdef extern c_parametric *pw_TE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *pw_EE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)

cdef extern c_parametric *powerlaw_free_emissivity_EE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *powerlaw_free_emissivity_TE_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *powerlaw_free_emissivity_XX_init(int ndet_T, int ndet_P, int *has_TEB, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)

cdef class pw_TE(parametric_pol):
  def __cinit__(self):
    self.initfunc = <void*>pw_TE_init

cdef class pw_EE(parametric_pol):
  def __cinit__(self):
    self.initfunc = <void*>pw_EE_init

cdef class powerlaw_free_emissivity_EE(parametric_pol):
  def __cinit__(self):
    self.initfunc = <void*>powerlaw_free_emissivity_EE_init

cdef class powerlaw_free_emissivity_TE(parametric_pol):
  def __cinit__(self):
    self.initfunc = <void*>powerlaw_free_emissivity_TE_init

cdef class powerlaw_free_emissivity_XX(parametric_pol):
  def __cinit__(self):
    self.initfunc = <void*>powerlaw_free_emissivity_XX_init

def galf_P_rename_func(v,rups):
  if v.startswith("galf"):
    rv = v.replace("galf","pwfe")
    rups[v]=rv

def galres_TT_rename_func(v,rups):
  if v.startswith("galres_TT"):
    rv = v.replace("galres_TT","pwfe")
    rups[v]=rv

def galres_TE_rename_func(v,rups):
  if v.startswith("galres_TE"):
    rv = v.replace("galres_TE","pwfe_XX")
    rups[v]=rv

def dalf_P_rename_func(v,rups):
  if v.startswith("dalf"):
    rv = v.replace("dalf","pwfe")
    rups[v]=rv

def mk_galx_P_rename_func(nn):
  def galx_P_rename_func(v,rups):
    if v.startswith("galf_%s"%nn):
      rv = v.replace("galf_%s"%nn,"pwfe_XX")
      rups[v]=rv
  return galx_P_rename_func


galf_TE = rename_machine(powerlaw_free_emissivity_TE,{"galf_TE_l2_norm":"1","galf_TE_l_pivot":"500","galf_TE_index":"-2.4"},galf_P_rename_func)
galf_EE = rename_machine(powerlaw_free_emissivity_EE,{"galf_EE_l2_norm":"1","galf_EE_l_pivot":"500","galf_EE_index":"-2.4"},galf_P_rename_func)
dalf_EE = rename_machine(powerlaw_free_emissivity_EE,{"dalf_EE_l2_norm":"1","dalf_EE_l_pivot":"500","dalf_EE_index":"-2."},dalf_P_rename_func)

galx_EE = rename_machine(powerlaw_free_emissivity_XX,{"galf_EE_l2_norm":"1","galf_EE_l_pivot":"500","galf_EE_index":"-2.4","galf_EE_kind":"1"},mk_galx_P_rename_func("EE"))
galx_BB = rename_machine(powerlaw_free_emissivity_XX,{"galf_BB_l2_norm":"1","galf_BB_l_pivot":"500","galf_BB_index":"-2.4","galf_BB_kind":"2"},mk_galx_P_rename_func("BB"))
galx_TE = rename_machine(powerlaw_free_emissivity_XX,{"galf_TE_l2_norm":"1","galf_TE_l_pivot":"500","galf_TE_index":"-2.4","galf_TE_kind":"3"},mk_galx_P_rename_func("TE"))
galx_TB = rename_machine(powerlaw_free_emissivity_XX,{"galf_TB_l2_norm":"1","galf_TB_l_pivot":"500","galf_TB_index":"-2.4","galf_TB_kind":"4"},mk_galx_P_rename_func("TB"))
galx_EB = rename_machine(powerlaw_free_emissivity_XX,{"galf_EB_l2_norm":"1","galf_EB_l_pivot":"500","galf_EB_index":"-2.4","galf_EB_kind":"5"},mk_galx_P_rename_func("EB"))

galres_TT = rename_machine(powerlaw_free_emissivity,{"galres_TT_l2_norm":"1","galres_TT_l_pivot":"150","galres_TT_index":"-2.5"},galres_TT_rename_func)
galres_TE = rename_machine(powerlaw_free_emissivity_XX,{"galres_TE_l2_norm":"1","galres_TE_l_pivot":"150","galres_TE_index":"-3.1","galres_TE_kind":"3"},galres_TE_rename_func)

def ps_P_rename_func(v,rups):
  if v.startswith("ps"):
    rv = v.replace("ps","pwfe")
    rups[v]=rv

ps_TE = rename_machine(powerlaw_free_emissivity_TE,{"ps_TE_l2_norm":"1","ps_TE_l_pivot":"3000","ps_TE_index":"0"},ps_P_rename_func)
ps_EE = rename_machine(powerlaw_free_emissivity_EE,{"ps_EE_l2_norm":"1","ps_EE_l_pivot":"3000","ps_EE_index":"0"},ps_P_rename_func)

component_list = ["galres_TT","galres_TE","powerlaw_free_emissivity_XX","pw_TE","pw_EE","powerlaw_free_emissivity_EE","powerlaw_free_emissivity_TE","galf_TE","galf_EE","ps_TE","ps_EE","dalf_EE","galx_EE","galx_BB","galx_TE","galx_TB","galx_EB"]
