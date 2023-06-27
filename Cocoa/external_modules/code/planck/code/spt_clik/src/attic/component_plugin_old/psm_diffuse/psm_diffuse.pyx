from clik.parametric cimport c_parametric, error, doError, parametric, parametric_template

cdef extern c_parametric *psm_diffuse_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)

cdef class psm_diffuse(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> psm_diffuse_init;
    self.template_name = "psm_diffuse.dat"
    self.plugin_name = "psm_diffuse"

component_list = ["psm_diffuse"]

