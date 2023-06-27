from clik.parametric cimport c_parametric, error, doError, parametric, parametric_template

cdef extern c_parametric *psm_compact_bydet_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)

cdef class psm_compact_bydet(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> psm_compact_bydet_init;
    self.template_name = "psm_compact_bydet.dat"
    self.plugin_name = "psm_compact_bydet"

component_list = ["psm_compact_bydet"]

