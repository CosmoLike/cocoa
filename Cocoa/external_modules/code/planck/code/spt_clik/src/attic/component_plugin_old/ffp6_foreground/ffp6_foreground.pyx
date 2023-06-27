from clik.parametric cimport c_parametric, error, doError, parametric, parametric_template

cdef extern c_parametric *ffp6_thermaldust_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)
cdef extern c_parametric *ffp6_firb_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_in, error **err)

cdef class ffp6_thermaldust(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> ffp6_thermaldust_init;
    self.template_name = "ffp6_thermaldust.dat"
    self.plugin_name = "ffp6_foreground"

cdef class ffp6_firb(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> ffp6_firb_init;
    self.template_name = "ffp6_firb.dat"
    self.plugin_name = "ffp6_foreground"


component_list = ["ffp6_thermaldust", "ffp6_firb"]

