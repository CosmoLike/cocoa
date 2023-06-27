from clik.parametric cimport c_parametric, error, doError, parametric, parametric_template

cdef extern c_parametric *radiogal_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *ir_poisson_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *ir_clustered_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)


cdef class radiogal(parametric):
  def __cinit__(self):
    self.initfunc = <void*>radiogal_init


cdef class ir_poisson(parametric):
  def __cinit__(self):
    self.initfunc = <void*>ir_poisson_init

cdef class ir_clustered(parametric):
  def __cinit__(self):
    self.initfunc = <void*>ir_clustered_init

    
component_list = ["radiogal","ir_clustered","ir_poisson"]

