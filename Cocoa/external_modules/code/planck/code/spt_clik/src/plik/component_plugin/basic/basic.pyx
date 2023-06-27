from clik.parametric cimport c_parametric, error, doError, parametric, parametric_template

cdef extern c_parametric *poisson_tensor_bydet_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *powerlaw_tensor_bydet_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *poisson_tensor_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *powerlaw_tensor_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *powerlaw_triangle_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *powerlaw_tanh_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)
cdef extern c_parametric *pointsource_init(int ndet, int *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, error **err)


cdef class poisson_tensor_bydet(parametric):
  def __cinit__(self):
    self.initfunc = <void*>poisson_tensor_bydet_init


cdef class powerlaw_tensor_bydet(parametric):
  def __cinit__(self):
    self.initfunc = <void*>powerlaw_tensor_bydet_init

cdef class poisson_tensor(parametric):
  def __cinit__(self):
    self.initfunc = <void*>poisson_tensor_init


cdef class powerlaw_tensor(parametric):
  def __cinit__(self):
    self.initfunc = <void*>powerlaw_tensor_init

cdef class powerlaw_triangle(parametric):
  def __cinit__(self):
    self.initfunc = <void*>powerlaw_triangle_init

cdef class powerlaw_tanh(parametric):
  def __cinit__(self):
    self.initfunc = <void*>powerlaw_tanh_init

cdef class pointsource(parametric):
  def __cinit__(self):
    self.initfunc = <void*>pointsource_init

    
component_list = ["pointsource","poisson_tensor","powerlaw_tensor","powerlaw_triangle","powerlaw_tanh","poisson_tensor_bydet","powerlaw_tensor_bydet"]

