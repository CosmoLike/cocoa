cdef extern from "errorlist.h":
  ctypedef struct error:
    pass
    
cdef object doError(error **err)

cdef extern from "clik_parametric.h":
  ctypedef struct c_parametric "parametric":
    int lmin,lmax,ndet,nfreq,nvar

ctypedef c_parametric* (*simple_init)(int , double *, int , char** , char **, int , char **, int , int , error **)
ctypedef c_parametric* (*template_init)(int , double *, int , char** , char **, int , char **, int , int , double*, error **)
ctypedef c_parametric* (*pol_init)(int , int , int*, double *, int , char** , char **, int , char **, int , int , error **)
ctypedef c_parametric* (*pol_template_init)(int , int , int*, double *, int , char** , char **, int , char **, int , int , double*, error **)

cdef class parametric:
  cdef c_parametric* celf
  cdef int nell
  cdef readonly object varpar,parvalues,defaults
  cdef void* initfunc
  cdef readonly object rename,emaner,voidmask,name
  cdef readonly object _template
  cdef readonly object ismul
  cdef object component_name
  
cdef class parametric_template(parametric):
  cdef object template_name
  cdef object plugin_name

cdef class parametric_pol(parametric):
  pass
cdef class parametric_pol_template(parametric_pol):
  cdef object template_name
  cdef object plugin_name

cdef class parametric_mul(parametric):
  pass

cdef class parametric_template_mul(parametric_template):
  pass
