from clik.parametric cimport c_parametric, error, doError, parametric, parametric_template

cdef extern c_parametric *ir_clustered_pep_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_clustered_in, error **err)
cdef extern c_parametric *ir_clustered_pep2_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_clustered_in, error **err)
cdef extern c_parametric *ir_poisson_pep_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_clustered_in, error **err)
cdef extern c_parametric *ir_poisson_guilaine_init(int ndet, double *detlist, int ndef, char** defkey, char **defvalue, int nvar, char **varkey, int lmin, int lmax, double* rq_clustered_in, error **err)

cdef class ir_clustered_pep(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> ir_clustered_pep_init;
    self.template_name = "ir_clustered_pep.dat"
    self.plugin_name = "pep_cib"

cdef class ir_clustered_pep2(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> ir_clustered_pep2_init;
    self.template_name = "ir_clustered_pep2.dat"
    self.plugin_name = "pep_cib"

cdef class ir_poisson_pep(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> ir_poisson_pep_init;
    self.template_name = "ir_poisson_pep.dat"
    self.plugin_name = "pep_cib"

cdef class ir_poisson_guilaine(parametric_template):
  def __cinit__(self):
    self.initfunc = <void*> ir_poisson_guilaine_init;
    self.template_name = "ir_poisson_guilaine.dat"
    self.plugin_name = "pep_cib"

component_list = ["ir_poisson_pep","ir_poisson_guilaine","ir_clustered_pep","ir_clustered_pep2"]
