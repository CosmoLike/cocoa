import autoinstall_lib as atl

version = "gsl-1.14"
tool = "gsl"

#print "-> loading %s autoinstall (using version %s)"%(tool,version)

def options(opt):
  atl.add_lib_option(tool,opt,install=True)
  
def configure(ctx):
  atl.conf_lib(ctx,tool,["gsl","gslcblas"],["gsl_ran_gaussian","gsl_blas_dgemm","gsl_linalg_cholesky_invert"],["gsl/gsl_randist.h","gsl/gsl_blas.h","gsl/gsl_linalg.h"],install=installgsl)
    
def installgsl(ctx):
  atl.installsmthg_pre(ctx,"ftp://ftp.gnu.org/gnu/gsl/"+version+".tar.gz",version+".tar.gz")
  atl.installsmthg_post(ctx,version,tool)
