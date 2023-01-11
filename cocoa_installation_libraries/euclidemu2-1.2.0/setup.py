from setuptools import setup, Extension, find_packages
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import os, platform
from distutils.sysconfig import get_python_lib
from site import getusersitepackages

# Setting the compiler to g++
os.environ["CC"] = "x86_64-conda-linux-gnu-g++"
os.environ["CXX"] = "x86_64-conda-linux-gnu-g++"

if 'LDFLAGS' in os.environ.keys():
    ldfl=os.environ['LDFLAGS']
    new_ldfl=ldfl.replace('-Wl,-dead_strip_dylibs ','')
    os.environ['LDFLAGS']=new_ldfl

# Getting the two possible locations for the installed files
pathtopythonlib1=get_python_lib()
pathtopythonlib2=getusersitepackages()

with open("README.md", 'r') as f:
    long_description = f.read()

pack_name='euclidemu2'

extensions=Extension(name=pack_name,
                           sources=["src/euclidemu2.pyx","src/cosmo.cxx","src/emulator.cxx"],
                           include_dirs=["../src/"],
                           libraries=["gsl","gslcblas"],
                           language="c++",
                           extra_compile_args=['-std=c++11',
                                               '-D PRINT_FLAG=0',
                                               '-D PATH_TO_EE2_DATA_FILE1="'+pathtopythonlib1+'/'+pack_name+'/ee2_bindata.dat"',
                                               '-D PATH_TO_EE2_DATA_FILE2="'+pathtopythonlib2+'/euclidemu2/ee2_bindata.dat"']
                           )


setup(name=pack_name,
      version="1.2.0",
      author="Pedro Carrilho,  Mischa Knabenhans",
      author_email="pedromgcarrilho@gmail.com",
      description="Python wrapper for EuclidEmulator2",
      long_description=long_description,
      long_description_content_type='text/markdown',
      url="https://github.com/PedroCarrilho/EuclidEmulator2/tree/pywrapper",
      cmdclass={'build_ext': build_ext},
      ext_modules = cythonize(extensions,language_level = 3),
      packages=[pack_name],
      package_dir={pack_name: 'src'},
      package_data={pack_name: ["ee2_bindata.dat","cosmo.h","emulator.h","units_and_constants.h"]},
      include_package_data=True,
      )
