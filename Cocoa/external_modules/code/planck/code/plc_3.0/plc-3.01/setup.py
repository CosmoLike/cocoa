from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy as nm
import os
import re

version = open("svnversion").read()+" MAKEFILE"
prefix = "./"
includes = [nm.get_include(),prefix+"/include","src/","src/python/clik","src/plik/"]
defines = [("HAS_LAPACK",None),("LAPACK_CLIK",None),("NOHEALPIX",None),("CLIK_LENSING",None)]
options = {}
options["include_dirs"] = includes
options["define_macros"] = defines
import os
link_clik = os.environ["LINK_CLIK"]
link_path = re.findall("-L(.+?)\s",link_clik)
link_cmd = re.findall("-l(.+?)\s",link_clik)

options["libraries"] = link_cmd
options["library_dirs"] = link_path
setup(name='clik',
      version=version,
      cmdclass = {'build_ext': build_ext},
      ext_modules = [Extension("clik.lkl", ["src/python/clik/lkl.pyx"],**options),
                     Extension("clik.lkl_lensing", ["src/python/clik/lkl_lensing.pyx"],**options),
                     Extension("clik.parametric", ["src/python/clik/parametric.pyx"],**options),
                     Extension("clik.rel2015", ["src/plik/component_plugin/rel2015/rel2015.pyx"],**options)
                     ],
      package_dir = {'clik': 'src/python/clik'},
      packages = ['clik'],
)