from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os.path as osp

import numpy as nm
import os

# remplace la chaine de charactere par le resultat de 
# gfortran -print-file-name=libgfortran.so 
# (ou sur mac
# gfortran -print-file-name=libgfortran.dylib )

fruntime = "/usr/local/gfortran/lib/gcc/i386-apple-darwin8.10.1/4.5.0/../../../libgfortran.dylib"
fruntime = "/usr/lib/gcc/x86_64-linux-gnu/4.4.5/libgfortran.so"

Lfruntime = osp.split(fruntime)[0]
Lfruntime = osp.realpath(Lfruntime)
lfruntime = "gfortran"

setup(
    name='egfs',
    version='1.0',
    packages=['egfs'],
    package_dir={'egfs': '.'},
    package_data={'egfs': ['data/*.dat']},
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("egfs", ["egfs.pyx","clik_egfs.c","distribution.c","io.c","errorlist.c", "mvdens.c"],
                             include_dirs = [nm.get_include()],
                             library_dirs=['.',Lfruntime],
                             libraries=["egfs",lfruntime, "gsl", "gslcblas" ]
                             ),
    			   Extension("parametric", ["parametric.pyx","clik_parametric.c","io.c","errorlist.c"],
                             include_dirs = [nm.get_include()],
                             library_dirs=['.',Lfruntime]
                             ),
    			   Extension("parametric", ["parametric.pyx","clik_parametric.c","io.c","errorlist.c"],
                             include_dirs = [nm.get_include()],
                             library_dirs=['.',Lfruntime]
                             )],
)
