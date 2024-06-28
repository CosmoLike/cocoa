#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

# --------------------------------------------------------------------------
# IF SET, THEN COCOA ADOPTS FFTW10. OTHERWISE, COCOA ADOPTS FFTW8
# --------------------------------------------------------------------------
#export FFTW_NEW_VERSION=1

# --------------------------------------------------------------------------
# IF SET, COCOA DOES NOT USE SYSTEM PIP PACKAGES 
# --------------------------------------------------------------------------
export DONT_USE_SYSTEM_PIP_PACKAGES=1

# --------------------------------------------------------------------------
# WE USE COLASLIM ENV WITH JUST PYTHON AND GCC TO TEST MANUAL INSTALLATION
# --------------------------------------------------------------------------
#conda create --name cocoalitepy38 python=3.8 --quiet --yes \
#   && conda install -n cocoalitepy38 --quiet --yes  \
#   'conda-forge::libgcc-ng=12.3.0' \
#   'conda-forge::libstdcxx-ng=12.3.0' \
#   'conda-forge::libgfortran-ng=12.3.0' \
#   'conda-forge::gxx_linux-64=12.3.0' \
#   'conda-forge::gcc_linux-64=12.3.0' \
#   'conda-forge::gfortran_linux-64=12.3.0' \
#   'conda-forge::openmpi=4.1.5' \
#   'conda-forge::sysroot_linux-64=2.17' \
#   'conda-forge::git=2.40.0' \
#   'conda-forge::git-lfs=3.3.0'
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# USER NEEDS TO SPECIFY THE FLAGS BELOW SO COCOA CAN FIND PYTHON/GCC/HDF5...
# --------------------------------------------------------------------------

export PYTHON_VERSION=3.8

export GLOBALPYTHON3="${CONDA_PREFIX:?}"/bin/python${PYTHON_VERSION:?}

export GLOBAL_PACKAGES_LOCATION="${CONDA_PREFIX:?}"

export GLOBALPIP3="${CONDA_PREFIX:?}"/bin/pip3

export GIT="${CONDA_PREFIX:?}/bin/git"

export WGET="${CONDA_PREFIX:?}/bin/wget"

export CURL="${CONDA_PREFIX:?}/bin/curl"

# --------------------------------------------------------------------------

INT_INCL="${CONDA_PREFIX:?}/include"

INT_LIB="${CONDA_PREFIX:?}/lib"

INT_INCL_PY="${INT_INCL:?}/python${PYTHON_VERSION:?}"

INT_INCL_PY_SP="${INT_INCL_PY}/site-packages"

INT_INCL_PY_SP_NP="${INT_INCL_PY_SP:?}/numpy/core/include"

# --------------------------------------------------------------------------

export PATH="${CONDA_PREFIX:?}"/bin:$PATH

export CFLAGS="${CFLAGS} -I${CONDA_PREFIX:?}/include"

export LDFLAGS="${LDFLAGS} -L${CONDA_PREFIX:?}/lib"

export C_INCLUDE_PATH="${INT_INCL:?}/":$C_INCLUDE_PATH

export C_INCLUDE_PATH="${INT_INCL_PY:?}m/":$C_INCLUDE_PATH

export C_INCLUDE_PATH="${INT_INCL_PY_SP_NP:?}/":$C_INCLUDE_PATH

export CPLUS_INCLUDE_PATH="${INT_INCL:?}/":$CPLUS_INCLUDE_PATH

export CPLUS_INCLUDE_PATH="${INT_INCL_PY:?}m/":$CPLUS_INCLUDE_PATH

export CPLUS_INCLUDE_PATH="${INT_INCL_PY_SP_NP:?}/":$CPLUS_INCLUDE_PATH

export PYTHONPATH="${INT_INCL_PY_SP:?}/":$PYTHONPATH

export PYTHONPATH="${INT_LIB:?}/":$PYTHONPATH

export LD_RUN_PATH="${INT_INCL_PY_SP}/":$LD_RUN_PATH

export LD_RUN_PATH="${INT_LIB:?}/":$LD_RUN_PATH

export LIBRARY_PATH="${INT_INCL_PY_SP}/":$LIBRARY_PATH

export LIBRARY_PATH="${INT_LIB:?}/":$LIBRARY_PATH

export CMAKE_INCLUDE_PATH="${INT_INCL:?}/":$CMAKE_INCLUDE_PATH

export CMAKE_INCLUDE_PATH="${INT_INCL_PY:?}m/":$CMAKE_INCLUDE_PATH    

export CMAKE_LIBRARY_PATH="${INT_INCL_PY_SP:?}/":$CMAKE_LIBRARY_PATH

export CMAKE_LIBRARY_PATH="${INT_LIB:?}":$CMAKE_LIBRARY_PATH

export INCLUDE_PATH="${INT_INCL:?}/":$INCLUDE_PATH

export INCLUDEPATH="${INT_INCL:?}/":$INCLUDEPATH

export INCLUDE="${CONDA_PREFIX:?}"/x86_64-conda-linux-gnu/include:$INCLUDE

export INCLUDE="${INT_INCL:?}/":$INCLUDE

export CPATH="${INT_INCL:?}/":${CPATH}

export OBJC_INCLUDE_PATH="${INT_INCL:?}/":$OBJC_INCLUDE_PATH

export OBJC_PATH="${CONDA_PREFIX:?}/include/":$OBJC_PATH

# --------------------------------------------------------------------------

unset -v INT_INCL INT_LIB INT_INCL_PY INT_INCL_PY_SP INT_INCL_PY_SP_NP

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
export C_COMPILER="${CONDA_PREFIX:?}"/bin/x86_64-conda-linux-gnu-cc

export CXX_COMPILER="${CONDA_PREFIX:?}"/bin/x86_64-conda-linux-gnu-g++

export FORTRAN_COMPILER="${CONDA_PREFIX:?}"/bin/x86_64-conda-linux-gnu-gfortran

export MPI_FORTRAN_COMPILER="${CONDA_PREFIX:?}"/bin/mpif90

export MPI_CC_COMPILER="${CONDA_PREFIX:?}"/bin/mpicc

export MPI_CXX_COMPILER="${CONDA_PREFIX:?}"/bin/mpicxx

# --------------------------------------------------------------------------
# FINE-TUNNING OVER THE USE OF SYSTEM-WIDE PACKAGES
# --------------------------------------------------------------------------
#export IGNORE_XZ_INSTALLATION=1

#export IGNORE_DISTUTILS_INSTALLATION=1

#export IGNORE_C_GSL_INSTALLATION=1

#export IGNORE_C_CFITSIO_INSTALLATION=1

#export IGNORE_C_FFTW_INSTALLATION=1

#export IGNORE_CPP_BOOST_INSTALLATION=1

#export IGNORE_CMAKE_INSTALLATION=1

#export IGNORE_OPENBLAS_INSTALLATION=1

#export IGNORE_FORTRAN_LAPACK_INSTALLATION=1

#export IGNORE_CPP_ARMA_INSTALLATION=1

#export IGNORE_HDF5_INSTALLATION=1

#export IGNORE_EXPAT_CORE_PACKAGE=1

#export IGNORE_PIP_CORE_PACKAGES=1

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------