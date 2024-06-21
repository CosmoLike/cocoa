#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

export PYTHON_VERSION=3.8

export GLOBALPYTHON3="${CONDA_PREFIX:?}"/bin/python${PYTHON_VERSION:?}

export GLOBAL_PACKAGES_LOCATION="${CONDA_PREFIX:?}"

export GLOBALPIP3="${CONDA_PREFIX:?}"/bin/pip3

export GIT="${CONDA_PREFIX:?}"/bin/git

export PATH="${CONDA_PREFIX:?}"/bin:$PATH

export INT_INCL="${CONDA_PREFIX:?}/include"

export INT_LIB="${CONDA_PREFIX:?}/lib"

export INT_INCL_PY="${INT_INCL:?}/python${PYTHON_VERSION}"

export INT_INCL_PY_SP="${INT_INCL_PY}/site-packages"

export CFLAGS="${CFLAGS} -I${CONDA_PREFIX:?}/include"

export LDFLAGS="${LDFLAGS} -L${CONDA_PREFIX:?}/lib"

export C_INCLUDE_PATH="${INT_INCL:?}":$C_INCLUDE_PATH

export C_INCLUDE_PATH="${INT_INCL_PY:?}m/":$C_INCLUDE_PATH

export C_INCLUDE_PATH="${INT_INCL_PY_SP:?}/numpy/core/include/":$C_INCLUDE_PATH

export CPLUS_INCLUDE_PATH="${INT_INCL:?}":$CPLUS_INCLUDE_PATH

export CPLUS_INCLUDE_PATH="${INT_INCL_PY:?}m/":$CPLUS_INCLUDE_PATH

export CPLUS_INCLUDE_PATH="${INT_INCL_PY_SP:?}/numpy/core/include/":$CPLUS_INCLUDE_PATH

export PYTHONPATH="${INT_INCL_PY_SP:?}/":$PYTHONPATH

export PYTHONPATH="${INT_LIB:?}/":$PYTHONPATH

export LD_RUN_PATH="${INT_INCL_PY:?}/site-packages":$LD_RUN_PATH

export LD_RUN_PATH="${INT_LIB:?}/":$LD_RUN_PATH

export LIBRARY_PATH="${INT_INCL_PY:?}"/site-packages:$LIBRARY_PATH

export LIBRARY_PATH="${INT_LIB:?}/":$LIBRARY_PATH

export CMAKE_INCLUDE_PATH="${INT_INCL:?}"/:$CMAKE_INCLUDE_PATH

export CMAKE_INCLUDE_PATH="${INT_INCL_PY:?}m/":$CMAKE_INCLUDE_PATH    

export CMAKE_LIBRARY_PATH="${INT_INCL_PY:?}"/site-packages:$CMAKE_LIBRARY_PATH

export CMAKE_LIBRARY_PATH="${INT_LIB:?}":$CMAKE_LIBRARY_PATH

export INCLUDE_PATH="${INT_INCL:?}/":$INCLUDE_PATH

export INCLUDEPATH="${INT_INCL:?}/":$INCLUDEPATH

export INCLUDE="${CONDA_PREFIX:?}"/x86_64-conda-linux-gnu/include:$INCLUDE

export INCLUDE="${INT_INCL:?}/":$INCLUDE

export CPATH="${INT_INCL:?}/":${CPATH}

export OBJC_INCLUDE_PATH="${INT_INCL:?}/":OBJC_INCLUDE_PATH

export OBJC_PATH="${CONDA_PREFIX:?}"/include/:OBJC_PATH

# --------------------------------------------------------------------------
# COMPILER
# --------------------------------------------------------------------------
export C_COMPILER="${CONDA_PREFIX:?}"/bin/x86_64-conda-linux-gnu-cc

export CXX_COMPILER="${CONDA_PREFIX:?}"/bin/x86_64-conda-linux-gnu-g++

export FORTRAN_COMPILER="${CONDA_PREFIX:?}"/bin/x86_64-conda-linux-gnu-gfortran

export MPI_FORTRAN_COMPILER="${CONDA_PREFIX:?}"/bin/mpif90

export MPI_CC_COMPILER="${CONDA_PREFIX:?}"/bin/mpicc

export MPI_CXX_COMPILER="${CONDA_PREFIX:?}"/bin/mpicxx

# --------------------------------------------------------------------------
# IGNORE MOST PACKAGES (ALREADY ON CONDA)
# --------------------------------------------------------------------------
export IGNORE_XZ_INSTALLATION=1

export IGNORE_DISTUTILS_INSTALLATION=1

export IGNORE_C_GSL_INSTALLATION=1

export IGNORE_C_CFITSIO_INSTALLATION=1

export IGNORE_C_FFTW_INSTALLATION=1

export IGNORE_CPP_BOOST_INSTALLATION=1

export IGNORE_CMAKE_INSTALLATION=1

export IGNORE_OPENBLAS_INSTALLATION=1

export IGNORE_FORTRAN_LAPACK_INSTALLATION=1

export IGNORE_CPP_ARMA_INSTALLATION=1

export IGNORE_HDF5_INSTALLATION=1

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------