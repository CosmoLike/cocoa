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

export PYTHON_VERSION=3.8

export GLOBAL_PACKAGES_LOCATION="${CONDA_PREFIX:?}"

export GLOBALPYTHON3="${CONDA_PREFIX:?}"/bin/python${PYTHON_VERSION}

export GLOBALPIP3="${CONDA_PREFIX:?}"/bin/pip3

export GIT="${CONDA_PREFIX:?}"/bin/git

# --------------------------------------------------------------------------
# USER NEEDS TO SPECIFY THE FLAGS BELOW SO COCOA CAN FIND PYTHON/GCC/HDF5...
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
export PATH="${CONDA_PREFIX:?}"/bin:$PATH

export CFLAGS="${CFLAGS} -I"${CONDA_PREFIX:?}"/include"

export LDFLAGS="${LDFLAGS} -L"${CONDA_PREFIX:?}"/lib"

export C_INCLUDE_PATH="${CONDA_PREFIX:?}"/include:$C_INCLUDE_PATH

export C_INCLUDE_PATH="${CONDA_PREFIX:?}"/include/python${PYTHON_VERSION}m/:$C_INCLUDE_PATH

export CPLUS_INCLUDE_PATH="${CONDA_PREFIX:?}"/include:$CPLUS_INCLUDE_PATH

export CPLUS_INCLUDE_PATH="${CONDA_PREFIX:?}"/include/python${PYTHON_VERSION}m/:$CPLUS_INCLUDE_PATH

export PYTHONPATH="${CONDA_PREFIX:?}"/lib/python$PYTHON_VERSION/site-packages:$PYTHONPATH

export PYTHONPATH="${CONDA_PREFIX:?}"/lib:$PYTHONPATH

export LD_RUN_PATH="${CONDA_PREFIX:?}"/lib/python$PYTHON_VERSION/site-packages:$LD_RUN_PATH

export LD_RUN_PATH="${CONDA_PREFIX:?}"/lib:$LD_RUN_PATH

export LIBRARY_PATH="${CONDA_PREFIX:?}"/lib/python$PYTHON_VERSION/site-packages:$LIBRARY_PATH

export LIBRARY_PATH="${CONDA_PREFIX:?}"/lib:$LIBRARY_PATH

export CMAKE_INCLUDE_PATH="${CONDA_PREFIX:?}"/include/:$CMAKE_INCLUDE_PATH

export CMAKE_INCLUDE_PATH="${CONDA_PREFIX:?}"/include/python${PYTHON_VERSION}m/:$CMAKE_INCLUDE_PATH    

export CMAKE_LIBRARY_PATH="${CONDA_PREFIX:?}"/lib/python$PYTHON_VERSION/site-packages:$CMAKE_LIBRARY_PATH

export CMAKE_LIBRARY_PATH="${CONDA_PREFIX:?}"/lib:$CMAKE_LIBRARY_PATH

export INCLUDE_PATH="${CONDA_PREFIX:?}"/include/:$INCLUDE_PATH

export INCLUDEPATH="${CONDA_PREFIX:?}"/include/:$INCLUDEPATH

export INCLUDE="${CONDA_PREFIX:?}"/x86_64-conda-linux-gnu/include:$INCLUDE

export INCLUDE="${CONDA_PREFIX:?}"/include/:$INCLUDE

export CPATH="${CONDA_PREFIX:?}"/include/:$CPATH

export OBJC_INCLUDE_PATH="${CONDA_PREFIX:?}"/include/:$OBJC_INCLUDE_PATH

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
# FINE-TUNNING OVER THE USE OF SYSTEM-WIDE PACKAGES
# --------------------------------------------------------------------------
#export IGNORE_XZ_INSTALLATION=1

#export IGNORE_DISTUTILS_INSTALLATION=1

#export IGNORE_CMAKE_INSTALLATION=1

#export IGNORE_HDF5_INSTALLATION=1

#export IGNORE_C_GSL_INSTALLATION=1

#export IGNORE_C_CFITSIO_INSTALLATION=1

#export IGNORE_C_FFTW_INSTALLATION=1

#export IGNORE_OPENBLAS_INSTALLATION=1

#exportIGNORE_FORTRAN_LAPACK_INSTALLATION

#export IGNORE_CPP_BOOST_INSTALLATION=1

#export IGNORE_CPP_ARMA_INSTALLATION=1

# ------------------------------------------------------------------------------
# ------------------------------- PACKAGES LOCATION ----------------------------
# ------------------------------------------------------------------------------
export COCOA_ARMADILLO_DIR=armadillo-12.8.2/

export COCOA_BOOST_DIR=boost_1_81_0/

export COCOA_EXPAT_DIR=expat-2.5.0/

export COCOA_XZ_DIR=xz-5.2.5/

export COCOA_XZ_FILE=xz-5.2.5.tar.gz

export COCOA_CMAKE_DIR=cmake-3.26.4/

export COCOA_BINUTILS_DIR=binutils-2.37/

export COCOA_TEXINFO_DIR=texinfo-7.0.3/

export COCOA_OPENBLAS_DIR=OpenBLAS-0.3.23/

export COCOA_LAPACK_DIR=lapack-3.11.0/

export COCOA_HDF5_DIR=hdf5-1.12.3/

export COCOA_CFITSIO_DIR=cfitsio-4.0.0/

export COCOA_FFTW_DIR=fftw-3.3.10/

export COCOA_GSL_DIR=gsl-2.7/

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------