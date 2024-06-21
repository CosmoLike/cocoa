source ./installation_scripts/impl_unset_keys.sh
export ROOTDIR=$(pwd -P)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ----------------------- HOW COCOA SHOULD BE INSTALLED? -----------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
export MINICONDA_INSTALLATION=1
#export MANUAL_INSTALLATION=1

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ---------------------------------- GLOBAL SETTINGS ---------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# --------------------------- VERBOSE AS DEBUG TOOL ----------------------------
# ------------------------------------------------------------------------------
#export COCOA_OUTPUT_VERBOSE=1

# ------------------------------------------------------------------------------
# SKIP DOWNLOADING/DECOMPRESSING DATASETs (CAN SAVE TIME WHEN INSTALLING COCOA)
# ------------------------------------------------------------------------------
# export SKIP_DECOMM_ACT
# export SKIP_DECOMM_SPT
# export SKIP_DECOMM_PLANCK
# export SKIP_DECOMM_BICEP
# export SKIP_DECOMM_STRONG_LENSING
# export SKIP_DECOMM_SN
# export SKIP_DECOMM_BAO
export SKIP_DECOMM_SIMONS_OBSERVATORY=1
export SKIP_DECOMM_CAMSPEC=1
export SKIP_DECOMM_LIPOP=1

# ------------------------------------------------------------------------------
# ----- IF DEFINED, COSMOLIKE WILL BE COMPILED WITH DEBUG FLAG -----------------
# ----- DEBUG FLAG = ALL COMPILER WARNINGS + NO MATH OPTIMIZATION + NO OPENMP --
# ------------------------------------------------------------------------------
export COSMOLIKE_DEBUG_MODE=1

# ------------------------------------------------------------------------------
# ----- IF TRUE, THEN COCOA USES CLIK FROM https://github.com/benabed/clik -----
# ------------------------------------------------------------------------------
export USE_SPT_CLIK_PLANCK=1

# ------------------------------------------------------------------------------
# ----------------- CONTROL OVER THE COMPILATION OF EXTERNAL CODES -------------
# ------------------------------------------------------------------------------
#export IGNORE_CAMB_COMPILATION=1
#export IGNORE_CLASS_COMPILATION=1
#export IGNORE_COSMOLIKE_COMPILATION=1
#export IGNORE_POLYCHORD_COMPILATION=1
#export IGNORE_PLANCK_COMPILATION=1
#export IGNORE_ACT_COMPILATION=1
#export IGNORE_COBAYA_INSTALLATION=1
#export IGNORE_CAMSPEC_INSTALLATION=1
#export IGNORE_LIPOP_INSTALLATION=1
#export IGNORE_SO_INSTALLATION=1

# ------------------------------------------------------------------------------
# THREADING COMPILATION/INSTALLATION OF LIBRARIES
# ------------------------------------------------------------------------------
export MAKE_NUM_THREADS=4

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# PACKAGE VERSIONS. BE CAREFUL, CHANGES IN THE COMMIT ID MAY BREAK COCOA
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
export COBAYA_GIT_COMMIT="2636ea9ed399c35c5d276de1acb15aaafbcab10c"
export HILLIPOP_GIT_COMMIT="cc9cbe31991d4662522241543a46d44d2cdec251"
export LOLLIPOP_GIT_COMMIT="280a9c93d33bc6a058d6bf769ec82d9f7fdbd2b6"
export LIPOP_DATA_VERSION=4.2
export SPT3G_DATA_GIT_COMMIT="66da8e9e2f325024566fe13245788bf8ede897bc"
export HOLICOW_DATA_GIT_COMMIT="f792647d1fd6c09d9e052fef526669cbd702ab82"
export POLYCHORD_GIT_COMMIT="daba49d1385d065122db76a2b384050f9e95d278"
export CAMB_GIT_COMMIT="45d1c3d27e7480c0f9a82c98522c17ed422dd408"
export CLASS_GIT_COMMIT="8df566c1ff2d0b3e40e106567c435575aea337be"

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ----------------------------- END OF OPTIONS ---------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# SET ENV FLAGS FOR LINKING/RUNNING/COMPILATION OF PROGRAMS --------------------
# ------------------------------------------------------------------------------
if [ -n "${MANUAL_INSTALLATION}" ]; then
    
    # --------------------------------------------------------------------------
    # IF SET, THEN COCOA ADOPTS FFTW10. OTHERWISE, COCOA ADOPTS FFTW8
    # --------------------------------------------------------------------------
    #export FFTW_NEW_VERSION=1
    
    # --------------------------------------------------------------------------
    # IF SET, COCOA DOES NOT USE SYSTEM PIP PACKAGES 
    # --------------------------------------------------------------------------
    export DONT_USE_SYSTEM_PIP_PACKAGES=1
    
    # --------------------------------------------------------------------------
    # IF NOT SET, COCOA WILL INSTALL TENSORFLOW, KERAS, PYTORCH
    # --------------------------------------------------------------------------
    export IGNORE_EMULATOR_CPU_PIP_PACKAGES=1
    export IGNORE_EMULATOR_GPU_PIP_PACKAGES=1

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
    
    export OBJC_INCLUDE_PATH="${CONDA_PREFIX:?}"/include/:OBJC_INCLUDE_PATH
    
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
# ------------------------------------------------------------------------------
elif [ -n "${MINICONDA_INSTALLATION}" ]; then
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

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

    # --------------------------------------------------------------------------
    # IF NOT SET, COCOA WILL INSTALL TENSORFLOW, KERAS, AND PYTORCH 
    # --------------------------------------------------------------------------
    export IGNORE_EMULATOR_CPU_PIP_PACKAGES=1
    
    export IGNORE_EMULATOR_GPU_PIP_PACKAGES=1
fi

ulimit -s unlimited
export COBAYA_PACKAGES_PATH=$ROOTDIR/external_modules
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_CMAKE_INSTALLATION}" ]; then
    export CMAKE_ROOT=${ROOTDIR:?}/.local/bin/cmake
    export CMAKE=${ROOTDIR:?}/.local/bin/cmake
else
    export CMAKE=cmake
fi
# ------------------------------------------------------------------------------
if [ -n "${IGNORE_CPP_INSTALLATION}" ]; then
  export IGNORE_CPP_BOOST_INSTALLATION=1
  export IGNORE_CPP_ARMA_INSTALLATION=1
  export IGNORE_CPP_SPDLOG_INSTALLATION=1
  export IGNORE_CPP_CARMA_INSTALLATION=1
fi
# ------------------------------------------------------------------------------
if [ -n "${IGNORE_C_INSTALLATION}" ]; then
  export IGNORE_C_CFITSIO_INSTALLATION=1
  export IGNORE_C_FFTW_INSTALLATION=1
  export IGNORE_C_GSL_INSTALLATION=1
fi
# ------------------------------------------------------------------------------
if [ -n "${IGNORE_FORTRAN_INSTALLATION}" ]; then
  export IGNORE_FORTRAN_LAPACK_INSTALLATION=1
fi
# ------------------------------------------------------------------------------
if [ -n "${GLOBAL_PACKAGES_LOCATION}" ]; then
  export GLOBAL_PACKAGES_INCLUDE=$GLOBAL_PACKAGES_LOCATION/include
  export GLOBAL_PACKAGES_LIB=$GLOBAL_PACKAGES_LOCATION/lib
fi
# ------------------------------------------------------------------------------
if [ -z "${THREAD_UNXZ}" ]; then
  export MAKE_NUM_THREADS=1
fi
# ------------------------------------------------------------------------------
export PYTHON3="${ROOTDIR}/.local/bin/python3"
export PIP3="${PYTHON3:?} -m pip"
export COBAYA_PACKAGES_PATH=external_modules
# ------------------------------------------------------------------------------
if [ -n "${COSMOLIKE_DEBUG_MODE}" ]; then
    export SPDLOG_LEVEL=debug
else
    export SPDLOG_LEVEL=info
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
export COCOA_RUN_EVALUATE="mpirun -n 1 --oversubscribe --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by core --map-by numa:pe=4 cobaya-run"

export COCOA_RUN_MCMC="mpirun -n 4 --oversubscribe --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by core --map-by numa:pe=4 cobaya-run"

# ------------------------------------------------------------------------------
# DEBUG THE COMPILATION OF PREREQUISITES PACKAGES. IF YOU NEED TO RUN ----------
# SETUP_COCOA_INSTALLATION_PACKAGES >1x AND WANT TO SKIP -----------------------
# PACKAGE FILE DECOMPRESSION ---------------------------------------------------
# ------------------------------------------------------------------------------
#export DEBUG_SKIP_FILE_DECOMPRESSION_SETUP_COCOA=1

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------- PACKAGES LOCATION ----------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
export COCOA_SPDLOG_DIR=spdlog/
export COCOA_ARMADILLO_DIR=armadillo-12.8.2/
export COCOA_BOOST_DIR=boost_1_81_0/
export COCOA_EXPAT_DIR=expat-2.5.0/
export COCOA_XZ_DIR=xz-5.2.5/
export COCOA_XZ_FILE=xz-5.2.5.tar.gz
export COCOA_CARMA_DIR=carma/
export COCOA_CMAKE_DIR=cmake-3.26.4/
export COCOA_BINUTILS_DIR=binutils-2.37/
export COCOA_TEXINFO_DIR=texinfo-7.0.3/
export COCOA_OPENBLAS_DIR=OpenBLAS-0.3.23/
export COCOA_LAPACK_DIR=lapack-3.11.0/
export COCOA_HDF5_DIR=hdf5-1.12.3/
export COCOA_CFITSIO_DIR=cfitsio-4.0.0/
export COCOA_FFTW_DIR=fftw-3.3.10/
export COCOA_GSL_DIR=gsl-2.7/
export CAMB_NAME='CAMB'
export POLY_NAME="PolyChordLite"
export CLASS_NAME="class_public"
export ACT_NAME="pyactlike"
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------