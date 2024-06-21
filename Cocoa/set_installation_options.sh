#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
source ./installation_scripts/impl_unset_keys.sh
export ROOTDIR=$(pwd -P)

# ------------------------------------------------------------------------------
# ----------------------- HOW COCOA SHOULD BE INSTALLED? -----------------------
# ------------------------------------------------------------------------------
export MINICONDA_INSTALLATION=1
#export MANUAL_INSTALLATION=1

# ------------------------------------------------------------------------------
# --------------------------- VERBOSE AS DEBUG TOOL ----------------------------
# ------------------------------------------------------------------------------
#export COCOA_OUTPUT_VERBOSE=1

# ------------------------------------------------------------------------------
# --- SKIP DOWNLOAD/DECOMPRESSING DATASET (SAVE TIME WHEN INSTALLING COCOA) ----
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
# ---- IF DEFINED, COSMOLIKE WILL BE COMPILED WITH DEBUG FLAG ------------------
# ---- DEBUG FLAG = ALL COMPILER WARNINGS + NO MATH OPTIMIZATION + NO OPENMP ---
# ------------------------------------------------------------------------------
export COSMOLIKE_DEBUG_MODE=1

# ------------------------------------------------------------------------------
# ----- IF TRUE, THEN COCOA USES CLIK FROM https://github.com/benabed/clik -----
# ------------------------------------------------------------------------------
export USE_SPT_CLIK_PLANCK=1

# ------------------------------------------------------------------------------
# --------------- CONTROL OVER THE COMPILATION OF EXTERNAL CODES ---------------
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
# THREADING COMPILATION/INSTALLATION OF LIBRARIES ------------------------------
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

# --------------------------------------------------------------------------
# IF NOT SET, COCOA WILL INSTALL TENSORFLOW, KERAS, AND PYTORCH 
# --------------------------------------------------------------------------
export IGNORE_EMULATOR_CPU_PIP_PACKAGES=1

export IGNORE_EMULATOR_GPU_PIP_PACKAGES=1

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
      
  # NEED TO ADJUST THE FLAGS IN THE FILE BELOW IF MANUAL INSTALLATION  
  source .flags_manual_installation.sh || return 1

elif [ -n "${MINICONDA_INSTALLATION}" ]; then

  source .flags_miniconda_installation.sh || return 1

fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

source .flags_derived.sh || return 1

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
# ------------------------------- PACKAGES LOCATION ----------------------------
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
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
export CAMB_NAME='CAMB'
export POLY_NAME="PolyChordLite"
export CLASS_NAME="class_public"
export ACT_NAME="pyactlike"
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------