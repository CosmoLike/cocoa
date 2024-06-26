#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

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
if [ -n "${DEBUG_SKIP_FILE_DECOMPRESSION_SETUP_COCOA}" ]; then
  export SKIP_DECOMM_CORE_PACKAGES=1
  export SKIP_DECOMM_ACT=1
  export SKIP_DECOMM_SPT=1
  export SKIP_DECOMM_PLANCK=1
  export SKIP_DECOMM_BICEP=1
  export SKIP_DECOMM_STRONG_LENSING=1
  export SKIP_DECOMM_SN=1
  export SKIP_DECOMM_BAO=1
  export SKIP_DECOMM_SIMONS_OBSERVATORY=1
  export SKIP_DECOMM_CAMSPEC=1
  export SKIP_DECOMM_LIPOP=1
fi

# ------------------------------------------------------------------------------
if [ -n "${GLOBAL_PACKAGES_LOCATION}" ]; then
  export GLOBAL_PACKAGES_INCLUDE=$GLOBAL_PACKAGES_LOCATION/include
  export GLOBAL_PACKAGES_LIB=$GLOBAL_PACKAGES_LOCATION/lib
fi

# ------------------------------------------------------------------------------
export PYTHON3="${ROOTDIR}/.local/bin/python3"
export PIP3="${PYTHON3:?} -m pip"
export COBAYA_PACKAGES_PATH=external_modules

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -n "${COSMOLIKE_DEBUG_MODE}" ]; then
    export SPDLOG_LEVEL=debug
else
    export SPDLOG_LEVEL=info
fi

# ----------------------------------------------------------------------------
# ---------------------------- GLOBAL FUNCTIONS ------------------------------
# ----------------------------------------------------------------------------

ptop() {
  echo -e "\033[1;34m  ${1} \033[0m"
}

pbottom() {
  echo -e "\033[1;34m  \e[4m${1} DONE\033[0m"
}

ptop2() {
  echo -e "\033[1;44m${1} \033[0m"
}

pbottom2() {
  echo -e "\033[1;44m\e[4m${1} DONE\033[0m"
}

pfail() {
  echo -e \
  "\033[0;31m       ERROR ENV VARIABLE ${1:-"empty arg"} NOT DEFINED \033[0m"
}

cdroot() {
  cd "${ROOTDIR:?}" 2>"/dev/null" || { echo -e \
    "\033[0;31m\t\t CD ROOTDIR (${ROOTDIR}) FAILED \033[0m"; return 1; }
}

if [ -z "${COCOA_OUTPUT_VERBOSE}" ]; then
  export OUT1="/dev/null"; export OUT2="/dev/null"
  export MNT="${MAKE_NUM_THREADS:-1}"; 
  [[ ${MNT} == +([0-9]) ]] || export MNT=1
else
  export OUT1="/dev/tty"; export OUT2="/dev/tty"
  export MNT=1
fi

fail_script_msg () {
  local MSG="\033[0;31m        (${1:-"empty arg"}) we cannot run \e[3m"
  local MSG2="\033[0m"
  echo -e "${MSG} ${2:-"empty arg"} ${MSG2}"
}

# ----------------------------------------------------------------------------
# ------------------------------- ERROR CODES --------------------------------
# ----------------------------------------------------------------------------
export EC1="PYTHON SETUP CLEAN"

export EC2="MAKE CLEAN"

export EC3="PIP3 INSTALL"

export EC4="PYTHON3 SETUP.PY BUILD"

export EC5="PYTHON3 WAF CONFIGURE"

export EC6="PYTHON3 WAF INSTALL"

export EC7="MAKE ALL"

export EC8="MAKE"

export EC9="PYTHON3 SETUP INSTALL"

export EC10="MAKE INSTALL"

export EC11="CONFIGURE"

export EC12="CMAKE"

export EC13="PIP INSTALL"

export EC14="MKDIR BUILD FOLDER"

export EC15="GIT CLONE"

export EC16="GIT CHECKOUT"

export EC17="SHELL SCRIPT FILE PATCH"

export EC18="PYTHON3 WAF DISTCLEAN"

export EC19="BOOTSTRAP"

export EC20="MKDIR"

export EC21="B2"

export EC22="SHELL SCRIPT"

export EC23="GIT RESET"

export EC24="WGET"

export EC25="TAR (DECOMPRESS)"

export EC26="UNZIP (DECOMPRESS)"

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

export COCOA_RUN_EVALUATE="mpirun -n 1 --oversubscribe --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by core --map-by numa:pe=4 cobaya-run"

export COCOA_RUN_MCMC="mpirun -n 4 --oversubscribe --mca btl vader,tcp,self --bind-to core:overload-allowed --rank-by core --map-by numa:pe=4 cobaya-run"

# ------------------------------------------------------------------------------
# DEBUG THE COMPILATION OF PREREQUISITES PACKAGES. IF YOU NEED TO RUN ----------
# SETUP_COCOA_INSTALLATION_PACKAGES >1x AND WANT TO SKIP FILE DECOMPRESSION ----
# ------------------------------------------------------------------------------
#export DEBUG_SKIP_FILE_DECOMPRESSION_SETUP_COCOA=1

# ------------------------------------------------------------------------------
# -------- PACKAGE LOCATION (DEBUG/RARELY MODIFIED FLAGS) ----------------------
# ------------------------------------------------------------------------------
export COCOA_SPDLOG_DIR=spdlog/

export COCOA_CARMA_DIR=carma/
