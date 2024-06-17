#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_OPENBLAS_INSTALLATION}" ]; then
  pfail() {
    echo -e "\033[0;31m ERROR ENV VARIABLE ${1} IS NOT DEFINED \033[0m"
    unset pfail
  }
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'
    return 1
  fi
  if [ -z "${CXX_COMPILER}" ]; then
    pfail 'CXX_COMPILER'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${C_COMPILER}" ]; then
    pfail 'C_COMPILER'
    cd $ROOTDIR
    return 1
  fi
   if [ -z "${FORTRAN_COMPILER}" ]; then
    pfail 'FORTRAN_COMPILER'
    cd $ROOTDIR
    return 1
  fi

  ptop2 'SETUP_OPENBLAS'
  ptop  'INSTALLING OPENBLAS LIBRARY'
  
  unset_env_vars () {
    cd $ROOTDIR
    unset MAKE_NB_JOBS
    unset OUT1
    unset OUT2
    unset OPENBLAS_MNT
  }

  if [ -z "${DEBUG_OPENBLAS_PACKAGE}" ]; then
    export OUT1="/dev/null"
    export OUT2="/dev/null"
    export OPENBLAS_MNT="${MAKE_NUM_THREADS}"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
    export OPENBLAS_MNT=1
  fi

  cd $ROOTDIR/../cocoa_installation_libraries/$COCOA_OPENBLAS_DIR
  
  export MAKE_NB_JOBS=$OPENBLAS_MNT
  
  make clean > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "MAKE CLEAN"
    return 1
  fi

  make CC=$C_COMPILER FC=$FORTRAN_COMPILER USE_OPENMP=1 > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "MAKE INSTALL"
    return 1
  fi

  make install PREFIX=$ROOTDIR/.local > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "MAKE INSTALL"
    return 1
  fi

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  unset_env_vars
  unset unset_env_vars
  pbottom  'INSTALLING OPENBLAS LIBRARY'
  pbottom2 'SETUP_OPENBLAS'
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------