#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_FORTRAN_INSTALLATION}" ]; then
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
  if [ -z "${CMAKE}" ]; then
    pfail 'CMAKE'
    cd $ROOTDIR
    return 1
  fi
  unset_env_vars () {
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset FMNT
    unset pfail
    unset unset_env_vars
  }
  fail () {
    export FAILMSG="\033[0;31m WE CANNOT RUN \e[3m"
    export FAILMSG2="\033[0m"
    echo -e "${FAILMSG} ${1} ${FAILMSG2}"
    unset_env_vars
    unset FAILMSG
    unset FAILMSG2
    unset fail
  }
  if [ -z "${DEBUG_FORTRAN_PACKAGES}" ]; then
    if [ -z "${MAKE_NUM_THREADS}" ]; then
      pfail 'MAKE_NUM_THREADS'
      cd $ROOTDIR
      return 1
    fi
    export OUT1="/dev/null"
    export OUT2="/dev/null"
    export FMNT="${MAKE_NUM_THREADS}"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
    export FMNT=1
  fi

  ptop2 'SETUP_FORTRAN_PACKAGES'

  # ----------------------------------------------------------------------------
  # ---------------------------------- FORTRAN LAPACK --------------------------
  # ----------------------------------------------------------------------------
  if [ -z "${IGNORE_FORTRAN_LAPACK_INSTALLATION}" ]; then
    ptop 'INSTALLING LAPACK FORTRAN LIBRARY'

    cd $ROOTDIR/../cocoa_installation_libraries/
    if [ $? -ne 0 ]; then
      fail "CD COCOA_INSTALLATION_LIBRARIES"
      return 1
    fi

    rm -rf lapack-build
    
    mkdir lapack-build
    if [ $? -ne 0 ]; then
      fail "MKDIR LAPACK-BUILD"
      return 1
    fi

    cd ./lapack-build

    $CMAKE -DBUILD_SHARED_LIBS=TRUE \
      -DCMAKE_INSTALL_PREFIX=$ROOTDIR/.local \
      -DCMAKE_C_COMPILER=$C_COMPILER \
      --log-level=ERROR \
      ../${COCOA_LAPACK_DIR} > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail "CMAKE"
      return 1
    fi

    make -j $FMNT all > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail "MAKE"
      return 1
    fi

    make install > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail "MAKE INSTALL"
      return 1
    fi
    
    cd $ROOTDIR    
    pbottom 'INSTALLING LAPACK FORTRAN LIBRARY'
  fi

  # ----------------------------------------------------------------------------
  # ---------------------------------------------------------------------------- 
  unset_env_vars
  pbottom2 'SETUP_FORTRAN_PACKAGES'
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------