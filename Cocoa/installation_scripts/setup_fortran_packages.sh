#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_FORTRAN_INSTALLATION}" ]; then
  pfail_sftrp() {
    echo -e "\033[0;31m ERROR ENV VARIABLE ${1} IS NOT DEFINED \033[0m"
    unset pfail_sftrp
  }
  if [ -z "${ROOTDIR}" ]; then
    pfail_sftrp 'ROOTDIR'
    return 1
  fi
  if [ -z "${CXX_COMPILER}" ]; then
    pfail_sftrp 'CXX_COMPILER'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${C_COMPILER}" ]; then
    pfail_sftrp 'C_COMPILER'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${FORTRAN_COMPILER}" ]; then
    pfail_sftrp 'FORTRAN_COMPILER'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${CMAKE}" ]; then
    pfail_sftrp 'CMAKE'
    cd $ROOTDIR
    return 1
  fi
  unset_env_vars_sftrp () {
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset FMNT
    unset pfail_sftrp
    unset unset_env_vars_sftrp
  }
  fail_sftrp () {
    export MSG="\033[0;31m (setup_fortran_packages.sh) WE CANNOT RUN \e[3m"
    export MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"
    unset_env_vars_sftrp
    unset MSG
    unset MSG2
    unset fail_sftrp
  }
  if [ -z "${DEBUG_FORTRAN_PACKAGES}" ]; then
    if [ -z "${MAKE_NUM_THREADS}" ]; then
      pfail_sftrp 'MAKE_NUM_THREADS'
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
  # ------------------------------- FORTRAN LAPACK -----------------------------
  # ----------------------------------------------------------------------------
  if [ -z "${IGNORE_FORTRAN_LAPACK_INSTALLATION}" ]; then
    ptop 'INSTALLING LAPACK FORTRAN LIBRARY'

    if [ -z "${COCOA_LAPACK_DIR}" ]; then
      pfail 'COCOA_LAPACK_DIR'
      cd $ROOTDIR
      return 1
    fi

    cd $ROOTDIR/../cocoa_installation_libraries/
    if [ $? -ne 0 ]; then
      fail_sftrp "CD COCOA_INSTALLATION_LIBRARIES"
      return 1
    fi

    rm -rf lapack-build
    
    mkdir lapack-build
    if [ $? -ne 0 ]; then
      fail_sftrp "MKDIR LAPACK-BUILD"
      return 1
    fi

    cd ./lapack-build

    $CMAKE -DBUILD_SHARED_LIBS=TRUE \
      -DCMAKE_INSTALL_PREFIX=$ROOTDIR/.local \
      -DCMAKE_C_COMPILER=$C_COMPILER \
      --log-level=ERROR \
      ../${COCOA_LAPACK_DIR} > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_sftrp "CMAKE"
      return 1
    fi

    make -j $FMNT all > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_sftrp "MAKE"
      return 1
    fi

    make install > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_sftrp "MAKE INSTALL"
      return 1
    fi
    
    cd $ROOTDIR    
    pbottom 'INSTALLING LAPACK FORTRAN LIBRARY'
  fi

  # ----------------------------------------------------------------------------
  # ---------------------------------------------------------------------------- 
  unset_env_vars_sftrp
  pbottom2 'SETUP_FORTRAN_PACKAGES'
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------