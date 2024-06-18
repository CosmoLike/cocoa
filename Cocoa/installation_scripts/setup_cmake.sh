#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_CMAKE_INSTALLATION}" ]; then
  pfail_scmk() {
    echo -e "\033[0;31m ERROR ENV VARIABLE ${1} IS NOT DEFINED \033[0m"
    unset pfail_scmk
  }
  if [ -z "${ROOTDIR}" ]; then
    pfail_scmk 'ROOTDIR'
    return 1
  fi
  if [ -z "${CXX_COMPILER}" ]; then
    pfail_scmk 'CXX_COMPILER'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${C_COMPILER}" ]; then
    pfail_scmk 'C_COMPILER'
    cd $ROOTDIR
    return 1
  fi
   if [ -z "${FORTRAN_COMPILER}" ]; then
    pfail_scmk 'FORTRAN_COMPILER'
    cd $ROOTDIR
    return 1
  fi 
  if [ -z "${DEBUG_CMAKE_PACKAGE}" ]; then
    if [ -z "${MAKE_NUM_THREADS}" ]; then
      pfail_scmk 'MAKE_NUM_THREADS'
      cd $ROOTDIR
      return 1
    fi 
    export OUT1="/dev/null"
    export OUT2="/dev/null"
    export CMAKE_MNT="${MAKE_NUM_THREADS}"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
    export CMAKE_MNT=1
  fi
  unset_env_vars_scmk () {
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset CMAKE_MNT
    unset pfail_scmk
    unset unset_env_vars_scmk
  }
  fail_scmk () {
    export MSG="\033[0;31m (setup_cmake.sh) WE CANNOT RUN \e[3m"
    export MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"
    unset_env_vars_scmk
    unset MSG
    unset MSG2
    unset fail_scmk
  }

  ptop2 'SETUP_CMAKE'
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  ptop 'INSTALLING CMAKE LIBRARY'

  if [ -z "${COCOA_CMAKE_DIR}" ]; then
    pfail_scmk 'COCOA_CMAKE_DIR'
    cd $ROOTDIR
    return 1
  fi

  cd $ROOTDIR/../cocoa_installation_libraries/$COCOA_CMAKE_DIR

  env CC=$C_COMPILER CXX=$CPP_COMPILER ./bootstrap \
    --prefix=$ROOTDIR/.local > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_scmk "BOOTSTRAP"
    return 1
  fi

  make -j $CMAKE_MNT > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_scmk "MAKE"
    return 1
  fi

  make install > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_scmk "INSTALL"
    return 1
  fi

  cd $ROOTDIR
  pbottom 'INSTALLING CMAKE LIBRARY'
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  unset_env_vars_scmk
  pbottom2 'SETUP_CMAKE'
fi
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------