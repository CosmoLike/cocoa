#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_CMAKE_INSTALLATION}" ]; then
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
  if [ -z "${DEBUG_CMAKE_PACKAGE}" ]; then
    if [ -z "${MAKE_NUM_THREADS}" ]; then
      pfail 'MAKE_NUM_THREADS'
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
  unset_env_vars () {
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset CMAKE_MNT
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

  ptop2 'SETUP_CMAKE'
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  ptop 'INSTALLING CMAKE LIBRARY'

  if [ -z "${COCOA_CMAKE_DIR}" ]; then
    pfail 'COCOA_CMAKE_DIR'
    cd $ROOTDIR
    return 1
  fi

  cd $ROOTDIR/../cocoa_installation_libraries/$COCOA_CMAKE_DIR

  env CC=$C_COMPILER CXX=$CPP_COMPILER ./bootstrap \
    --prefix=$ROOTDIR/.local > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "BOOTSTRAP"
    return 1
  fi

  make -j $CMAKE_MNT > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "MAKE"
    return 1
  fi

  make install > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "INSTALL"
    return 1
  fi

  cd $ROOTDIR
  pbottom 'INSTALLING CMAKE LIBRARY'
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  unset_env_vars
  pbottom2 'SETUP_CMAKE'
fi
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------