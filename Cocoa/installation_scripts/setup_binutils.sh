#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_DISTUTILS_INSTALLATION}" ]; then
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
  if [ -z "${MAKE_NUM_THREADS}" ]; then
    pfail 'MAKE_NUM_THREADS'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${ptop}" || -z "${ptop2}" || -z "${pbottom}" || -z "${pbottom2}" ]; then
    pfail "PTOP/PBOTTOM"
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${DEBUG_DISTUTILS_PACKAGE}" ]; then
    export OUT1="/dev/null"
    export OUT2="/dev/null"
    export DISTUTILS_MAKE_NUM_THREADS="${MAKE_NUM_THREADS}"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
    export DISTUTILS_MAKE_NUM_THREADS=1
  fi
  unset_env_vars () {
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset DISTUTILS_MAKE_NUM_THREADS
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

  ptop2 'SETUP_BINUTILS'

  # ----------------------------------------------------------------------------
  # ----------------------------- TEXINFO LIBRARY  -----------------------------
  # ----------------------------------------------------------------------------
  ptop 'INSTALLING TEXINFO LIBRARY'

  if [ -z "${COCOA_TEXINFO_DIR}" ]; then
    pfail 'COCOA_TEXINFO_DIR'
    cd $ROOTDIR
    return 1
  fi

  cd $ROOTDIR/../cocoa_installation_libraries/$COCOA_TEXINFO_DIR

  FC=$FORTRAN_COMPILER CC=$C_COMPILER ./configure --prefix=$ROOTDIR/.local \
    --disable-perl-xs > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "CONFIGURE"
    return 1
  fi

  make -j $DISTUTILS_MAKE_NUM_THREADS all > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "MAKE ALL"
    return 1
  fi

  make install > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "MAKE INSTALL"
    return 1
  fi

  cd $ROOTDIR/
  pbottom 'INSTALLING TEXINFO LIBRARY'

  # ----------------------------------------------------------------------------
  # ----------------------------- DISTUTILS LIBRARY  ---------------------------
  # ----------------------------------------------------------------------------
  ptop 'INSTALLING BINUTILS LIBRARY'

  if [ -z "${COCOA_BINUTILS_DIR}" ]; then
    pfail 'COCOA_BINUTILS_DIR'
    cd $ROOTDIR
    return 1
  fi
  
  cd $ROOTDIR/../cocoa_installation_libraries/$COCOA_BINUTILS_DIR

  FC=$FORTRAN_COMPILER CC=$C_COMPILER ./configure --prefix=$ROOTDIR/.local > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "CONFIGURE"
    return 1
  fi

  make -j $DISTUTILS_MAKE_NUM_THREADS > ${OUT1} 2> ${OUT2}
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
  pbottom 'INSTALLING BINUTILS LIBRARY'

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  unset_env_vars
  pbottom2 'SETUP_BINUTILS'
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------