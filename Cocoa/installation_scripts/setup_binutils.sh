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
  if [ -z "${DEBUG_DISTUTILS_PACKAGE}" ]; then
    export OUT1="/dev/null"
    export OUT2="/dev/null"
    export DISTUTILS_MAKE_NUM_THREADS="${MAKE_NUM_THREADS}"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
    export DISTUTILS_MAKE_NUM_THREADS=1
  fi
  unset_env_vars_sbinutils () {
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset DISTUTILS_MAKE_NUM_THREADS
    unset pfail
    unset unset_env_vars_sbinutils
  }
  fail_sbinutils () {
    export MSG="\033[0;31m (setup_binutils.sh) WE CANNOT RUN \e[3m"
    export MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"
    unset_env_vars_sbinutils
    unset MSG
    unset MSG2
    unset fail_sbinutils
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
    fail_sbinutils "CONFIGURE"
    return 1
  fi

  make -j $DISTUTILS_MAKE_NUM_THREADS all > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_sbinutils "MAKE ALL"
    return 1
  fi

  make install > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_sbinutils "MAKE INSTALL"
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
    fail_sbinutils "CONFIGURE"
    return 1
  fi

  make -j $DISTUTILS_MAKE_NUM_THREADS > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_sbinutils "MAKE"
    return 1
  fi

  make install > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_sbinutils "MAKE INSTALL"
    return 1
  fi

  cd $ROOTDIR
  pbottom 'INSTALLING BINUTILS LIBRARY'
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  unset_env_vars_sbinutils
  pbottom2 'SETUP_BINUTILS'
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------