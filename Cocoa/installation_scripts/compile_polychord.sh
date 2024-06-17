#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_POLYCHORD_COMPILATION}" ]; then
  # ----------------------------------------------------------------------------
  source $ROOTDIR/installation_scripts/clean_polychord.sh
  # ----------------------------------------------------------------------------
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
  if [ -z "${PIP3}" ]; then
    pfail 'PIP3'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${PYTHON3}" ]; then
    pfail "PYTHON3"
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${MAKE_NUM_THREADS}" ]; then
    pfail 'MAKE_NUM_THREADS'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${POLY_NAME}" ]; then
    pfail 'POLY_NAME'
    cd $ROOTDIR
    return 1
  fi
  unset_env_vars () {
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset POLY_MAKE_NUM_THREADS
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
  if [ -z "${DEBUG_POLY_OUTPUT}" ]; then
    export OUT1="/dev/null"
    export OUT2="/dev/null"
    export POLY_MAKE_NUM_THREADS="${MAKE_NUM_THREADS}"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
    export POLY_MAKE_NUM_THREADS=1
  fi

  # ---------------------------------------------------------------------------
  ptop2 'COMPILING POLYCHORD'

  cd $ROOTDIR/external_modules/code/$POLY_NAME/

  make all > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "MAKE ALL"
    return 1
  fi

  make -j $POLY_MAKE_NUM_THREADS pypolychord > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "MAKE PYPOLYCHORD"
    return 1
  fi

  CC=$MPI_CC_COMPILER CXX=$MPI_CXX_COMPILER $PYTHON3 setup.py install \
      --prefix $ROOTDIR/.local  > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "PYTHON3 SETUP INSTALL"
    return 1
  fi

  unset_env_vars
  pbottom2 'COMPILING POLYCHORD'
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------