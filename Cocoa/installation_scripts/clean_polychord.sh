#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_POLYCHORD_COMPILATION}" ]; then
  pfail() {
    echo -e "\033[0;31m ERROR ENV VARIABLE ${1} IS NOT DEFINED \033[0m"
    unset pfail
  }
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'
    return 1
  fi
  if [ -z "${PYTHON3}" ]; then
    pfail "PYTHON3"
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${PYTHON_VERSION}" ]; then
    pfail "PYTHON3"
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${POLY_NAME}" ]; then
    pfail 'POLY_NAME'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${ptop}" || -z "${ptop2}" || -z "${pbottom}" || -z "${pbottom2}" ]; then
    pfail "PTOP/PBOTTOM"
    cd $ROOTDIR
    return 1
  fi
  unset_env_vars () {
    cd $ROOTDIR
    unset OUT1
    unset OUT2
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
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
  fi

  # ---------------------------------------------------------------------------
  ptop 'CLEANING POLYCHORD'

  rm -rf $ROOTDIR/.local/lib/python$PYTHON_VERSION/site-packages/pypolychord-*

  cd $ROOTDIR/external_modules/code/$POLY_NAME

  make clean > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "MAKE CLEAN"
    return 1
  fi

  rm -rf ./lib/*.a
  rm -rf ./lib/*.so

  unset_env_vars
  pbottom 'CLEANING POLYCHORD'
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------