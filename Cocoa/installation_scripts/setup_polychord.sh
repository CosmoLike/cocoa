#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_POLYCHORD_COMPILATION}" ]; then
  echo -e '\033[1;34m''\tINSTALLING CAMB VIA SETUP_CAMB''\033[0m'

  if [ -z "${ROOTDIR}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE ROOTDIR IS NOT DEFINED''\033[0m'
    return 1
  fi
  if [ -z "${GIT}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE GIT IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${CAMB_GIT_COMMIT}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE POLYCHORD_GIT_COMMIT IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
  fi

  unset_env_vars () {
    cd $ROOTDIR
    unset OUT_POLY_1
    unset OUT_POLY_2
    unset POLYCHORD_URL
    unset POLY_NAME
    unset FAILMSG
    unset FAILMSG2
    unset fail
    unset unset_env_vars
  }

  fail () {
    export FAILMSG="\033[0;31m SETUP POLYCHORD COULD NOT RUN \e[3m"
    export FAILMSG2="\033[0m"
    
    echo -e "${FAILMSG} ${ARG} ${FAILMSG2}"
    
    unset_env_vars
  }

  if [ -z "${DEBUG_POLY_OUTPUT}" ]; then
    export OUT_POLY_1="/dev/null"
    export OUT_POLY_2="/dev/null"
  else
    export OUT_POLY_1="/dev/tty"
    export OUT_POLY_2="/dev/tty"
  fi

  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  export POLY_NAME="PolyChordLite"
  export POLYCHORD_URL="https://github.com/PolyChord/PolyChordLite.git"
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------

  rm -rf $ROOTDIR/external_modules/code/$POLY_NAME

  cd $ROOTDIR/external_modules/code/

  $GIT clone $POLYCHORD_URL $POLY_NAME > ${OUT_POLY_1} 2> ${OUT_POLY_2}
  if [ $? -ne 0 ]; then
    fail "GIT CLONE"
    return 1
  fi

  cd $ROOTDIR/external_modules/code/$POLY_NAME

  $GIT checkout $POLYCHORD_GIT_COMMIT > ${OUT_POLY_1} 2> ${OUT_POLY_2}
  if [ $? -ne 0 ]; then
    fail "GIT CHECKOUT"
    return 1
  fi

  cp $ROOTDIR/../cocoa_installation_libraries/polychord_changes/Makefile.patch .
  if [ $? -ne 0 ]; then
    fail "CP FILE PATCH"
    return 1
  fi

  patch -u Makefile -i Makefile.patch > ${OUT_POLY_1} 2> ${OUT_POLY_2}
  if [ $? -ne 0 ]; then
    fail "PATCH FILE"
    return 1
  fi

  cp $ROOTDIR/../cocoa_installation_libraries/polychord_changes/setup.patch .
  if [ $? -ne 0 ]; then
    fail "CP FILE PATCH"
    return 1
  fi

  patch -u setup.py -i setup.patch > ${OUT_POLY_1} 2> ${OUT_POLY_2}
  if [ $? -ne 0 ]; then
    fail "PATCH FILE"
    return 1
  fi

  unset_env_vars
  echo -e '\033[1;34m''\t\e[4mSETUP POLYCHORD DONE''\033[0m'
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------