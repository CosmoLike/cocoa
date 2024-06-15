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
  if [ -z "${GIT}" ]; then
    pfail 'GIT'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${POLYCHORD_GIT_COMMIT}" ]; then
    pfail 'POLYCHORD_GIT_COMMIT'
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
    unset POLY_URL
    unset POLY_CHANGES
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
  ptop 'INSTALLING POLYCHORD'

  export POLY_URL="https://github.com/PolyChord/PolyChordLite.git"
  export POLY_CHANGES="${ROOTDIR}/../cocoa_installation_libraries/polychord_changes"

  rm -rf $ROOTDIR/external_modules/code/$POLY_NAME

  cd $ROOTDIR/external_modules/code/

  $GIT clone $POLY_URL $POLY_NAME > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "GIT CLONE"
    return 1
  fi

  cd $ROOTDIR/external_modules/code/$POLY_NAME

  $GIT checkout $POLYCHORD_GIT_COMMIT > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "GIT CHECKOUT"
    return 1
  fi

  cp $POLY_CHANGES/Makefile.patch .
  if [ $? -ne 0 ]; then
    fail "CP FILE PATCH (MAKEFILE)"
    return 1
  fi

  patch -u Makefile -i Makefile.patch > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "PATCH FILE (MAKEFILE)"
    return 1
  fi

  cp $POLY_CHANGES/setup.patch .
  if [ $? -ne 0 ]; then
    fail "CP FILE PATCH (SETUP)"
    return 1
  fi

  patch -u setup.py -i setup.patch > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "PATCH FILE (SETUP)"
    return 1
  fi

  unset_env_vars
  pbottom 'INSTALLING POLYCHORD'
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------