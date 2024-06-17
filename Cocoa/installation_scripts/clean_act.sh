#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_ACT_COMPILATION}" ]; then
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
  if [ -z "${DEBUG_ACT_OUTPUT}" ]; then
    export OUT1="/dev/null"
    export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
  fi

  # --------------------------------------------------------------------------- 
  ptop2 'CLEANING ACT'

  cd $ROOTDIR/external_modules/code/pyactlike/

  rm -rf ./build/
  rm -rf ./pyactlike.egg-info/

  $PYTHON3 setup.py clean > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "PYTHON SETUP CLEAN"
    return 1
  fi

  unset_env_vars
  pbottom2 'CLEANING ACT'
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------