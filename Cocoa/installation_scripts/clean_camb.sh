#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_CAMB_COMPILATION}" ]; then
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
  if [ -z "${CAMB_NAME}" ]; then
    pfail 'CAMB_NAME'
    cd $ROOTDIR
    return 1
  fi
  unset_env_vars_clean_camb () {
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset pfail
    unset unset_env_vars_clean_camb
  }
  fail_clean_camb () {
    export MSG="\033[0;31m (clean_camb.sh) WE CANNOT RUN \e[3m"
    export MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"
    unset_env_vars_clean_camb
    unset MSG
    unset MSG2
    unset fail_clean_camb
  }
  if [ -z "${DEBUG_CAMB_OUTPUT}" ]; then
    export OUT1="/dev/null"
    export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
  fi

  # ---------------------------------------------------------------------------
  ptop 'CLEANING CAMB'

  cd $ROOTDIR/external_modules/code/$CAMB_NAME/

  rm -rf ./build/

  rm -rf ./camb/__pycache__/

  rm -f ./camb/camblib.so

  rm -rf ./forutils/Releaselib/

  $PYTHON3 setup.py clean > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_clean_camb "PYTHON SETUP CLEAN"
    return 1
  fi

  unset_env_vars_clean_camb
  pbottom 'CLEANING CAMB'
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------