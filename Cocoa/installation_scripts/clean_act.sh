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
  unset_env_vars_clean_act () {
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset pfail
    unset unset_env_vars_clean_act
  }
  fail_clean_act () {
    export MSG="\033[0;31m (clean_act.sh) WE CANNOT RUN \e[3m"
    export MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"
    unset_env_vars_clean_act
    unset MSG
    unset MSG2
    unset fail_clean_act
  }
  if [ -z "${DEBUG_ACT_OUTPUT}" ]; then
    export OUT1="/dev/null"
    export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
  fi

  # --------------------------------------------------------------------------- 
  ptop 'CLEANING ACT'

  cd $ROOTDIR/external_modules/code/pyactlike/

  rm -rf ./build/
  rm -rf ./pyactlike.egg-info/

  $PYTHON3 setup.py clean > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_clean_act "PYTHON SETUP CLEAN"
    return 1
  fi

  unset_env_vars_clean_act
  pbottom 'CLEANING ACT'
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------