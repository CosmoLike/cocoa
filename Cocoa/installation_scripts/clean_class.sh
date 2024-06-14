#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_CLASS_COMPILATION}" ]; then
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
  if [ -z "${CLASS_NAME}" ]; then
    pfail 'CAMB_NAME'
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
  if [ -z "${DEBUG_CLASS_OUTPUT}" ]; then
    export OUT1="/dev/null"
    export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
  fi

  # ---------------------------------------------------------------------------
  ptop 'CLEANING CLASS'

  cd $ROOTDIR/external_modules/code/$CLASS_NAME/

  rm $ROOTDIR/.local/lib/python$PYTHON_VERSION/site-packages/classy*

  #Workaround Cocoa .gitignore entry on /include
  rm -rf ./include > ${OUT1} 2> ${OUT2}

  make clean > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "MAKE CLEAN"
    return 1
  fi

  rm -f class

  cd ./python

  $PYTHON3 setup.py clean > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "PYTHON SETUP CLEAN"
    return 1
  fi

  rm -rf ./build/
  rm -rf ./classy.egg-info

  unset_env_vars
  pbottom 'CLEANING CLASS'
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------