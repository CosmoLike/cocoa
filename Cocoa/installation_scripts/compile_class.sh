#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_CLASS_COMPILATION}" ]; then
  # ---------------------------------------------------------------------------
  source $ROOTDIR/installation_scripts/clean_class.sh
  # ---------------------------------------------------------------------------
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
  if [ -z "${CLASS_NAME}" ]; then
    pfail 'CLASS_NAME'
    cd $ROOTDIR
    return 1
  fi
  unset_env_vars_comp_class () {
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset pfail
    unset unset_env_vars_comp_class
  }
  fail_comp_class () {
    export MSG="\033[0;31m (compile_class.sh) WE CANNOT RUN \e[3m"
    export MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"
    unset_env_vars_comp_class
    unset MSG
    unset MSG2
    unset fail_comp_class
  }
  if [ -z "${DEBUG_CLASS_OUTPUT}" ]; then
    export OUT1="/dev/null"
    export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
  fi

  # ---------------------------------------------------------------------------
  ptop 'COMPILING CLASS'

  cd $ROOTDIR/external_modules/code/$CLASS_NAME/
  if [ $? -ne 0 ]; then
    fail_comp_class "CD CLASS FOLDER"
    return 1
  fi

  # ---------------------------------------------------------------------------
  # historical: workaround Cocoa .gitignore entry on /include 
  # also good because class compilation changes files on /include
  # in case this script is called twice
  # ---------------------------------------------------------------------------  
  rm -rf ./include

  # ---------------------------------------------------------------------------
  # historical: workaround Cocoa .gitignore entry on /include
  # also good because class compilation changes files on /include
  # ---------------------------------------------------------------------------
  cp -r  ./include2 ./include > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_comp_class "CP INCLUDE2 FOLDER"
    return 1
  fi

  CC=$C_COMPILER PYTHON=$PYTHON3 make all > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_comp_class "MAKE ALL"
    return 1
  fi
   
  cd ./python

  CC=$C_COMPILER $PYTHON3 setup.py build > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_comp_class "PYTHON3 SETUP.PY BUILD"
    return 1
  fi

  unset_env_vars_comp_class
  pbottom 'COMPILING CLASS'
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------