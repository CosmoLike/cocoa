#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_ACT_COMPILATION}" ]; then
  # ---------------------------------------------------------------------------
  source $ROOTDIR/installation_scripts/clean_act.sh
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
  unset_env_vars_compile_class () {
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset pfail
    unset unset_env_vars_compile_class
  }
  fail_comp_class () {
    export MSG="\033[0;31m (compile_act.sh) WE CANNOT RUN \e[3m"
    export MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"
    unset_env_vars_compile_class
    unset MSG
    unset MSG2
    unset fail_comp_class
  }
  if [ -z "${DEBUG_ACT_OUTPUT}" ]; then
    export OUT1="/dev/null"
    export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
  fi

  # ---------------------------------------------------------------------------  
  ptop 'COMPILING ACT'

  cd $ROOTDIR/external_modules/code/pyactlike/
 
  $PIP3 install . --prefix=$ROOTDIR/.local > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_comp_class "PIP3 INSTALL ."
    return 1
  fi

  unset_env_vars_compile_class
  pbottom 'COMPILING ACT'
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------