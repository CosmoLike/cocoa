#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_CAMB_COMPILATION}" ]; then
  # ----------------------------------------------------------------------------
  source $ROOTDIR/installation_scripts/clean_camb.sh
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
  if [ -z "${CAMB_NAME}" ]; then
    pfail 'CAMB_NAME'
    cd $ROOTDIR
    return 1
  fi
  unset_env_vars_comp_camb () {
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset CAMB_NAME
    unset pfail
    unset unset_env_vars_comp_camb
  }
  fail_comp_camb () {
    export MSG="\033[0;31m (compile_camb.sh) WE CANNOT RUN \e[3m"
    export MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"
    unset_env_vars_comp_camb
    unset MSG
    unset MSG2
    unset fail_comp_camb
  }
  if [ -z "${DEBUG_CAMB_OUTPUT}" ]; then
    export OUT1="/dev/null"
    export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
  fi
  
  # --------------------------------------------------------------------------- 
  ptop 'COMPILING CAMB'

  cd $ROOTDIR/external_modules/code/$CAMB_NAME/
  if [ $? -ne 0 ]; then
    fail_comp_camb "CD CAMB FOLDER"
    return 1
  fi

  COMPILER=$FORTRAN_COMPILER F90C=$FORTRAN_COMPILER $PYTHON3 setup.py build > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_comp_camb "PYTHON3 SETUP.PY BUILD"
    return 1
  fi

  unset_env_vars_comp_camb
  pbottom 'COMPILING CAMB'
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------