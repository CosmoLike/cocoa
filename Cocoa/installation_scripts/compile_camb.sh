#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_CAMB_COMPILATION}" ]; then
  
  if [ -z "${ROOTDIR}" ]; then
    echo -e "\033[0;31m\t\t ERROR ENV VARIABLE ${ROOTDIR} IS NOT DEFINED \033[0m"
    return 1
  fi

  # ----------------------------------------------------------------------------
  # Clean any previous compilation
  source "${ROOTDIR:?}/installation_scripts/clean_camb.sh"
  # ----------------------------------------------------------------------------
  
  pfail() {
    echo -e \
    "\033[0;31m\t\t ERROR ENV VARIABLE ${1:-"empty arg"} NOT DEFINED \033[0m"
    unset pfail
  }

  cdroot() {
    cd "${ROOTDIR:?}" 2>"/dev/null" || { echo -e \
      "\033[0;31m\t\t CD ROOTDIR (${ROOTDIR}) FAILED \033[0m"; return 1; }
    unset cdroot
  }
  
  if [ -z "${CXX_COMPILER}" ]; then
    pfail 'CXX_COMPILER'; cdroot; return 1;
  fi
  
  if [ -z "${C_COMPILER}" ]; then
    pfail 'C_COMPILER'; cdroot; return 1;
  fi
  
  if [ -z "${PIP3}" ]; then
    pfail 'PIP3'; cdroot; return 1;
  fi
  
  if [ -z "${PYTHON3}" ]; then
    pfail "PYTHON3"; cdroot; return 1;
  fi
    
  if [ -z "${CAMB_NAME}" ]; then
    pfail 'CAMB_NAME'; cdroot; return 1;
  fi
  
  unset_env_vars_comp_camb () {
    unset OUT1
    unset OUT2
    unset pfail
    unset unset_env_vars_comp_camb
    cdroot || return 1;
  }
  
  fail_comp_camb () {
    local MSG="\033[0;31m\t\t (compile_camb.sh) we cannot run \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG} ${1:-"empty arg"} ${MSG2}"
    unset fail_comp_camb
    unset_env_vars_comp_camb
  }
  
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { fail_comp_camb "CD FOLDER: ${1}"; return 1; }
  }
  
  if [ -z "${DEBUG_CAMB_OUTPUT}" ]; then
    export OUT1="/dev/null"; export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"; export OUT2="/dev/tty"
  fi

  # --------------------------------------------------------------------------- 
  # --------------------------------------------------------------------------- 
  
  ptop 'COMPILING CAMB'

  cdfolder "${ROOTDIR:?}/external_modules/code/${CAMB_NAME:?}/" || return 1

  COMPILER=${FORTRAN_COMPILER:?} F90C=${FORTRAN_COMPILER:?} ${PYTHON3:?} \
    setup.py build >${OUT1:?} 2>${OUT2:?} || 
    { fail_comp_camb "PYTHON3 SETUP.PY BUILD"; return 1; }

  unset_env_vars_comp_camb || return 1
  
  pbottom 'COMPILING CAMB'
  
  # --------------------------------------------------------------------------- 
  # --------------------------------------------------------------------------- 
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------