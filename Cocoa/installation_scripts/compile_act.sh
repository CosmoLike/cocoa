#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_ACT_COMPILATION}" ]; then
  
  if [ -z "${ROOTDIR}" ]; then
    echo -e "\033[0;31m\t\t ERROR ENV VARIABLE ROOTDIR NOT DEFINED \033[0m"
    return 1
  fi
  
  # ---------------------------------------------------------------------------
  # Clean any previous compilation
  source "${ROOTDIR:?}/installation_scripts/clean_act.sh" || return 1
  # ---------------------------------------------------------------------------
  
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
      
  unset_env_vars_comp_act () {
    unset OUT1
    unset OUT2
    unset pfail
    unset PACKDIR
    unset unset_env_vars_comp_act
    cdroot || return 1;
  }
  
  fail_comp_act () {
    local MSG="\033[0;31m\t\t (compile_act.sh) we cannot run \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"
    unset fail_comp_act
    unset_env_vars_comp_act
  }
  
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { fail_comp_act "CD FOLDER: ${1}"; return 1; }
  }
  
  if [ -z "${COCOA_OUTPUT_VERBOSE}" ]; then
    export OUT1="/dev/null"; export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"; export OUT2="/dev/tty"
  fi

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------  
  
  ptop 'COMPILING ACT'

  export PACKDIR="${ROOTDIR:?}/external_modules/code/${ACT_NAME:-"pyactlike"}"

  cdfolder "${PACKDIR}" || return 1
 
  ${PIP3:?} install . --prefix="${ROOTDIR:?}/.local" \
    >${OUT1:?} 2>${OUT2:?} #|| { fail_comp_act "PIP3 INSTALL ."; return 1; }

  #unset_env_vars_comp_act || return 1
  
  pbottom 'COMPILING ACT'
  
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------