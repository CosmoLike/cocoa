#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_ACT_COMPILATION}" ]; then
  
  pfail() {
    echo -e \
    "\033[0;31m\t\t ERROR ENV VARIABLE ${1:-"empty arg"} NOT DEFINED \033[0m"
    unset pfail
  }
  
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1
  fi
  
  cdroot() {
    cd "${ROOTDIR:?}" 2>"/dev/null" || { echo -e \
      "\033[0;31m\t\t CD ROOTDIR (${ROOTDIR}) FAILED \033[0m"; return 1; }
    unset cdroot
  }
  
  if [ -z "${PYTHON3}" ]; then
    pfail "PYTHON3"; cdroot; return 1;
  fi
    
  unset_env_vars_clean_act () {
    unset pfail
    unset OUT1
    unset OUT2
    unset cdfolder
    unset PACKDIR
    unset unset_env_vars_clean_act
    cdroot || return 1;
  }
  
  fail_clact () {
    local MSG="\033[0;31m\t\t (clean_act.sh) we cannot run \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG}${1:-"empty arg"}${MSG2}"
    unset fail_clact
    unset_env_vars_clean_act
  }
  
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { fail_clact "CD FOLDER: ${1}"; return 1; }
  }
  
  if [ -z "${COCOA_OUTPUT_VERBOSE}" ]; then
    export OUT1="/dev/null"; export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"; export OUT2="/dev/tty"
  fi

  # ---------------------------------------------------------------------------
  # --------------------------------------------------------------------------- 
  
  ptop 'CLEANING ACT'

  export PACKDIR="${ROOTDIR:?}/external_modules/code/${ACT_NAME:-"pyactlike"}";

  rm -rf "${PACKDIR:?}/build/"
  rm -rf "${PACKDIR:?}/pyactlike.egg-info/"

  cdfolder "${PACKDIR}" || return 1

  "${PYTHON3:?}" setup.py clean >${OUT1:?} 2>${OUT2:?} ||
    { fail_clact "PYTHON SETUP CLEAN"; return 1; }

  unset_env_vars_clean_act || return 1
  
  pbottom 'CLEANING ACT'
  
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------

fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------