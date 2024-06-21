#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_CAMB_COMPILATION}" ]; then
  
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
    pfail "PYTHON3"; cdroot; return 1
  fi
  
  unset_env_vars_clean_camb () {
    unset OUT1
    unset OUT2
    unset pfail
    unset PACKDIR
    unset unset_env_vars_clean_camb
    cdroot || return 1;
  }

  fail_clcb () {
    local MSG="\033[0;31m\t\t (clean_camb.sh) WE CANNOT RUN \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG} ${1:-"empty arg"} ${MSG2}"
    unset fail_clcb
    unset_env_vars_clean_camb
  }

  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { fail_clcb "CD FOLDER: ${1}"; return 1; }
  }
  
  if [ -z "${COCOA_OUTPUT_VERBOSE}" ]; then
    export OUT1="/dev/null"; export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"; export OUT2="/dev/tty"
  fi

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  ptop 'CLEANING CAMB'

  export PACKDIR="${ROOTDIR:?}/external_modules/code/${CAMB_NAME:-"CAMB"}"

  rm -rf "${PACKDIR:?}/build/"
  rm -rf "${PACKDIR:?}/camb/__pycache__/"
  rm -f  "${PACKDIR:?}/camb/camblib.so"
  rm -rf "${PACKDIR:?}/forutils/Releaselib/"

  cdfolder "${PACKDIR}" || return 1

  "${PYTHON3:?}" setup.py clean >${OUT1:?} 2>${OUT2:?} ||
    { fail_clcb "PYTHON SETUP CLEAN"; return 1; }

  unset_env_vars_clean_camb || return 1
  
  pbottom 'CLEANING CAMB'
  
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------

fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------