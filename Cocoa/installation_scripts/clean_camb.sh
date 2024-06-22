#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_CAMB_COMPILATION}" ]; then
    
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1
  fi
  
  source "${ROOTDIR:?}/installation_scripts/.check_flags.sh" || return 1;
  
  unset_env_vars () {
    unset -v PACKDIR
    cdroot || return 1;
  }

  error () {
    fail_script_msg "clean_camb.sh" "${1}"
    unset -f error
    unset_env_vars || return 1;
  }

  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
  }
  
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  ptop 'CLEANING CAMB' || return 1

  unset_env_vars || return 1

  PACKDIR="${ROOTDIR:?}/external_modules/code/${CAMB_NAME:-"CAMB"}"

  rm -rf "${PACKDIR:?}/build/"
  rm -rf "${PACKDIR:?}/camb/__pycache__/"
  rm -f  "${PACKDIR:?}/camb/camblib.so"
  rm -rf "${PACKDIR:?}/forutils/Releaselib/"

  cdfolder "${PACKDIR}" || return 1

  "${PYTHON3:?}" setup.py clean >${OUT1:?} 2>${OUT2:?} || 
    { error "${EC1:?}"; return 1; }
  
  unset_env_vars || return 1
  
  pbottom 'CLEANING CAMB' || return 1
  
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------