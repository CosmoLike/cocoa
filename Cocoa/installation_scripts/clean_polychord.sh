#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_POLYCHORD_COMPILATION}" ]; then
  
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
  
  if [ -z "${PYTHON_VERSION}" ]; then
    pfail "PYTHON_VERSION"; cdroot; return 1;
  fi
    
  unset_env_vars_clean_poly () {
    unset OUT1
    unset OUT2
    unset pfail
    unset PLIB
    unset PACKDIR
    unset unset_env_vars_clean_poly
    cdroot || return 1;
  }
  
  fail_cl_poly () {
    local MSG="\033[0;31m\t\t (clean_polychord.sh) we cannot run \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG} ${1:-"empty arg"} ${MSG2}"
    unset fail_cl_poly
    unset_env_vars_clean_poly
    return 1
  }
  
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { fail_cl_poly "CD FOLDER: ${1}"; return 1; }
  }
  
  if [ -z "${COCOA_OUTPUT_VERBOSE}" ]; then
    export OUT1="/dev/null"; export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"; export OUT2="/dev/tty"
  fi
  
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  
  ptop 'CLEANING POLYCHORD'

  export PLIB="${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages"
  export PACKDIR="${ROOTDIR:?}/external_modules/code/${POLY_NAME:-"PolyChordLite"}"

  cdfolder "${PACKDIR}" || return 1

  make clean >${OUT1:?} 2>${OUT2:?} || { fail_cl_poly "MAKE CLEAN"; return 1; }
  
  rm -rf "${PLIB:?}"/pypolychord-*
  rm -rf "${ROOTDIR:?}/external_modules/code/${POLYF:?}/lib/*.a"
  rm -rf "${ROOTDIR:?}/external_modules/code/${POLYF:?}/lib/*.so"

  unset_env_vars_clean_poly || return 1
  
  pbottom 'CLEANING POLYCHORD'

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------