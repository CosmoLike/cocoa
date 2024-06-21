#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_PLANCK_COMPILATION}" ]; then
  
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
  
  unset_env_vars_clean_planck () { 
    unset OUT1
    unset OUT2
    unset pfail
    unset ECF
    unset unset_env_vars_clean_planck
    cdroot || return 1;
  }
  
  fail_clplk () {
    local MSG="\033[0;31m\t\t (clean_planck.sh) we cannot run \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG} ${1:-"empty arg"} ${MSG2}"
    unset fail_clplk
    unset_env_vars_clean_planck
  }
  
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { fail_clplk "CD FOLDER: ${1}"; return 1; }
  }
  
  if [ -z "${COCOA_OUTPUT_VERBOSE}" ]; then
    export OUT1="/dev/null"; export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"; export OUT2="/dev/tty"
  fi
  
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  
  ptop 'CLEANING PLANCK LIKELIHOOD' || return 1

  export ECF="external_modules/code/planck/code"
  
  if [ -z "${USE_SPT_CLIK_PLANCK}" ]; then
    cdfolder "${ROOTDIR:?}/${ECF:?}/plc_3.0/plc-3.1/" || return 1
  else
    cdfolder "${ROOTDIR:?}/${ECF:?}/spt_clik/" || return 1
  fi

  rm -f  "${ROOTDIR:?}/.local"/bin/clik*
  rm -f  "${ROOTDIR:?}/.local"/lib/libclik_f90.so
  rm -f  "${ROOTDIR:?}/.local"/lib/libclik.so
  rm -rf "${ROOTDIR:?}/.local"/lib/python/site-packages/clik
  rm -rf "${ROOTDIR:?}/.local"/share/clik
  rm -f  "${ROOTDIR:?}/.local"/include/clik*
  rm -f  .lock-waf_*

  "${PYTHON3:?}" waf distclean >${OUT1:?} 2>${OUT2:?} ||
    { fail_clplk "PYTHON WAF DISTCLEAN"; return 1; }

  unset_env_vars_clean_planck || return 1
  
  pbottom 'CLEANING PLANCK LIKELIHOOD' || return 1
  
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------

fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------