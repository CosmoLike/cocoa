#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_PLANCK_COMPILATION}" ]; then
  
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1
  fi
  
  # parenthesis = run in a subshell
  ( source "${ROOTDIR:?}/installation_scripts/.check_flags.sh" ) || return 1;

  unset_env_vars () { 
    unset -v ECF
    cdroot || return 1;
  }

  unset_env_funcs () {
    unset -f cdfolder cpfolder error
    unset -f unset_env_funcs
    cdroot || return 1;
  }

  unset_all () {
    unset_env_vars
    unset_env_funcs
    unset -f unset_all
    cdroot || return 1;
  }
  
  error () {
    fail_script_msg "clean_planck.sh" "${1}"
    unset_all || return 1;
  }
  
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
  }
    
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  
  ptop 'CLEANING PLANCK LIKELIHOOD' || return 1

  ECF="external_modules/code/planck/code"
  
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

  "${PYTHON3:?}" waf distclean \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC18\?}"; return 1; }

  unset_all || return 1
  
  pbottom 'CLEANING PLANCK LIKELIHOOD' || return 1

fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------