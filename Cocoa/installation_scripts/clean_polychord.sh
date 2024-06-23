#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_POLYCHORD_COMPILATION}" ]; then
    
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1
  fi
  
  source "${ROOTDIR:?}/installation_scripts/.check_flags.sh" || return 1;
    
  unset_env_vars () {
    unset -v PLIB PACKDIR
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
    fail_script_msg "clean_polychord.sh" "${1}"
    unset_all || return 1
  }
  
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
  }
    
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------

  ptop 'CLEANING POLYCHORD' || return 1

  unset_env_vars || return 1

  PLIB="${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages"
  
  PACKDIR="${ROOTDIR:?}/external_modules/code/${POLY_NAME:-"PolyChordLite"}"

  cdfolder "${PACKDIR}" || return 1

  make clean >${OUT1:?} 2>${OUT2:?} || { error "${EC2:?}"; return 1; }
  
  rm -rf "${PLIB:?}"/pypolychord-*
  rm -rf "${PACKDIR:?}/lib/*.a"
  rm -rf "${PACKDIR:?}/lib/*.so"

  unset_all || return 1
  
  pbottom 'CLEANING POLYCHORD' || return 1
  
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------