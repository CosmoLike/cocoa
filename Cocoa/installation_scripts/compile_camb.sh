#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_CAMB_COMPILATION}" ]; then
  
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1
  fi

  # ----------------------------------------------------------------------------
  # Clean any previous compilation
  source "${ROOTDIR:?}/installation_scripts/clean_camb.sh" || return 1
  # ----------------------------------------------------------------------------
  
  source "${ROOTDIR:?}/installation_scripts/.check_flags.sh" || return 1;

  unset_env_vars () {
    unset -v PACKDIR
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
    fail_script_msg "compile_camb.sh" "${1}"
    unset_all || return 1
  }
  
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER ${1}"; return 1; }
  }
  
  # --------------------------------------------------------------------------- 
  # --------------------------------------------------------------------------- 
  # ---------------------------------------------------------------------------

  ptop 'COMPILING CAMB' || return 1

  unset_env_vars || return 1

  PACKDIR="${ROOTDIR:?}/external_modules/code/${CAMB_NAME:-"CAMB"}"

  cdfolder "${PACKDIR}" || return 1

  COMPILER="${FORTRAN_COMPILER:?}" F90C="${FORTRAN_COMPILER:?}" \
    ${PYTHON3:?} setup.py build \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC4:?}"; return 1; }

  unset_all || return 1
  
  pbottom 'COMPILING CAMB' || return 1
  
  # --------------------------------------------------------------------------- 
  # --------------------------------------------------------------------------- 
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------