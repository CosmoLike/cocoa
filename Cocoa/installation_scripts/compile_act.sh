#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_ACT_COMPILATION}" ]; then
  
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1
  fi
    
  ( source "${ROOTDIR:?}/installation_scripts/.check_flags.sh" ) || return 1;
 
  unset_env_vars () {
    unset -v ECODEF ACTF PACKDIR
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
    fail_script_msg "$(basename ${BASH_SOURCE[0]})" "${1}"
    unset_all || return 1
  }
  
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER ${1}"; return 1; }
  }

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------  
  # ---------------------------------------------------------------------------  
  
  ptop 'COMPILING ACT' || return 1

  unset_env_vars || return 1

  # E = EXTERNAL, CODE, F=FODLER
  ECODEF="${ROOTDIR:?}/external_modules/code"

  ACTF=${ACT_NAME:-"pyactlike"}

  PACKDIR="${ECODEF:?}/${ACTF:?}"

  cdfolder "${PACKDIR}" || return 1

  # ---------------------------------------------------------------------------
  # cleaning any previous compilation

  rm -rf "${PACKDIR:?}/build/"
  rm -rf "${PACKDIR:?}/pyactlike.egg-info/"

  "${PYTHON3:?}" setup.py clean \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC1:?}"; return 1; }
  
  # ---------------------------------------------------------------------------  
 
  ${PIP3:?} install . --prefix="${ROOTDIR:?}/.local" \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC3:?}"; return 1; }

  unset_all || return 1
  
  pbottom 'COMPILING ACT' || return 1
  
fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------