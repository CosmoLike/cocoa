#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_FGSPECTRA_CODE}" ]; then
  
  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi
    
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;
 
  unset_env_vars () {
    unset -v ECODEF FOLDER PACKDIR PRINTNAME PLIB
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
    fail_script_msg "$(basename "${BASH_SOURCE[0]}")" "${1}"
    unset_all || return 1
  }
  
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER ${1}"; return 1; }
  }

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------  
  # ---------------------------------------------------------------------------  
  
  unset_env_vars || return 1

  # ---------------------------------------------------------------------------

  # E = EXTERNAL, CODE, F=FODLER
  ECODEF="${ROOTDIR:?}/external_modules/code"

  FOLDER=${FGSPECTRA_NAME:-"fgspectra"}

  PACKDIR="${ECODEF:?}/${FOLDER:?}"

  # Name to be printed on this shell script messages
  PRINTNAME="FGSPECTRA"

  ptop "COMPILING ${PRINTNAME:?}" || return 1

  cdfolder "${PACKDIR}" || return 1

  # ---------------------------------------------------------------------------
  # cleaning any previous compilation
  rm -rf "${PACKDIR:?}/build/"
  rm -rf "${PACKDIR:?}/fgspectra.egg-info/"
  
  "${PYTHON3:?}" setup.py clean >${OUT1:?} 2>${OUT2:?} || { error "${EC1:?}"; return 1; }

  PLIB="${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages"
  rm -rf  "${PLIB:?}"/fgspectra
  rm -rf  "${PLIB:?}"/fgspectra-*  
  # ---------------------------------------------------------------------------  
 
  (export LD_LIBRARY_PATH=${CONDA_PREFIX:?}/lib:$LD_LIBRARY_PATH && \
   export LD_LIBRARY_PATH=${ROOTDIR:?}/.local/lib:$LD_LIBRARY_PATH && \
   ${PIP3:?} install . --prefix="${ROOTDIR:?}/.local" \
   >${OUT1:?} 2>${OUT2:?} || { error "${EC3:?}"; return 1; })

  pbottom "COMPILING ${PRINTNAME:?}" || return 1

  # ---------------------------------------------------------------------------

  unset_all || return 1
fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------