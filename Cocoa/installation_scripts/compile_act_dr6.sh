#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_ACTDR6_CODE}" ]; then
  
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

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------  
  # ---------------------------------------------------------------------------- 

  unset_env_vars || return 1

  # E = EXTERNAL, CODE, F=FODLER
  ECODEF="${ROOTDIR:?}/external_modules/code"

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------  

  ptop "COMPILING ACT-DR6 (CMBONLY)" || return 1

  PACKDIR="${ECODEF:?}/${ACTDR6_CMBONLY_NAME:-"act_dr6_cmbonly"}"

  # ---------------------------------------------------------------------------- 
  # cleaning any previous compilation
  rm -rf "${PACKDIR:?}/build/"
  rm -rf "${PACKDIR:?}/syslibrary.egg-info/"
  PLIB="${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages"
  rm -rf  "${PLIB:?}/${ACTDR6_CMBONLY_NAME:-"act_dr6_cmbonly"}"
  rm -rf  "${PLIB:?}/${ACTDR6_CMBONLY_NAME:-"act_dr6_cmbonly"}"-*
  # ----------------------------------------------------------------------------

  (export LD_LIBRARY_PATH=${CONDA_PREFIX:?}/lib:$LD_LIBRARY_PATH && \
   export LD_LIBRARY_PATH=${ROOTDIR:?}/.local/lib:$LD_LIBRARY_PATH && \
   env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" ${PIP3:?} install "${PACKDIR:?}" \
   --prefix="${ROOTDIR:?}/.local" --no-index --no-deps --no-build-isolation \
   >${OUT1:?} 2>${OUT2:?} || { error "${EC13:?}"; return 1; })
    
  cdfolder "${ROOTDIR}" || return 1

  pbottom "COMPILING ACT-DR6 (CMBONLY)" || return 1  

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------  

  ptop "COMPILING ACT-DR6 (MFLIKE)" || return 1

  PACKDIR="${ECODEF:?}/${ACTDR6_MFLIKE_NAME:-"act_dr6_mflike"}"

  # ---------------------------------------------------------------------------- 
  # cleaning any previous compilation
  rm -rf "${PACKDIR:?}/build/"
  rm -rf "${PACKDIR:?}/syslibrary.egg-info/"
  PLIB="${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages"
  rm -rf  "${PLIB:?}/${ACTDR6_MFLIKE_NAME:-"act_dr6_mflike"}"
  rm -rf  "${PLIB:?}/${ACTDR6_MFLIKE_NAME:-"act_dr6_mflike"}"-*
  # ----------------------------------------------------------------------------

  (export LD_LIBRARY_PATH=${CONDA_PREFIX:?}/lib:$LD_LIBRARY_PATH && \
   export LD_LIBRARY_PATH=${ROOTDIR:?}/.local/lib:$LD_LIBRARY_PATH && \
   env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" ${PIP3:?} install "${PACKDIR:?}" \
   --prefix="${ROOTDIR:?}/.local" --no-index --no-deps --no-build-isolation \
   >${OUT1:?} 2>${OUT2:?} || { error "${EC13:?}"; return 1; })
    
  cdfolder "${ROOTDIR}" || return 1
  
  pbottom "COMPILING ACT-DR6 (MFLIKE)" || return 1

  # ---------------------------------------------------------------------------  
  # ---------------------------------------------------------------------------

  cdfolder "${ROOTDIR}" || return 1

  unset_all || return 1
  
fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------