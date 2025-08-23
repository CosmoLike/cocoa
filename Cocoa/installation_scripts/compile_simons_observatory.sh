#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_SIMONS_OBSERVATORY_LIKELIHOOD_CODE}" ]; then
  
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
 
  unset_env_vars || return 1

  # E = EXTERNAL, CODE, F=FODLER
  ECODEF="${ROOTDIR:?}/external_modules/code"

  PLIB="${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages"
  
  # ----------------------------------------------------------------------------  
  # ----------------------------------------------------------------------------
  
  ptop "COMPILING SIMONS OBSERVATORY SYSLIBRARY" || return 1

  PACKDIR="${ECODEF:?}/${SO_SYSLIB_NAME:-"syslibrary"}"

  # ---------------------------------------------------------------------------- 
  # cleaning any previous compilation
  rm -rf "${PACKDIR:?}/build/"
  rm -rf "${PACKDIR:?}/syslibrary.egg-info/"  
  rm -rf  "${PLIB:?}/syslibrary"
  rm -rf  "${PLIB:?}/syslibrary"-*
  # ----------------------------------------------------------------------------

  cdfolder "${PACKDIR}" || return 1

  #prevent all compile_XXX.sh from using the internet (run @compute nodes)
  #FROM: https://github.com/pypa/pip/issues/12050
  #That is why we use --no-dependencies --no-index --no-build-isolation
  (env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" \
    ${PIP3:?} install . \
      --prefix="${ROOTDIR:?}/.local" \
      --no-index \
      --no-deps \
      --no-build-isolation \
  ) >>${OUT1:?} 2>>${OUT2:?} || { error "${EC13:?}"; return 1; }
    
  cdfolder "${ROOTDIR}" || return 1

  pbottom "COMPILING SIMONS OBSERVATORY SYSLIBRARY" || return 1

  # ----------------------------------------------------------------------------  
  # ----------------------------------------------------------------------------

  ptop "COMPILING SIMONS OBSERVATORY MKLIKE" || return 1

  PACKDIR="${ECODEF:?}/mflike"

  # ----------------------------------------------------------------------------
  # cleaning any previous compilation
  rm -rf "${PACKDIR:?}/build/"
  rm -rf "${PACKDIR:?}/mflike.egg-info/"
  rm -rf  "${PLIB:?}/mflike"
  rm -rf  "${PLIB:?}/mflike"-*
  # ----------------------------------------------------------------------------

  (env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" \
    ${PIP3:?} install ${PACKDIR:?} \
      --prefix="${ROOTDIR:?}/.local" \
      --no-index \
      --no-dependencies \
      --no-build-isolation
  ) >>${OUT1:?} 2>>${OUT2:?} || { error "${EC13:?}"; return 1; }

  cdfolder "${ROOTDIR}" || return 1

  pbottom "COMPILING SIMONS OBSERVATORY MKLIKE" || return 1

  cdfolder "${ROOTDIR}" || return 1

  unset_all || return 1  
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------