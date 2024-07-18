#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_HYREC_CODE}" ]; then
  
  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi

  # parenthesis = run in a subshell  
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

  unset_env_vars () {
    unset -v ECODEF FOLDER PACKDIR PRINTNAME
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

  FOLDER="${HYREC_NAME:-"hyrec2"}"

  PACKDIR="${ECODEF:?}/${FOLDER:?}"

  # Name to be printed on this shell script messages
  PRINTNAME="HYREC2 RECOMBINATION CODE"

  ptop "COMPILING ${PRINTNAME:?}" || return 1

  # ---------------------------------------------------------------------------
  # cleaning any previous compilation
  
  rm -f  "${ROOTDIR:?}/.local/lib/libhyrec.a"
  rm -f  "${ROOTDIR:?}/.local/lib/libhyrec.so"
  rm -rf "${ROOTDIR:?}/.local/include/hyrec2"
  
  # ---------------------------------------------------------------------------
  
  cdfolder "${PACKDIR}" || return 1
  
  # ---------------------------------------------------------------------------
  # create .local/include/hyrec2 where headers will be located
  
  mkdir "${ROOTDIR:?}/.local/include/hyrec2"  \
    >${OUT1:?} 2>${OUT2:?}  || { error "${EC20:?}"; return 1; }

  CC="${C_COMPILER:?}" make install \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC10:?}"; return 1; }
  
  pbottom "COMPILING ${PRINTNAME:?}" || return 1

  # ---------------------------------------------------------------------------

  unset_all || return 1

fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------