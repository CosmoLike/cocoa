#!/bin/bash

if [ -z "${IGNORE_PLANCK_CMB_DATA}" ]; then

  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi

  # parenthesis = run in a subshell 
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;
    
  unset_env_vars () {
    unset -v EDATAF PRINTNAME
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
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
  }

  # --------------------------------------------------------------------------- 
  # --------------------------------------------------------------------------- 
  # ---------------------------------------------------------------------------

  unset_env_vars || return 1

  # E = EXTERNAL, DATA, F=FODLER
  EDATAF="${ROOTDIR:?}/external_modules/data"

  # Name to be printed on this shell script messages
  PRINTNAME="PLANCK-2018 (PLC 3.0)"

  # ---------------------------------------------------------------------------
 
  ptop "SETUP/UNXV ${PRINTNAME:?} DATA" || return 1

  # ---------------------------------------------------------------------------
  # note: in case this script is run twice
  if [ -n "${OVERWRITE_EXISTING_PLANCK_CMB_DATA}" ]; then
 
    rm -rf "${EDATAF:?}/planck/plc_3.0/low_l"
    rm -rf "${EDATAF:?}/planck/plc_3.0/lensing"
    rm -rf "${EDATAF:?}/planck/plc_3.0/hi_l/plik"
    rm -rf "${EDATAF:?}/planck/plc_3.0/hi_l/plik_lite"

  fi

  if [[ ! -d "${EDATAF:?}/planck/plc_3.0/low_l" && \
        ! -d "${EDATAF:?}/planck/plc_3.0/lensing" && \
        ! -d "${EDATAF:?}/planck/plc_3.0/hi_l/plik" && \
        ! -d "${EDATAF:?}/planck/plc_3.0/hi_l/plik_lite" ]]; then

    cdfolder "${EDATAF:?}/planck/plc_3.0" || return 1
    
    tar xf lensing.xz \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC25:?}"; return 1; }

    tar xf low_l.xz \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC25:?}"; return 1; }

    cdfolder "${EDATAF:?}/planck/plc_3.0/hi_l" || return 1
    
    tar xf plik.xz \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC25:?}"; return 1; }

    tar xf plik_lite.xz \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC25:?}"; return 1; }
    
  fi

  cdfolder "${ROOTDIR}" || return 1

  pbottom "SETUP/UNXV ${PRINTNAME:?} DATA" || return 1

  unset_all || return 1;

fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------