#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${SKIP_DECOMM_BICEP}" ]; then

  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1;
  fi

  unset_env_vars () {
    unset -v EDATAF FOLDER PACKDIR FILE
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
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
  }

  # --------------------------------------------------------------------------- 
  # --------------------------------------------------------------------------- 
  # ---------------------------------------------------------------------------

  unset_env_vars || return 1

  # E = EXTERNAL, DATA, F=FODLER
  EDATAF="${ROOTDIR:?}/external_modules/data"

  FOLDER="bicep_keck_2015"

  # PACK = PACKAGE, DIR = DIRECTORY
  PACKDIR="${EDATAF:?}/${FOLDER:?}"

  FILE="bicep_keck_2015.xz"

  # ---------------------------------------------------------------------------
  
  ptop 'DECOMPRESSING BICEP 2015 DATA' || return 1

  # ---------------------------------------------------------------------------
  # note: in case script run >1x w/ previous run stoped prematurely b/c error
  
  rm -rf "${PACKDIR:?}"
  
  # ---------------------------------------------------------------------------

  cdfolder "${EDATAF:?}" || return 1

  tar xf "${FILE:?}" >${OUT1:?} 2>${OUT2:?} || { error "${EC25:?}"; return 1; }

  # ---------------------------------------------------------------------------

  cdfolder "${ROOTDIR}" || return 1

  unset_all || return 1
  
  pbottom 'DECOMPRESSING BICEP 2015 DATA' || return 1

fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------