#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${SKIP_DECOMM_BAO}" ]; then

  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1;
  fi

  unset_env_vars () {
    unset -v EDATAF DATAF PACKDIR FILE
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
    fail_script_msg "unxv_bao.sh" "${1}"
    unset_all || return 1
  }
  
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER ${1}"; return 1; }
  }

  # --------------------------------------------------------------------------- 
  # --------------------------------------------------------------------------- 
  # ---------------------------------------------------------------------------

  ptop 'DECOMPRESSING BAO DATA' || return 1

  unset_env_vars || return 1

  # E = EXTERNAL, DATA, F=FODLER
  EDATAF="${ROOTDIR:?}/external_modules/data"

  DATAF="bao_data"

  # PACK = PACKAGE, DIR = DIRECTORY
  PACKDIR="${EDATAF:?}/${DATAF:?}"

  FILE="bao_data.xz"

  # ---------------------------------------------------------------------------
  # in case this script is called twice
  
  rm -rf "${PACKDIR:?}"

  # ---------------------------------------------------------------------------

  cdfolder "${EDATAF:?}" || return 1

  tar xf ${FILE:?} >${OUT1:?} 2>${OUT2:?} || { error "${EC25:?}"; return 1; }

  unset_all || return 1
  
  pbottom 'DECOMPRESSING BAO DATA' || return 1

fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------