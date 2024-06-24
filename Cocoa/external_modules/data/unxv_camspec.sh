#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${SKIP_DECOMM_CAMSPEC}" ]; then
  
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1;
  fi

  unset_env_vars () {
    unset -v EDATAF DATAF PACKDIR FILE URL_BASE URL
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
    fail_script_msg "unxv_camspec.sh" "${1}"
    unset_all || return 1
  }
  
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER ${1}"; return 1; }
  }

  # --------------------------------------------------------------------------- 
  # --------------------------------------------------------------------------- 
  # ---------------------------------------------------------------------------

  ptop "DECOMPRESSING CAMSPEC DATA" || return 1

  unset_env_vars || return 1

  # E = EXTERNAL, DATA, F=FODLER
  EDATAF="${ROOTDIR:?}/external_modules/data/planck/CamSpec"
  
  DATAF="CamSpec2021"

  # PACK = PACKAGE, DIR = DIRECTORY
  PACKDIR="${EDATAF:?}/${DATAF:?}"

  URL_BASE="https://github.com/CobayaSampler/planck_native_data/"

  FILE="CamSpec2021.zip"

  URL="${URL_BASE:?}/releases/download/v1/${FILE:?}"

  # ---------------------------------------------------------------------------
  # in case this script is called twice
  # ---------------------------------------------------------------------------
  rm -rf "${PACKDIR:?}"
  rm -rf "${EDATAF:?}/${FILE:?}"
  
  cdfolder "${EDATAF:?}" || return 1

  wget "${URL:?}" >${OUT1:?} 2>${OUT2:?} || { error "${EC24:?}"; return 1; }
  
  unzip "${FILE:?}" \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC26:?}"; return 1; }
  
  unset_all || return 1
  
  pbottom 'DECOMPRESSING CAMSPEC DATA' || return 1

fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------