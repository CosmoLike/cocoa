#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_CAMSPEC_CMB_DATA}" ]; then
  
  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi

  # parenthesis = run in a subshell 
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

  unset_env_vars () {
    unset -v EDATAF FOLDER PACKDIR FILE URL_BASE URL PRINTNAME
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

  unset_env_vars || return 1

  # E = EXTERNAL, DATA, F=FODLER
  EDATAF="${ROOTDIR:?}/external_modules/data/planck/CamSpec"
  
  FOLDER="CamSpec2021"

  # PACK = PACKAGE, DIR = DIRECTORY
  PACKDIR="${EDATAF:?}/${FOLDER:?}"

  FILE="CamSpec2021.zip"

  URL_BASE="https://github.com/CobayaSampler/planck_native_data"

  URL="${URL_BASE:?}/releases/download/v1/${FILE:?}"
  
  ptop "SETUP/UNXV CAMSPEC-2021 DATA"  || return 1

  if [ -n "${OVERWRITE_EXISTING_CAMPSPEC_CMB_DATA}" ]; then  
    rm -rf "${PACKDIR:?}"
    if [ -n "${REDOWNLOAD_EXISTING_CAMPSPEC_CMB_DATA}" ]; then  
      rm -rf "${EDATAF:?}/${FILE:?}"
    fi
  fi
  
  if [ ! -d "${PACKDIR:?}" ]; then
    cdfolder "${EDATAF:?}" || return 1

    if [ ! -e "${FILE:?}" ]; then
      "${WGET:?}" "${URL:?}" -q --show-progress \
        --progress=bar:force:noscroll || { error "${EC24:?}"; return 1; }
    fi

    unzip "${FILE:?}" >${OUT1:?} 2>${OUT2:?} || { error "${EC26:?}"; return 1; }
  fi
  
  pbottom "SETUP/UNXV CAMSPEC-2021 DATA" || return 1

  cdfolder "${ROOTDIR}" || return 1;

  unset_all || return 1

fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------