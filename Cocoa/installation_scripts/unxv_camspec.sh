#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${SKIP_DECOMM_CAMSPEC}" ]; then
  
  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi

  # parenthesis = run in a subshell 
  ( source "${ROOTDIR:?}/installation_scripts/.check_flags.sh" ) || return 1;

  unset_env_vars () {
    unset -v EDATAF FOLDER PACKDIR FILE URL_BASE URL
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
  EDATAF="${ROOTDIR:?}/external_modules/data/planck/CamSpec"
  
  FOLDER="CamSpec2021"

  # PACK = PACKAGE, DIR = DIRECTORY
  PACKDIR="${EDATAF:?}/${FOLDER:?}"

  FILE="CamSpec2021.zip"

  URL_BASE="https://github.com/CobayaSampler/planck_native_data/"

  URL="${URL_BASE:?}/releases/download/v1/${FILE:?}"
  
  # ---------------------------------------------------------------------------

  ptop "GETTING AND DECOMPRESSING CAMSPEC DATA (MAY TAKE A LONG TIME)" || return 1

  # ---------------------------------------------------------------------------
  # note: in case script run >1x w/ previous run stoped prematurely b/c error
  
  rm -rf "${PACKDIR:?}"

  rm -rf "${EDATAF:?}/${FILE:?}"
  
  # ---------------------------------------------------------------------------

  cdfolder "${EDATAF:?}" || return 1

  "${WGET:?}" "${URL:?}" -q --show-progress \
    --progress=bar:force:noscroll || { error "${EC24:?}"; return 1; }
  
  unzip "${FILE:?}" >${OUT1:?} 2>${OUT2:?} || { error "${EC26:?}"; return 1; }

  # ---------------------------------------------------------------------------

  unset_all || return 1
  
  pbottom 'GETTING AND DECOMPRESSING CAMSPEC DATA' || return 1

fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------