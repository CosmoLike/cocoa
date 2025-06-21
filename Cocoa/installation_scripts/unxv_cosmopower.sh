#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_COSMOPOWER_DATA}" ]; then

  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi

  # parenthesis = run in a subshell 
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;
  
  unset_env_vars () {
    unset -v EDATAF FOLDER URL PACKDIR
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
  EDATAF="${ROOTDIR:?}/external_modules/data/"

  URL="${COSMOPOWER_URL_DATA:-"https://github.com/cosmopower-organization/jense_2024_emulators.git"}"

  FOLDER="cosmopower"
  
  # PACK = PACKAGE, DIR = DIRECTORY
  PACKDIR="${EDATAF:?}/${FOLDER:?}"

  # ---------------------------------------------------------------------------
  ptop "SETUP/UNXV EMULATOR CMB COSMOPOWER DATA" || return 1

  # ---------------------------------------------------------------------------
  # in case this script is called twice
  if [ -n "${OVERWRITE_EXISTING_COSMOPOWER_DATA}" ]; then
    rm -rf "${PACKDIR:?}"
  fi
  # ---------------------------------------------------------------------------
  if [[ ! -d "${PACKDIR:?}" ]]; then
    cdfolder "${EDATAF:?}" || return 1

    ${GIT:?} clone "${URL:?}" "${FOLDER:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC15:?}"; return 1; }

    if [ -n "${COSMOPOWER_URL_DATA_COMMIT}" ]; then
      cdfolder "${PACKDIR:?}" || return 1

      ${GIT:?} checkout "${COSMOPOWER_URL_DATA_COMMIT:?}" \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi
  fi
  # ---------------------------------------------------------------------------
  
  pbottom "SETUP/UNXV EMULATOR CMB COSMOPOWER DATA" || return 1

  unset_all || return 1

fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------