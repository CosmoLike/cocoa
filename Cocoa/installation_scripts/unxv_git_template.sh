#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_SETUP_XXX_DATA}" ]; then

  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi

  # parenthesis = run in a subshell 
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;
  
  unset_env_vars () {
    unset -v EDATAF FOLDER URL PACKDIR PRINTNAME
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
  
  URL="${XXX_DATA_URL:-"https://github.com/XXX"}"

  # FOLDER = the directory name of the dataset
  FOLDER="XXX"
  
  # Name to be printed on this shell script messages
  PRINTNAME=XXX

  # PACK = PACKAGE, DIR = DIRECTORY
  PACKDIR="${EDATAF:?}/${FOLDER:?}"

  # ---------------------------------------------------------------------------

  ptop "GETTING AND DECOMPRESSING ${PRINTNAME:?} DATA" || return 1

  # ---------------------------------------------------------------------------
  # in case this script is called twice
  
  rm -rf "${EDATAF:?}/${FOLDER:?}"
  
  # ---------------------------------------------------------------------------

  cdfolder "${EDATAF:?}" || return 1

  ${GIT:?} clone "${URL:?}" "${FOLDER:?}" \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC15:?}"; return 1; }

  # XXX_DATA_GIT_COMMIT = commit hash
  if [ -n "${XXX_DATA_GIT_COMMIT}" ]; then
    
    cdfolder "${PACKDIR:?}" || return 1

    ${GIT:?} checkout "${XXX_DATA_GIT_COMMIT:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC16:?}"; return 1; }
  
  fi

  # ---------------------------------------------------------------------------
  
  pbottom "GETTING AND DECOMPRESSING ${PRINTNAME:?} DATA" || return 1

  unset_all || return 1

fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------