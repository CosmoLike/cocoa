#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_SETUP_SPT_CMB_DATA}" ]; then

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
  EDATAF="${ROOTDIR:?}/external_modules/data"

  FOLDER="${SPT_3G_NAME:-"spt_3g"}"

  # PACK = PACKAGE, DIR = DIRECTORY
  PACKDIR="${EDATAF:?}/${FOLDER:?}"

  URL="${SPT3G_DATA_URL:-"https://github.com/SouthPoleTelescope/spt3g_y1_dist.git"}"

  # ---------------------------------------------------------------------------

  ptop "GETTING AND DECOMPRESSING SPT-3G Y1 DATA" || return 1
  
  # ---------------------------------------------------------------------------
  # note: in case this script is run twice
  
  rm -rf "${PACKDIR:?}"
  
  # ---------------------------------------------------------------------------

  cdfolder "${EDATAF:?}" || return 1

  ${GIT:?} clone "${URL:?}" "${FOLDER:?}" \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC15:?}"; return 1; }

  if [ -n "${SPT3G_DATA_GIT_COMMIT}" ]; then
    
    cdfolder "${PACKDIR:?}"

    ${GIT:?} checkout "${SPT3G_DATA_GIT_COMMIT:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC16:?}"; return 1; }
  
  fi

  # ---------------------------------------------------------------------------

  unset_all || return 1; 

  pbottom 'GETTING AND DECOMPRESSING SPT-3G Y1 DATA' || return 1
  
fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------