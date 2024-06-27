#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${SKIP_DECOMM_ACT}" ]; then
  
  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi

  # parenthesis = run in a subshell 
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

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
  EDATAF="${ROOTDIR:?}/external_modules/data"
  
  FOLDER="act"

  # PACK = PACKAGE, DIR = DIRECTORY
  PACKDIR="${EDATAF:?}/${FOLDER:?}"

  FILE="${ACT_DR6_DATA_FILE:-"ACT_dr6_likelihood_v1.2.tgz"}"

  URL_BASE="https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr6/likelihood/data"

  URL="${ACT_DR6_DATA_URL:-"${URL_BASE:?}"}/${FILE:?}"

  # ---------------------------------------------------------------------------

  ptop 'DECOMPRESSING ACT-DR6 DATA (MAY TAKE A LONG TIME)' || return 1

  # ---------------------------------------------------------------------------
  # note: in case script run >1x w/ previous run stoped prematurely b/c error
  
  rm -rf "${PACKDIR:?}"

  # ---------------------------------------------------------------------------

  mkdir -p "${PACKDIR:?}" || { error "${EC20:?}"; return 1; }
  
  mkdir -p "${PACKDIR:?}/lensing" || { error "${EC20:?}"; return 1; }
    
  cdfolder "${PACKDIR:?}/lensing" || return 1

  "${WGET:?}" "${URL:?}" -q --show-progress --progress=bar:force || 
    { error "${EC24:?}"; return 1; }

  tar -zxvf "${FILE:?}" \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC25:?}"; return 1; }

  # ---------------------------------------------------------------------------

  unset_all || return 1
  
  pbottom 'DECOMPRESSING ACT-DR6 DATA' || return 1

fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------