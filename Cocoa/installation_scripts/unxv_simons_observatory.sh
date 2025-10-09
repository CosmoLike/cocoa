#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_SIMONS_OBSERVATORY_CMB_DATA}" ]; then

  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi

  # parenthesis = run in a subshell 
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;
    
  unset_env_vars () {
    unset -v EDATAF FOLDER URL VER PACKDIR FILE PRINTNAME
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
  EDATAF="${ROOTDIR:?}/external_modules/data"

  FOLDER="simons_observatory"

  URL=${SO_DATA_URL:-"https://portal.nersc.gov/cfs/sobs/users/MFLike_data"}

  VER="${SO_DATA_VERSION:-"v0.8"}"

  # PACK = PACKAGE, DIR = DIRECTORY
  PACKDIR="${EDATAF:?}/${FOLDER:?}"

  ptop "SETUP/UNXV SIMONS OBSERVATORY DATA" || return 1

  if [ -n "${OVERWRITE_EXISTING_SIMONS_OBSERVATORY_CMB_DATA}" ]; then
    rm -rf "${PACKDIR:?}"
    if [ -n "${REDOWNLOAD_EXISTING_SIMONS_OBSERVATORY_CMB_DATA}" ]; then
      for x in $(echo "${VER:?}")
      do
        FILE="${x}.tar.gz"
        rm -f "${EDATAF:?}/${FILE:?}"
      done
    fi
  fi
  
  if [ ! -d "${PACKDIR:?}" ]; then
    mkdir -p "${PACKDIR:?}" >${OUT1:?} 2>>${OUT2:?} || { error "${EC20:?}"; return 1; }
    
    cdfolder "${EDATAF:?}" || return 1
  
    # note: users can download multiple versions (reproduce existing work)
    # note: For example, SO_DATA_VERSION="v0.7.1 v0.8"
    for x in $(echo "${VER:?}")
    do
      FILE="${x}.tar.gz"

      if [ ! -e "${FILE:?}" ]; then
        "${WGET:?}" "${URL}/${FILE:?}" -q --show-progress --progress=bar:force || { error "${EC24:?}"; return 1; }
      fi

      TMP=$(tar -tf "${FILE:?}" | head -1 | cut -f1 -d"/")
      tar -zxvf "${FILE:?}" >${OUT1:?} 2>>${OUT2:?} || { error "${EC25:?}"; return 1; }
      mv "${TMP:?}" "${PACKDIR:?}"
    done
  fi

  pbottom "SETUP/UNXV SIMONS OBSERVATORY DATA" || return 1
  
  cdfolder "${ROOTDIR}" || return 1;

  unset_all || return 1
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------