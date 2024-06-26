#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${SKIP_DECOMM_SIMONS_OBSERVATORY}" ]; then

  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1;
  fi

  # parenthesis = run in a subshell 
  ( source "${ROOTDIR:?}/installation_scripts/.check_flags.sh" ) || return 1;
    
  unset_env_vars () {
    unset -v EDATAF FOLDER URL VER PACKDIR FILE
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

  FOLDER="simons_observatory"

  URL=${SO_DATA_URL:-"https://portal.nersc.gov/cfs/sobs/users/MFLike_data"}

  VER="${SO_DATA_VERSION:-"v0.8"}"

  # PACK = PACKAGE, DIR = DIRECTORY
  PACKDIR="${EDATAF:?}/${FOLDER:?}"

  # --------------------------------------------------------------------------

  ptop "GETTING AND DECOMPRESSING S.O. DATA (MAY TAKE A LONG TIME)" || return 1

  # --------------------------------------------------------------------------
  # note: in case script run >1x w/ previous run stoped prematurely b/c error
  
  rm -rf "${PACKDIR:?}"
  
  # ---------------------------------------------------------------------------

  mkdir -p "${PACKDIR:?}"
  
  cdfolder "${PACKDIR:?}" || return 1

  # note: users can download multiple versions (reproduce existing work)
  # note: For example, SO_DATA_VERSION="v0.7.1 v0.8"
  # note: This is only possible because each ver is saved on a separated folder
  for x in $(echo "${VER:?}")
  do
    FILE="${x}.tar.gz"

    "${WGET:?}" "${URL}/${FILE:?}" -q --show-progress --progress=bar:force \
      || { error "${EC24:?}"; return 1; }

    tar -zxvf "${FILE:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC25:?}"; return 1; }
  
  done

  # ---------------------------------------------------------------------------

  unset_all || return 1

  pbottom "GETTING AND DECOMPRESSING S.O. DATA" || return 1

fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------