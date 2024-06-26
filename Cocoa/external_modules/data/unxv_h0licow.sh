#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${SKIP_DECOMM_STRONG_LENSING}" ]; then

  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1;
  fi

  unset_env_vars () {
    unset -v EDATAF DATAF PACKDIR URL TDATA
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
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER ${1}"; return 1; }
  }

  # --------------------------------------------------------------------------- 
  # --------------------------------------------------------------------------- 
  # ---------------------------------------------------------------------------

  unset_env_vars || return 1

  # E = EXTERNAL, DATA, F=FODLER
  EDATAF="${ROOTDIR:?}/external_modules/data/"
  
  TDATA="holicow_tmp" # = TMP

  DATAF="h0licow_distance_chains"

  # PACK = PACKAGE, DIR = DIRECTORY
  PACKDIR="${EDATAF:?}/${TDATA:?}"

  URL="${HOLICOW_DATA_URL:-"https://github.com/shsuyu/H0LiCOW-public.git"}"

  # ---------------------------------------------------------------------------

  ptop 'DECOMPRESSING H0LICOW DATA' || return 1

  # ---------------------------------------------------------------------------
  # in case this script is called twice
  
  rm -rf "${PACKDIR:?}"
  
  # ---------------------------------------------------------------------------

  cdfolder "${EDATAF:?}" || return 1

  ${GIT:?} clone "${URL:?}" "${TDATA:?}" \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC15:?}"; return 1; }

  cdfolder "${PACKDIR:?}" || return 1

  if [ -n "${HOLICOW_DATA_GIT_COMMIT}" ]; then
    ${GIT:?} checkout "${HOLICOW_DATA_GIT_COMMIT:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC16:?}"; return 1; }
  fi

  mv "${DATAF:?}" "${EDATAF:?}"

  rm -rf "${PACKDIR:?}"

  # ---------------------------------------------------------------------------

  cdfolder "${ROOTDIR}" || return 1

  unset_all || return 1
  
  pbottom 'DECOMPRESSING H0LICOW DATA' || return 1

fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------