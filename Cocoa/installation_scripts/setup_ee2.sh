#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_EUCLID_EMULATOR_V2_CODE}" ]; then

  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi

  # parenthesis = run in a subshell
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

  unset_env_vars () {
    unset -v URL CCIL ECODEF FOLDER PACKDIR CHANGES TFOLDER TFILE TFILEP AL
    cdroot || return 1;
  }

  unset_env_funcs () {
    unset -f cdfolder cpfolder error cpfile
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

  cpfolder() {
    cp -r "${1:?}" "${2:?}"  \
      2>"/dev/null" || { error "CP FOLDER ${1} on ${2}"; return 1; }
  }

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------

  unset_env_vars || return 1

  CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

  # ---------------------------------------------------------------------------
  
  URL="${EE2_URL:-"https://github.com/miknab/EuclidEmulator2.git"}"

  # E = EXTERNAL, CODE, F=FODLER
  ECODEF="${ROOTDIR:?}/external_modules/code"

  FOLDER="${EE2_NAME:-"euclidemu2"}"

  PACKDIR="${ECODEF:?}/${FOLDER:?}"

  # Name to be printed on this shell script messages
  PRINTNAME="EUCLID EMULATOR V2"

  ptop "INSTALLING ${PRINTNAME:?}" || return 1;

  # ---------------------------------------------------------------------------
  # In case this script is called twice ---------------------------------------
  # ---------------------------------------------------------------------------
  rm -rf "${PACKDIR:?}"

  # ---------------------------------------------------------------------------
  # Clone from original repo --------------------------------------------------
  # ---------------------------------------------------------------------------
  cdfolder "${ECODEF:?}" || { cdroot; return 1; }

  "${CURL:?}" -fsS "${URL:?}" \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC27:?} (URL=${URL:?})"; return 1; }

  "${GIT:?}" clone "${URL:?}" "${FOLDER:?}" \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC15:?}"; return 1; }
  
  cdfolder "${PACKDIR}" || { cdroot; return 1; }

  if [ -n "${EE2_GIT_COMMIT}" ]; then
    "${GIT:?}" checkout "${EE2_GIT_COMMIT:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC16:?}"; return 1; }
  fi
  
  cdfolder "${ROOTDIR}" || return 1
  
  pbottom "INSTALLING ${PRINTNAME:?}" || return 1
  
  # ---------------------------------------------------------------------------

  unset_all || return 1
  
fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------