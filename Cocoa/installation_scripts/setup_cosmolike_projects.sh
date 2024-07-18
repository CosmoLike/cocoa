#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${ROOTDIR}" ]; then
  pfail 'ROOTDIR'; return 1
fi

# parenthesis = run in a subshell
( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" )  || return 1;

unset_env_vars () {
  unset -v PRINTNAME FOLDER URL PACKAGE_VERSION
  cdroot || return 1;
}

unset_env_funcs () {
  unset -f cdfolder cpfolder cpfile error gitact gitact1 gitact2
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

cpfile() {
  cp "${1:?}" "${2:?}" \
    2>"/dev/null" || { error "CP FILE ${1} on ${2}"; return 1; }
}

gitact() {
  #ARGUMENTS: FOLDER, VERSION, URL
  
  local PROJECT="${ROOTDIR}/projects"

  local PACKDIR="${PROJECT:?}/${1:?}"

  cdfolder "${PROJECT:?}" || return 1;

  # ---------------------------------------------------------------------------
  # In case this script runs twice --------------------------------------------
  # ---------------------------------------------------------------------------
  if [ -n "${OVERWRITE_EXISTING_COSMOLIKE_CODE}" ]; then
    
    rm -rf "${PACKDIR:?}"
  
  fi

  # ---------------------------------------------------------------------------
  # clone from original repo --------------------------------------------------
  # ---------------------------------------------------------------------------
  if [ ! -d "${PACKDIR:?}" ]; then

    "${CURL:?}" -fsS "${3:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC27:?} (URL=${3:?})"; return 1; }

    "${GIT:?}" clone "${3:?}" --branch "${2:?}" "${1:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error "GIT CLONE"; return 1; }

  fi
    
  cdfolder "${ROOTDIR}" || return 1;
}

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

unset_env_vars || return 1

# ----------------------------------------------------------------------------
# -------------------------------- LSST-Y1 -----------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_COSMOLIKE_LSSTY1_CODE}" ]; then 
  
  # Name to be printed on this shell script messages
  PRINTNAME="LSST_Y1"

  ptop "GETTING ${PRINTNAME:?}" || return 1

  FOLDER="${LSSTY1_NAME:-"lsst_y1"}"

  URL="${LSSTY1_URL:-"https://github.com/CosmoLike/cocoa_lsst_y1.git"}"

  PACKAGE_VERSION=${LSSTY1_COMMIT:?}

  gitact "${FOLDER:?}" "${URL:?}" "${PACKAGE_VERSION:?}"

  pbottom "GETTING ${PRINTNAME:?}" || return 1

fi

# TODO:

# ----------------------------------------------------------------------------
# -------------------------------- DES-Y3 ------------------------------------
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# ----------------------------- DES-Y1 x Planck ------------------------------
# ----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------