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
  unset -v PRINTNAME FOLDER URL PACKAGE_VERSION PACKAGE_BRANCH
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

gitact0() {  
  local PROJECT="${ROOTDIR:?}/projects"
  local PACKDIR="${PROJECT:?}/${1:?}"
  cdfolder "${PROJECT:?}" || return 1;

  # ---------------------------------------------------------------------------
  # In case this script runs twice --------------------------------------------
  # ---------------------------------------------------------------------------
  if [ -n "${OVERWRITE_EXISTING_COSMOLIKE_CODE}" ]; then
    rm -rf "${PACKDIR:?}"
  fi

  if [ ! -d "${PACKDIR:?}" ]; then
    "${GIT:?}" clone "${2:?}" "${1:?}" >${OUT1:?} 2>${OUT2:?} || { error "${EC15:?}"; return 1; }
  fi
    
  cdfolder "${ROOTDIR:?}" || return 1;
}

gitact1() {  
  local PACKDIR="${ROOTDIR:?}/projects/${1:?}"
  local TAG="${2:?}"

  if [ -d "${PACKDIR:?}" ]; then
    cdfolder "${PACKDIR:?}" || return 1;
    if "${GIT:?}" show-ref --quiet refs/heads/${TAG}; then
      echo "git branch exists" >${OUT1:?} 2>${OUT2:?} || { return 1; }
    else
      "${GIT:?}" checkout -b ${TAG} origin/${TAG} \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi
  fi
    
  cdfolder "${ROOTDIR:?}" || return 1;
}

gitact2() {  
  local PACKDIR="${ROOTDIR:?}/projects/${1:?}"
  local TAG="${2}"

  if [ -d "${PACKDIR:?}" ]; then
    cdfolder "${PACKDIR:?}" || return 1;
    "${GIT:?}" checkout ${TAG} >${OUT1:?} 2>${OUT2:?} || { error "${EC16:?}"; return 1; }
  fi
    
  cdfolder "${ROOTDIR}" || return 1;
}

gitact3() {  
  local PACKDIR="${ROOTDIR:?}/projects/${1:?}"
  local TAG="${2}"

  if [ -d "${PACKDIR:?}" ]; then
    cdfolder "${PACKDIR:?}" || return 1;
    "${GIT:?}" fetch --all --tags --prune
    "${GIT:?}" checkout tags/${TAG} -b ${TAG} \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC16:?}"; return 1; }
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

if [ -z "${IGNORE_COSMOLIKE_LSST_Y1_CODE}" ]; then 
  
  # Name to be printed on this shell script messages
  PRINTNAME="LSST_Y1"

  ptop "GETTING ${PRINTNAME:?}" || return 1

  FOLDER="${LSST_Y1_NAME:-"lsst_y1"}"

  URL="${LSST_Y1_URL:-"https://github.com/CosmoLike/cocoa_lsst_y1.git"}"

  gitact0 "${FOLDER:?}" "${URL:?}"

  if [ -n "${LSST_Y1_COMMIT}" ]; then
    gitact2 "${FOLDER:?}" "${LSST_Y1_COMMIT:?}"  || return 1
  elif [ -n "${LSST_Y1_BRANCH}" ]; then 
    gitact1 "${FOLDER:?}" "${LSST_Y1_BRANCH:?}" || return 1
  elif [ -n "${LSST_Y1_TAG}" ]; then 
    gitact3 "${FOLDER:?}" "${LSST_Y1_TAG:?}" || return 1
  fi

  pbottom "GETTING ${PRINTNAME:?}" || return 1

fi

# ----------------------------------------------------------------------------
# -------------------------------- DES_Y3 ------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_COSMOLIKE_DES_Y3_CODE}" ]; then 
  
  # Name to be printed on this shell script messages
  PRINTNAME="DES_Y3"

  ptop "GETTING ${PRINTNAME:?}" || return 1

  FOLDER="${DES_Y3_NAME:-"des_y3"}"

  URL="${DES_Y3_URL:-"https://github.com/CosmoLike/cocoa_des_y3.git"}"

  gitact0 "${FOLDER:?}" "${URL:?}"

  if [ -n "${DES_Y3_COMMIT}" ]; then
    gitact2 "${FOLDER:?}" "${DES_Y3_COMMIT:?}"  || return 1
  elif [ -n "${DES_Y3_BRANCH}" ]; then 
    gitact1 "${FOLDER:?}" "${DES_Y3_BRANCH:?}" || return 1
  elif [ -n "${DES_Y3_TAG}" ]; then 
    gitact3 "${FOLDER:?}" "${DES_Y3_TAG:?}" || return 1
  fi

  pbottom "GETTING ${PRINTNAME:?}" || return 1

fi

# ----------------------------------------------------------------------------
# -------------------------- ROMAN_FOURIER ------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_COSMOLIKE_ROMAN_FOURIER_CODE}" ]; then 
  
  # Name to be printed on this shell script messages
  PRINTNAME="ROMAN_FOURIER"

  ptop "GETTING ${PRINTNAME:?}" || return 1

  FOLDER="${ROMAN_FOURIER_NAME:-"roman_fourier"}"

  URL="${ROMAN_FOURIER_URL:-"git@github.com:CosmoLike/cocoa_roman_fourier.git"}"

  gitact0 "${FOLDER:?}" "${URL:?}"

  if [ -n "${ROMAN_FOURIER_COMMIT}" ]; then
    gitact2 "${FOLDER:?}" "${ROMAN_FOURIER_COMMIT:?}"  || return 1
  elif [ -n "${ROMAN_FOURIER_BRANCH}" ]; then 
    gitact1 "${FOLDER:?}" "${ROMAN_FOURIER_BRANCH:?}" || return 1
  elif [ -n "${ROMAN_FOURIER_TAG}" ]; then 
    gitact3 "${FOLDER:?}" "${ROMAN_FOURIER_TAG:?}" || return 1
  fi

  pbottom "GETTING ${PRINTNAME:?}" || return 1

fi

# ----------------------------------------------------------------------------
# ---------------------------- ROMAN_REAL ------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_COSMOLIKE_ROMAN_REAL_CODE}" ]; then 
  
  # Name to be printed on this shell script messages
  PRINTNAME="ROMAN_REAL"

  ptop "GETTING ${PRINTNAME:?}" || return 1

  FOLDER="${ROMAN_REAL_NAME:-"roman_fourier"}"

  URL="${ROMAN_REAL_URL:-"git@github.com:CosmoLike/cocoa_roman_real.git"}"

  gitact0 "${FOLDER:?}" "${URL:?}"

  if [ -n "${ROMAN_REAL_COMMIT}" ]; then
    gitact2 "${FOLDER:?}" "${ROMAN_REAL_COMMIT:?}"  || return 1
  elif [ -n "${ROMAN_REAL_BRANCH}" ]; then 
    gitact1 "${FOLDER:?}" "${ROMAN_REAL_BRANCH:?}" || return 1
  elif [ -n "${ROMAN_REAL_TAG}" ]; then 
    gitact3 "${FOLDER:?}" "${ROMAN_REAL_TAG:?}" || return 1
  fi

  pbottom "GETTING ${PRINTNAME:?}" || return 1

fi

# ----------------------------------------------------------------------------
# ----------------------------- DES_Y1 x Planck ------------------------------
# ----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------