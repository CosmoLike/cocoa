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
  unset -f cdfolder cpfolder cpfile error gitact0 gitact1 gitact2 gitact3
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
  local NAME="${1:?}"
  local URL="${2:?}"
  local PACKDIR="${PROJECT:?}/${NAME}"
  cdfolder "${PROJECT:?}" || return 1;

  # ---------------------------------------------------------------------------
  # In case this script runs twice --------------------------------------------
  # ---------------------------------------------------------------------------
  if [ -n "${OVERWRITE_EXISTING_COSMOLIKE_CODE}" ]; then
    rm -rf "${PACKDIR:?}"
  fi

  if [ ! -d "${PACKDIR:?}" ]; then
    "${GIT:?}" clone "${URL}" "${NAME}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }
  fi
    
  cdfolder "${ROOTDIR:?}" || return 1;
}

gitact1() { 
  local PROJECT="${ROOTDIR:?}/projects" 
  local NAME="${1:?}"
  local PACKDIR="${ROOTDIR:?}/projects/${NAME}"
  local URL="${2:?}"
  local TAG="${3:?}"
  
  cdfolder "${PROJECT:?}" || return 1;
  
  # ---------------------------------------------------------------------------
  # In case this script runs twice --------------------------------------------
  # ---------------------------------------------------------------------------
  if [ -n "${OVERWRITE_EXISTING_COSMOLIKE_CODE}" ]; then
    rm -rf "${PACKDIR:?}"
  fi

  if [ ! -d "${PACKDIR:?}" ]; then
    "${GIT:?}" clone "${URL}" "${NAME}" --branch "${TAG}" --single-branch \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }
  fi
    
  cdfolder "${ROOTDIR:?}" || return 1;
}

gitact2() {  
  local PACKDIR="${ROOTDIR:?}/projects/${1:?}"
  local TAG="${2}"

  if [ -d "${PACKDIR:?}" ]; then
    cdfolder "${PACKDIR:?}" || return 1;
    "${GIT:?}" checkout ${TAG} \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
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
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  fi
    
  cdfolder "${ROOTDIR}" || return 1;
}


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

unset_env_vars || return 1

# ----------------------------------------------------------------------------
# -------------------------------- AXIONS -----------------------------------
# ----------------------------------------------------------------------------

if [ -n "${INSTALL_PRIVATE_AXIONS_PROJECT}" ]; then 
  
  # Name to be printed on this shell script messages
  PRINTNAME="AXIONS PROJECT"

  ptop "GETTING ${PRINTNAME:?}" || return 1

  FOLDER="${AXIONS_PROJECT_NAME:-"axions"}"

  URL="${AXIONS_PROJECT_URL:?}"

  if [ -n "${AXIONS_PROJECT_GIT_COMMIT}" ]; then
    gitact0 "${FOLDER:?}" "${URL:?}"
    gitact2 "${FOLDER:?}" "${AXIONS_PROJECT_GIT_COMMIT:?}"  || return 1
  elif [ -n "${AXIONS_PROJECT_GIT_BRANCH}" ]; then 
    gitact1 "${FOLDER:?}" "${URL:?}" "${AXIONS_PROJECT_GIT_BRANCH:?}" || return 1
  elif [ -n "${AXIONS_PROJECT_GIT_TAG}" ]; then 
    gitact0 "${FOLDER:?}" "${URL:?}"
    gitact3 "${FOLDER:?}" "${AXIONS_PROJECT_GIT_TAG:?}" || return 1
  else
    gitact0 "${FOLDER:?}" "${URL:?}"
  fi

  pbottom "GETTING ${PRINTNAME:?}" || return 1

fi

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------