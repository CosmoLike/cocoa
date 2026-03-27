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
  
  cdfolder "${PROJECT:?}" || { unset_all; return 1; }

  # ---------------------------------------------------------------------------
  # In case this script runs twice --------------------------------------------
  # ---------------------------------------------------------------------------
  if [ -n "${OVERWRITE_EXISTING_COSMOLIKE_CODE}" ]; then
    rm -rf "${PACKDIR:?}"
  fi

  if [ ! -d "${PACKDIR:?}" ]; then
  
    "${GIT:?}" clone "${URL}" --depth ${GIT_CLONE_MAXIMUM_DEPTH:-1000} \
      --recursive "${NAME}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }
  
  fi
    
  cdfolder "${ROOTDIR:?}" || { unset_all; return 1; }
}

gitact1() { 
  local PROJECT="${ROOTDIR:?}/projects" 
  local NAME="${1:?}"
  local PACKDIR="${ROOTDIR:?}/projects/${NAME}"
  local URL="${2:?}"
  local TAG="${3:?}"
  
  cdfolder "${PROJECT:?}" || { unset_all; return 1; }
  
  # ---------------------------------------------------------------------------
  # In case this script runs twice --------------------------------------------
  # ---------------------------------------------------------------------------
  if [ -n "${OVERWRITE_EXISTING_COSMOLIKE_CODE:-}" ]; then
  
    rm -rf "${PACKDIR:?}"
  
  fi

  if [ ! -d "${PACKDIR:?}" ]; then
  
    "${GIT:?}" clone "${URL}" "${NAME}" --branch "${TAG}" --single-branch \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }
  
  fi
    
  cdfolder "${ROOTDIR:?}" || { unset_all; return 1; }
}

gitact2() {  
  local PACKDIR="${ROOTDIR:?}/projects/${1:?}"
  local TAG="${2}"

  if [ -d "${PACKDIR:?}" ]; then
    cdfolder "${PACKDIR:?}" || { unset_all; return 1; }

    if [ "$("${GIT:?}" rev-parse --is-shallow-repository)" = "true" ]; then
    
      "${GIT:?}" fetch --unshallow --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    
    else
    
      "${GIT:?}" fetch --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    
    fi

    "${GIT:?}" checkout ${TAG} \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  fi
    
  cdfolder "${ROOTDIR}" || { unset_all; return 1; }
}

gitact3() {  
  local PACKDIR="${ROOTDIR:?}/projects/${1:?}"
  local TAG="${2}"

  if [ -d "${PACKDIR:?}" ]; then
    cdfolder "${PACKDIR:?}" || { unset_all; return 1; }
    
    if [ "$("${GIT:?}" rev-parse --is-shallow-repository)" = "true" ]; then
   
      "${GIT:?}" fetch --unshallow --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
   
    else
   
      "${GIT:?}" fetch --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
   
    fi
    
    if "${GIT:?}" show-ref --verify --quiet "refs/heads/${TAG:?}"; then
   
      "${GIT:?}" checkout "${TAG:?}" \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
   
    else
   
      "${GIT:?}" checkout "tags/${TAG:?}" -b "${TAG:?}" \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
   
    fi

  fi
    
  cdfolder "${ROOTDIR}" || { unset_all; return 1; }
}


# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

unset_env_vars || return 1

# ----------------------------------------------------------------------------
# -------------------------------- AXIONS -----------------------------------
# ----------------------------------------------------------------------------

if [ -n "${INSTALL_PRIVATE_AXIONS_PROJECT:-}" ]; then 
  
  # Name to be printed on this shell script messages
  PRINTNAME="AXIONS PROJECT"

  ptop "GETTING ${PRINTNAME:?}" || { unset_all; return 1; }

  FOLDER="${AXIONS_PROJECT_NAME:-"axions"}"

  URL="${AXIONS_PROJECT_URL:?}"

  if [ -n "${AXIONS_PROJECT_GIT_COMMIT:-}" ]; then

    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }
  
    gitact2 "${FOLDER:?}" "${AXIONS_PROJECT_GIT_COMMIT:?}"  || { unset_all; return 1; }
  
  elif [ -n "${AXIONS_PROJECT_GIT_BRANCH:-}" ]; then 
  
    gitact1 "${FOLDER:?}" "${URL:?}" "${AXIONS_PROJECT_GIT_BRANCH:?}" || { unset_all; return 1; }
  
  elif [ -n "${AXIONS_PROJECT_GIT_TAG:-}" ]; then 
  
    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }
  
    gitact3 "${FOLDER:?}" "${AXIONS_PROJECT_GIT_TAG:?}" || { unset_all; return 1; }
  
  else
  
    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }
  
  fi

  pbottom "GETTING ${PRINTNAME:?}" || { unset_all; return 1; }

fi

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------