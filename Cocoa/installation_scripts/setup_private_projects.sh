#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${ROOTDIR:-}" ]; then
  pfail 'ROOTDIR'; return 1
fi

# parenthesis = run in a subshell
( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" )  || return 1;

unset_env_vars () {
  unset -v PRINTNAME FOLDER URL PACKAGE_VERSION PACKAGE_BRANCH
  cdroot || return 1;
}

unset_env_funcs () {
  unset -f cdfolder cpfolder cpfile error gitact0 gitact1 gitact2 gitact3 devurl
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
  if [ -n "${OVERWRITE_EXISTING_PRIVATE_PROJECTS_CODE}" ]; then
    rm -rf "${PACKDIR:?}"
  fi

  if [ ! -d "${PACKDIR:?}" ]; then
  
    "${GIT:?}" clone "${URL}" --depth ${GIT_CLONE_MAXIMUM_DEPTH:-1000} \
      --no-single-branch --recursive "${NAME}" \
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
  if [ -n "${OVERWRITE_EXISTING_PRIVATE_PROJECTS_CODE:-}" ]; then
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
  local COMMIT="${2}"
  local ARGS="--all --tags --prune"
  
  if [ -d "${PACKDIR:?}" ]; then
    cdfolder "${PACKDIR:?}" || { unset_all; return 1; }
    
    # unshallow the repo if necessary before check out a specific commit/tag
    if [ "$("${GIT:?}" rev-parse --is-shallow-repository)" = "true" ]; then
      ARGS="--unshallow ${ARGS}"
    fi
    "${GIT:?}" fetch ${ARGS:?} \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }

    # in case it this script runs twice:
    # check if the current commit = wanted commit
    local CURRENT=$("${GIT:?}" rev-parse HEAD 2>/dev/null)
    local TARGET=$("${GIT:?}" rev-parse "${COMMIT}" 2>/dev/null)
    if [ "${CURRENT}" != "${TARGET}" ]; then
      # first check unstagged, second check staged but uncommitted work
      # only checkout on a second run if there isn't uncommitted work
      if "${GIT:?}" diff --quiet HEAD 2>/dev/null && \
         "${GIT:?}" diff --cached --quiet HEAD 2>/dev/null; then
        "${GIT:?}" checkout "${COMMIT}" \
          >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
      else
        local tmp="skipping checkout: uncommitted changes in ${PACKDIR}"
        warning_script_msg "$(basename "${BASH_SOURCE[0]}")" "${tmp:?}"
      fi
    fi
  
  fi
    
  cdfolder "${ROOTDIR}" || { unset_all; return 1; }
}

gitact3() {  
  local PACKDIR="${ROOTDIR:?}/projects/${1:?}"
  local TAG="${2}"
  local ARGS="--all --tags --prune"

  if [ -d "${PACKDIR:?}" ]; then
    cdfolder "${PACKDIR:?}" || { unset_all; return 1; }

    # unshallow the repo if necessary before check out a specific commit/tag
    if [ "$("${GIT:?}" rev-parse --is-shallow-repository)" = "true" ]; then
      ARGS="--unshallow ${ARGS}"
    fi
    "${GIT:?}" fetch ${ARGS:?} \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    
    if "${GIT:?}" show-ref --verify --quiet "refs/heads/${TAG:?}TMP"; then
      # in case it this script runs twice (CB = CURRENT BRANCH)
      # check if the branch = tag name. 
      # why this works? --branch creates a local branch with the TAG name
      local CB=$("${GIT:?}" rev-parse --abbrev-ref HEAD 2>/dev/null)
      CB="${CB##*/}"
      if [ "${CB}" != "${TAG}TMP" ]; then
        # first check unstagged, second check staged but uncommitted work
        # only checkout on a second run if there isn't uncommitted work
        if "${GIT:?}" diff --quiet HEAD 2>/dev/null && \
           "${GIT:?}" diff --cached --quiet HEAD 2>/dev/null; then
          "${GIT:?}" checkout "${TAG:?}TMP" \
              >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
        else
          local tmp="skipping checkout: uncommitted changes in ${PACKDIR}"
          warning_script_msg "$(basename "${BASH_SOURCE[0]}")" "${tmp:?}"
        fi
      fi
    else
      "${GIT:?}" checkout "tags/${TAG:?}" -b "${TAG:?}TMP" \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi
  fi
    
  cdfolder "${ROOTDIR}" || { unset_all; return 1; }
}

devurl() {
  # SWITCH_TO_DEV_MODE=1: rewrite GitHub https URLs to their ssh form
  local U="${1:?}"
  if [ -n "${SWITCH_TO_DEV_MODE:-}" ]; then
    case "${U}" in
      https://*github.com/*) U="git@github.com:${U#*github.com/}" ;;
    esac
  fi
  echo "${U}"
}

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

unset_env_vars || return 1

# ----------------------------------------------------------------------------
# -------------------------------- AXIONS -----------------------------------
# ----------------------------------------------------------------------------
if [ -n "${INSTALL_AXIONS_2025_PROJECT:-}" ]; then 

  # Name to be printed on this shell script messages
  PRINTNAME="AXIONS PROJECT"

  ptop "GETTING ${PRINTNAME:?}" || { unset_all; return 1; }

  FOLDER="${AXIONS_PROJECT_NAME:-"axions"}"

  URL="${AXIONS_2025_PROJECT_URL:?}"
  URL=$(devurl "${URL:?}")

  if [ -n "${AXIONS_2025_PROJECT_GIT_COMMIT:-}" ]; then

    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }
    
    gitact2 "${FOLDER:?}" "${AXIONS_2025_PROJECT_GIT_COMMIT:?}" || { unset_all; return 1; }
  
  elif [ -n "${AXIONS_2025_PROJECT_GIT_BRANCH:-}" ]; then 
  
    gitact1 "${FOLDER:?}" "${URL:?}" "${AXIONS_2025_PROJECT_GIT_BRANCH:?}" || { unset_all; return 1; }
  
  elif [ -n "${AXIONS_2025_PROJECT_GIT_TAG:-}" ]; then 
  
    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }
    
    gitact3 "${FOLDER:?}" "${AXIONS_2025_PROJECT_GIT_TAG:?}" || { unset_all; return 1; }
  
  else
  
    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }
  fi

  pbottom "GETTING ${PRINTNAME:?}" || { unset_all; return 1; }

fi

# ----------------------------------------------------------------------------
# -------------------------------- INFPC -------------------------------------
# ----------------------------------------------------------------------------

if [ -n "${INSTALL_PRIVATE_INFPC_PROJECT:-}" ]; then 
 
  # Name to be printed on this shell script messages
  PRINTNAME="INFPC PROJECT"

  ptop "GETTING ${PRINTNAME:?}" || { unset_all; return 1; }

  FOLDER="${INFPC_PROJECT_NAME:-"InPCReiPC"}"

  URL="${INFPC_PROJECT_URL:-"git@github.com:SBU-COSMOLIKE/InPCReiPCNewCAMB.git"}"
  URL=$(devurl "${URL:?}")

  if [ -n "${INFPC_PROJECT_GIT_COMMIT:-}" ]; then

    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }

    gitact2 "${FOLDER:?}" "${INFPC_PROJECT_GIT_COMMIT:?}"  || { unset_all; return 1; }

  elif [ -n "${INFPC_PROJECT_GIT_BRANCH:-}" ]; then 

    gitact1 "${FOLDER:?}" "${URL:?}" "${INFPC_PROJECT_GIT_BRANCH:?}" || { unset_all; return 1; }

  elif [ -n "${INFPC_PROJECT_GIT_TAG:-}" ]; then 

    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }

    gitact3 "${FOLDER:?}" "${INFPC_PROJECT_GIT_TAG:?}" || { unset_all; return 1; }

  else

    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }

  fi

  pbottom "GETTING ${PRINTNAME:?}" || { unset_all; return 1; }

fi
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

if [ -n "${INSTALL_AXIE_CAMB_2026_PROJECT:-}" ]; then 
  # Name to be printed on this shell script messages
  PRINTNAME="AXIE CAMB 2026 PROJECT"

  ptop "GETTING ${PRINTNAME:?}" || { unset_all; return 1; }

  FOLDER="${AXIE_CAMB_2026_PROJECT_NAME:-"axicambv2"}"

  URL="${AXIE_CAMB_2026_PROJECT_URL:?}"
  URL=$(devurl "${URL:?}")

  if [ -n "${AXIE_CAMB_2026_PROJECT_GIT_COMMIT:-}" ]; then

    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }

    gitact2 "${FOLDER:?}" "${AXIE_CAMB_2026_PROJECT_GIT_COMMIT:?}"  || { unset_all; return 1; }
  
  elif [ -n "${AXIE_CAMB_2026_PROJECT_GIT_BRANCH:-}" ]; then 
  
    gitact1 "${FOLDER:?}" "${URL:?}" "${AXIE_CAMB_2026_PROJECT_GIT_BRANCH:?}" || { unset_all; return 1; }
  
  elif [ -n "${AXIE_CAMB_2026_PROJECT_GIT_TAG:-}" ]; then 
  
    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }
  
    gitact3 "${FOLDER:?}" "${AXIE_CAMB_2026_PROJECT_GIT_TAG:?}" || { unset_all; return 1; }
  
  else
  
    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }
  
  fi

  pbottom "GETTING ${PRINTNAME:?}" || { unset_all; return 1; }
fi

unset_all || return 1;

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

return 55; # why this odd number? Setup_cocoa will cache this installation only
           #   if this script runs entirely. What if the user close the terminal 
           #   or the system shuts down in the middle of a git clone?  
           #   In this case, PACKDIR would exists, but it is corrupted

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------