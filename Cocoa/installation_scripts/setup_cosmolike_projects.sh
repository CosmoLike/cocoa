#!/bin/bash
# ------------------------------------------------------------------------------
# setup_cosmolike_projects.sh: download the Cosmolike project repositories
# selected in set_installation_options.sh.
#
# Sourced by setup_cocoa.sh, never executed directly. Setup scripts MAY use
# the internet (git clone, wget, pip index); the matching compile step never
# does, so anything that must be downloaded happens here.
#
# Reruns are safe: existing results are kept unless
# OVERWRITE_EXISTING_COSMOLIKE_CODE is set. On failure error() names the
# failing step and unsets everything defined here; on success the script
# returns 55, the value the runner caches so a finished step is not repeated.
# ------------------------------------------------------------------------------
if [ -z "${ROOTDIR:-}" ]; then
  pfail 'ROOTDIR'; return 1
fi

# parenthesis = run in a subshell
( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" )  || return 1;

# unset_env_vars: forget every variable the body defines. Sourced scripts
# share the caller's shell, so anything not unset here would leak into the
# user's environment. cdroot returns to the directory the user launched from.
unset_env_vars () {
  unset -v PRINTNAME FOLDER URL PACKAGE_VERSION PACKAGE_BRANCH
  cdroot || return 1;
}

# unset_env_funcs: forget every helper function defined below (functions leak
# into the user's shell exactly like variables).
unset_env_funcs () {
  unset -f cdfolder cpfolder cpfile error gitact0 gitact1 gitact2 gitact3 devurl
  unset -f unset_env_funcs
  cdroot || return 1;
}

# unset_all: the single cleanup entry point: variables, then functions, then
# itself.
unset_all () {
  unset_env_vars
  unset_env_funcs
  unset -f unset_all
  cdroot || return 1;
}

# error: print which step failed (fail_script_msg names this script and the
# message) and run the FULL cleanup. A call site whose helper already routes
# through error() therefore uses plain `|| return 1`; adding another unset_all
# there would run the cleanup twice.
error () {
  fail_script_msg "$(basename "${BASH_SOURCE[0]}")" "${1}"
  unset_all || return 1
} 

# cdfolder: cd that reports the target and cleans up on failure.
cdfolder() {
  cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
}

# cpfolder: recursive copy that reports source/target and cleans up on
# failure.
cpfolder() {
  cp -r "${1:?}" "${2:?}"  \
    2>"/dev/null" || { error "CP FOLDER ${1} on ${2}"; return 1; }
}

# cpfile: file copy that reports source/target and cleans up on failure.
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
  if [ -n "${OVERWRITE_EXISTING_COSMOLIKE_CODE:-}" ]; then
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
  if [ -n "${OVERWRITE_EXISTING_COSMOLIKE_CODE:-}" ]; then
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
  local COMMIT="${2}"
  local ARGS="--all --tags --prune"
  
  if [ -d "${PACKDIR:?}" ]; then
    cdfolder "${PACKDIR:?}" || return 1;
    
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
    
  cdfolder "${ROOTDIR}" || return 1;
}

gitact3() {  
  local PACKDIR="${ROOTDIR:?}/projects/${1:?}"
  local TAG="${2}"
  local ARGS="--all --tags --prune"

  if [ -d "${PACKDIR:?}" ]; then
    cdfolder "${PACKDIR:?}" || return 1;

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
    
  cdfolder "${ROOTDIR}" || return 1;
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
# -------------------------------- LSST-Y1 -----------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_COSMOLIKE_LSST_Y1_CODE}" ]; then 
  
  # Name to be printed on this shell script messages
  PRINTNAME="LSST_Y1"

  ptop "GETTING ${PRINTNAME:?}" || { unset_all; return 1; }

  FOLDER="${LSST_Y1_NAME:-"lsst_y1"}"

  URL="${LSST_Y1_URL:-"https://github.com/CosmoLike/cocoa_lsst_y1.git"}"
  URL=$(devurl "${URL:?}")

  if [ -n "${LSST_Y1_GIT_COMMIT:-}" ]; then
  
    gitact0 "${FOLDER:?}" "${URL:?}" || return 1;
  
    gitact2 "${FOLDER:?}" "${LSST_Y1_GIT_COMMIT:?}" || return 1;
  
  elif [ -n "${LSST_Y1_GIT_BRANCH:-}" ]; then 
  
    gitact1 "${FOLDER:?}" "${URL:?}" "${LSST_Y1_GIT_BRANCH:?}" || return 1;
  
  elif [ -n "${LSST_Y1_GIT_TAG:-}" ]; then 
  
    gitact0 "${FOLDER:?}" "${URL:?}" || return 1;
  
    gitact3 "${FOLDER:?}" "${LSST_Y1_GIT_TAG:?}" || return 1;
  
  else
  
    gitact0 "${FOLDER:?}" "${URL:?}" || return 1;
  
  fi

  pbottom "GETTING ${PRINTNAME:?}" || { unset_all; return 1; }

fi

# ----------------------------------------------------------------------------
# -------------------------------- DES_Y3 ------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_COSMOLIKE_DES_Y3_CODE}" ]; then 
  
  # Name to be printed on this shell script messages
  PRINTNAME="DES_Y3"

  ptop "GETTING ${PRINTNAME:?}" || { unset_all; return 1; }

  FOLDER="${DES_Y3_NAME:-"des_y3"}"

  URL="${DES_Y3_URL:-"https://github.com/CosmoLike/cocoa_des_y3.git"}"
  URL=$(devurl "${URL:?}")

  if [ -n "${DES_Y3_GIT_COMMIT:-}" ]; then
  
    gitact0 "${FOLDER:?}" "${URL:?}" || return 1;
  
    gitact2 "${FOLDER:?}" "${DES_Y3_GIT_COMMIT:?}" || return 1;
  
  elif [ -n "${DES_Y3_GIT_BRANCH:-}" ]; then 
  
    gitact1 "${FOLDER:?}" "${URL:?}" "${DES_Y3_GIT_BRANCH:?}" || return 1;
  
  elif [ -n "${DES_Y3_GIT_TAG:-}" ]; then 
  
    gitact0 "${FOLDER:?}" "${URL:?}" || return 1;
  
    gitact3 "${FOLDER:?}" "${DES_Y3_GIT_TAG:?}" || return 1;
  
  else
  
    gitact0 "${FOLDER:?}" "${URL:?}" || return 1;
  
  fi

  pbottom "GETTING ${PRINTNAME:?}" || { unset_all; return 1; }

fi

# ----------------------------------------------------------------------------
# -------------------------- ROMAN_FOURIER ------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_COSMOLIKE_ROMAN_FOURIER_CODE}" ]; then 
  
  # Name to be printed on this shell script messages
  PRINTNAME="ROMAN_FOURIER"

  ptop "GETTING ${PRINTNAME:?}" || { unset_all; return 1; }

  FOLDER="${ROMAN_FOURIER_NAME:-"roman_fourier"}"

  URL="${ROMAN_FOURIER_URL:-"https://github.com/CosmoLike/cocoa_roman_fourier.git"}"
  URL=$(devurl "${URL:?}")

  if [ -n "${ROMAN_FOURIER_GIT_COMMIT:-}" ]; then
  
    gitact0 "${FOLDER:?}" "${URL:?}" || return 1;
   
    gitact2 "${FOLDER:?}" "${ROMAN_FOURIER_GIT_COMMIT:?}" || return 1;
  
  elif [ -n "${ROMAN_FOURIER_GIT_BRANCH:-}" ]; then 
  
    gitact1 "${FOLDER:?}" "${URL:?}" "${ROMAN_FOURIER_GIT_BRANCH:?}" || return 1;
  
  elif [ -n "${ROMAN_FOURIER_GIT_TAG:-}" ]; then 
  
    gitact0 "${FOLDER:?}" "${URL:?}" || return 1;
   
    gitact3 "${FOLDER:?}" "${ROMAN_FOURIER_GIT_TAG:?}" || return 1;
  
  else
  
    gitact0 "${FOLDER:?}" "${URL:?}" || return 1;
  
  fi

  pbottom "GETTING ${PRINTNAME:?}" || { unset_all; return 1; }

fi

# ----------------------------------------------------------------------------
# ---------------------------- ROMAN_REAL ------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_COSMOLIKE_ROMAN_REAL_CODE}" ]; then 
  
  # Name to be printed on this shell script messages
  PRINTNAME="ROMAN_REAL"

  ptop "GETTING ${PRINTNAME:?}" || { unset_all; return 1; }

  FOLDER="${ROMAN_REAL_NAME:-"roman_real"}"

  URL="${ROMAN_REAL_URL:-"https://github.com/CosmoLike/cocoa_roman_real.git"}"
  URL=$(devurl "${URL:?}")

  if [ -n "${ROMAN_REAL_GIT_COMMIT:-}" ]; then
  
    gitact0 "${FOLDER:?}" "${URL:?}" || return 1;
  
    gitact2 "${FOLDER:?}" "${ROMAN_REAL_GIT_COMMIT:?}" || return 1;
  
  elif [ -n "${ROMAN_REAL_GIT_BRANCH:-}" ]; then 
  
    gitact1 "${FOLDER:?}" "${URL:?}" "${ROMAN_REAL_GIT_BRANCH:?}" || return 1;
  
  elif [ -n "${ROMAN_REAL_GIT_TAG:-}" ]; then 
  
    gitact0 "${FOLDER:?}" "${URL:?}" || return 1;
  
    gitact3 "${FOLDER:?}" "${ROMAN_REAL_GIT_TAG:?}" || return 1;
  
  else
  
    gitact0 "${FOLDER:?}" "${URL:?}" || return 1;
  
  fi

  pbottom "GETTING ${PRINTNAME:?}" || { unset_all; return 1; }

fi

# ----------------------------------------------------------------------------
# ----------------------------- DES x Planck ---------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_COSMOLIKE_DESXPLANCK_CODE}" ]; then 
  
  # Name to be printed on this shell script messages
  PRINTNAME="DESXPLANCK"

  ptop "GETTING ${PRINTNAME:?}" || { unset_all; return 1; }

  FOLDER="${DESXPLANCK_GIT_NAME:-"desy1xplanck"}"

  URL="${DESXPLANCK_URL:-"https://github.com/CosmoLike/cocoa_desy1xplanck.git"}"
  URL=$(devurl "${URL:?}")

  if [ -n "${DESXPLANCK_GIT_COMMIT:-}" ]; then

    gitact0 "${FOLDER:?}" "${URL:?}" || return 1;
    
    gitact2 "${FOLDER:?}" "${DESXPLANCK_GIT_COMMIT:?}"  || return 1;

  elif [ -n "${DESXPLANCK_GIT_BRANCH:-}" ]; then 

    gitact1 "${FOLDER:?}" "${URL:?}" "${DESXPLANCK_GIT_BRANCH:?}" || return 1;

  elif [ -n "${DESXPLANCK_GIT_TAG:-}" ]; then 

    gitact0 "${FOLDER:?}" "${URL:?}" || return 1;
    
    gitact3 "${FOLDER:?}" "${DESXPLANCK_GIT_TAG:?}" || return 1;

  else

    gitact0 "${FOLDER:?}" "${URL:?}" || return 1;

  fi

  pbottom "GETTING ${PRINTNAME:?}" || { unset_all; return 1; }

fi

# ----------------------------------------------------------------------------
# ----------------------------- ROMAN KL -------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_COSMOLIKE_ROMAN_KL_CODE}" ]; then 
  
  # Name to be printed on this shell script messages
  PRINTNAME="ROMAN_KL"

  ptop "GETTING ${PRINTNAME:?}" || { unset_all; return 1; }

  FOLDER="${ROMAN_KL_NAME:-"roman_kl"}"

  URL="${ROMAN_KL_URL:-"https://github.com/CosmoLike/cocoa_roman_kl.git"}"
  URL=$(devurl "${URL:?}")

  if [ -n "${ROMAN_KL_GIT_COMMIT:-}" ]; then

    gitact0 "${FOLDER:?}" "${URL:?}" || return 1;
    
    gitact2 "${FOLDER:?}" "${ROMAN_KL_GIT_COMMIT:?}"  || return 1;

  elif [ -n "${ROMAN_KL_GIT_BRANCH:-}" ]; then 

    gitact1 "${FOLDER:?}" "${URL:?}" "${ROMAN_KL_GIT_BRANCH:?}" || return 1;

  elif [ -n "${ROMAN_KL_GIT_TAG:-}" ]; then 

    gitact0 "${FOLDER:?}" "${URL:?}" || return 1;
    
    gitact3 "${FOLDER:?}" "${ROMAN_KL_GIT_TAG:?}" || return 1;

  else

    gitact0 "${FOLDER:?}" "${URL:?}" || return 1;

  fi

  pbottom "GETTING ${PRINTNAME:?}" || { unset_all; return 1; }

fi

# ----------------------------------------------------------------------------

unset_all || return 1;

#-----------------------------------------------------------------------------

return 55; # why this odd number? Setup_cocoa will cache this installation only
           #   if this script runs entirely. What if the user close the terminal 
           #   or the system shuts down in the middle of a git clone?  
           #   In this case, PACKDIR would exists, but it is corrupted


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
