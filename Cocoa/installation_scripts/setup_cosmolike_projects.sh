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
  if [ -n "${OVERWRITE_EXISTING_COSMOLIKE_CODE:-}" ]; then
    rm -rf "${PACKDIR:?}"
  fi

  if [ ! -d "${PACKDIR:?}" ]; then
  
    "${GIT:?}" clone "${URL}" --depth ${GIT_CLONE_MAXIMUM_DEPTH:-1000} \
      --recursive "${NAME}" \
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

    "${GIT:?}" checkout "${TAG}" \
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
# -------------------------------- LSST-Y1 -----------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_COSMOLIKE_LSST_Y1_CODE}" ]; then 
  
  # Name to be printed on this shell script messages
  PRINTNAME="LSST_Y1"

  ptop "GETTING ${PRINTNAME:?}" || { unset_all; return 1; }

  FOLDER="${LSST_Y1_NAME:-"lsst_y1"}"

  URL="${LSST_Y1_URL:-"https://github.com/CosmoLike/cocoa_lsst_y1.git"}"

  if [ -n "${LSST_Y1_GIT_COMMIT:-}" ]; then
  
    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }
  
    gitact2 "${FOLDER:?}" "${LSST_Y1_GIT_COMMIT:?}" || { unset_all; return 1; }
  
  elif [ -n "${LSST_Y1_GIT_BRANCH:-}" ]; then 
  
    gitact1 "${FOLDER:?}" "${URL:?}" "${LSST_Y1_GIT_BRANCH:?}" || { unset_all; return 1; }
  
  elif [ -n "${LSST_Y1_GIT_TAG:-}" ]; then 
  
    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }
  
    gitact3 "${FOLDER:?}" "${LSST_Y1_GIT_TAG:?}" || { unset_all; return 1; }
  
  else
  
    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }
  
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

  if [ -n "${DES_Y3_GIT_COMMIT:-}" ]; then
  
    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }
  
    gitact2 "${FOLDER:?}" "${DES_Y3_GIT_COMMIT:?}" || { unset_all; return 1; }
  
  elif [ -n "${DES_Y3_GIT_BRANCH:-}" ]; then 
  
    gitact1 "${FOLDER:?}" "${URL:?}" "${DES_Y3_GIT_BRANCH:?}" || { unset_all; return 1; }
  
  elif [ -n "${DES_Y3_GIT_TAG:-}" ]; then 
  
    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }
  
    gitact3 "${FOLDER:?}" "${DES_Y3_GIT_TAG:?}" || { unset_all; return 1; }
  
  else
  
    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }
  
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

  if [ -n "${ROMAN_FOURIER_GIT_COMMIT:-}" ]; then
  
    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }
   
    gitact2 "${FOLDER:?}" "${ROMAN_FOURIER_GIT_COMMIT:?}" || { unset_all; return 1; }
  
  elif [ -n "${ROMAN_FOURIER_GIT_BRANCH:-}" ]; then 
  
    gitact1 "${FOLDER:?}" "${URL:?}" "${ROMAN_FOURIER_GIT_BRANCH:?}" || { unset_all; return 1; }
  
  elif [ -n "${ROMAN_FOURIER_GIT_TAG:-}" ]; then 
  
    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }
   
    gitact3 "${FOLDER:?}" "${ROMAN_FOURIER_GIT_TAG:?}" || { unset_all; return 1; }
  
  else
  
    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }
  
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

  if [ -n "${ROMAN_REAL_GIT_COMMIT:-}" ]; then
  
    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }
  
    gitact2 "${FOLDER:?}" "${ROMAN_REAL_GIT_COMMIT:?}" || { unset_all; return 1; }
  
  elif [ -n "${ROMAN_REAL_GIT_BRANCH:-}" ]; then 
  
    gitact1 "${FOLDER:?}" "${URL:?}" "${ROMAN_REAL_GIT_BRANCH:?}" || { unset_all; return 1; }
  
  elif [ -n "${ROMAN_REAL_GIT_TAG:-}" ]; then 
  
    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }
  
    gitact3 "${FOLDER:?}" "${ROMAN_REAL_GIT_TAG:?}" || { unset_all; return 1; }
  
  else
  
    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }
  
  fi

  pbottom "GETTING ${PRINTNAME:?}" || { unset_all; return 1; }

fi

# ----------------------------------------------------------------------------
# ----------------------------- DES x Planck ------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_COSMOLIKE_DESXPLANCK_CODE}" ]; then 
  
  # Name to be printed on this shell script messages
  PRINTNAME="DESXPLANCK"

  ptop "GETTING ${PRINTNAME:?}" || { unset_all; return 1; }

  FOLDER="${DESXPLANCK_NAME:-"desy1xplanck"}"

  URL="${DESXPLANCK_URL:-"https://github.com/CosmoLike/cocoa_desy1xplanck.git"}"

  if [ -n "${DESXPLANCK_GIT_COMMIT:-}" ]; then

    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }
    
    gitact2 "${FOLDER:?}" "${DESXPLANCK_GIT_COMMIT:?}"  || { unset_all; return 1; }

  elif [ -n "${DESXPLANCK_GIT_BRANCH:-}" ]; then 

    gitact1 "${FOLDER:?}" "${URL:?}" "${DESXPLANCK_GIT_BRANCH:?}" || { unset_all; return 1; }

  elif [ -n "${DESXPLANCK_GIT_TAG:-}" ]; then 

    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }
    
    gitact3 "${FOLDER:?}" "${DESXPLANCK_GIT_TAG:?}" || { unset_all; return 1; }

  else

    gitact0 "${FOLDER:?}" "${URL:?}" || { unset_all; return 1; }

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