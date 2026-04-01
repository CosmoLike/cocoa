#!/bin/bash

# IGNORE_EMULTRF_CODE EMULTRF_URL OVERWRITE_EXISTING_EMULTRF_CODE EMULTRF_GIT_COMMIT

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -n "${IGNORE_EMULTRF_CODE:-}" ]; then
  return 99
fi

if [ -z "${ROOTDIR:-}" ]; then
  source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
fi

# parenthesis = run in a subshell
( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

unset_env_vars () {
  unset -v URL CCIL ECODEF FOLDER PACKDIR  
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

# E = EXTERNAL, CODE, F=FODLER
ECODEF="${ROOTDIR:?}/external_modules/code"

# ----------------------------------------------------------------------------
# In case this script is called twice ----------------------------------------
# ----------------------------------------------------------------------------

if [ ! -d "${ECODEF:?}/emulators" ]; then
  mkdir -p "${ECODEF:?}/emulators" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC20:?}"; return 1; }
fi

URL="${EMULTRF_URL:-"https://github.com/CosmoLike/emulators_code.git"}"

FOLDER="${EMULTRF_NAME:-"emultrf"}"

PACKDIR="${ECODEF:?}/emulators/${FOLDER:?}"

ptop "INSTALLING EMULATOR CMB TRF THEORY (COBAYA)" || { unset_all; return 1; }

# ----------------------------------------------------------------------------
# In case this script is called twice ----------------------------------------
# ----------------------------------------------------------------------------
if [ -n "${OVERWRITE_EXISTING_EMULTRF_CODE:-}" ]; then

  rm -rf "${PACKDIR:?}"

fi

if [ ! -d "${PACKDIR:?}" ]; then

  cdfolder "${ECODEF:?}/emulators" || { unset_all; return 1; }

  "${GIT:?}" clone "${URL:?}" \
    --depth "${GIT_CLONE_MAXIMUM_DEPTH:-1000}" \
    --recurse-submodules \
    "${PACKDIR:?}" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }

  cdfolder "${PACKDIR:?}" || { unset_all; return 1; }

  if [[ -n "${EMULTRF_GIT_COMMIT:-}" || 
        -n "${EMULTRF_GIT_BRANCH:-}" || 
        -n "${EMULTRF_GIT_TAG:-}" ]]; then

    if [ "$("${GIT:?}" rev-parse --is-shallow-repository)" = "true" ]; then
      "${GIT:?}" fetch --unshallow --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    else
      "${GIT:?}" fetch --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi

  fi

  if [ -n "${EMULTRF_GIT_COMMIT:-}" ]; then

    "${GIT:?}" checkout "${EMULTRF_GIT_COMMIT:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }

  elif [ -n "${EMULTRF_GIT_BRANCH}" ]; then
    
    "${GIT:?}" checkout "${EMULTRF_GIT_BRANCH:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }

  elif [ -n "${EMULTRF_GIT_TAG:-}" ]; then 

    "${GIT:?}" checkout tags/${EMULTRF_GIT_TAG} -b ${EMULTRF_GIT_TAG} \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  
  fi  

fi

cdfolder "${ROOTDIR:?}" || { unset_all; return 1; }

pbottom "INSTALLING EMULATOR CMB TRF THEORY (COBAYA)" || { unset_all; return 1; }

# ---------------------------------------------------------------------------

unset_all || return 1
  
#-------------------------------------------------------------------------------

return 55; # why this odd number? Setup_cocoa will cache this installation only
           #   if this script runs entirely. What if the user close the terminal 
           #   or the system shuts down in the middle of a git clone?  
           #   In this case, PACKDIR would exists, but it is corrupted

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------