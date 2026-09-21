#!/bin/bash
# ------------------------------------------------------------------------------
# setup_simons_observatory.sh: download the Simons Observatory likelihood code
# (mflike).
#
# Sourced by setup_cocoa.sh, never executed directly. Setup scripts MAY use
# the internet (git clone, wget, pip index); the matching compile step never
# does, so anything that must be downloaded happens here.
#
# Skipped (return 99) when IGNORE_SIMONS_OBSERVATORY_LIKELIHOOD_CODE is set.
# Reruns are safe: existing results are kept unless
# OVERWRITE_EXISTING_SIMONS_OBSERVATORY_CODE is set. On failure error() names
# the failing step and unsets everything defined here; on success the script
# returns 55, the value the runner caches so a finished step is not repeated.
# ------------------------------------------------------------------------------
if [ -n "${IGNORE_SIMONS_OBSERVATORY_LIKELIHOOD_CODE:-}" ]; then
  return 99
fi

if [ -z "${ROOTDIR:-}" ]; then
  source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
fi

# parenthesis = run in a subshell
( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

# unset_env_vars: forget every variable the body defines. Sourced scripts
# share the caller's shell, so anything not unset here would leak into the
# user's environment. cdroot returns to the directory the user launched from.
unset_env_vars () {
  unset -v COB CCCOB COBLIKE URL FOLDER PACKDIR PRINTNAME ECODEF URL COBTH
  cdroot || return 1;
}

# unset_env_funcs: forget every helper function defined below (functions leak
# into the user's shell exactly like variables).
unset_env_funcs () {
  unset -f cdfolder error
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

unset_env_vars || return 1

# E = EXTERNAL, CODE, F=FODLER
ECODEF="${ROOTDIR:?}/external_modules/code"

CCIL="${ROOTDIR:?}/../cocoa_installation_libraries" # IL = installation lib

COB="${ROOTDIR:?}/cobaya"        # COB = Cobaya

CCCOB="${CCIL:?}/cobaya_changes" # CC = CoCoA, COB = Cobaya (Cocoa Cobaya)

COBLIKE="cobaya/likelihoods"     # COB = Cobaya, LIKE = likelihoods

COBTH="cobaya/theories"          # COB = Cobaya, TH = theories

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

ptop "SETUP SIMONS OBSERVATORY SYSLIBRARY" || { unset_all; return 1; }

URL="${SO_SYSLIB_URL:-"https://github.com/simonsobs/syslibrary.git"}"
  
FOLDER="${SO_SYSLIB_NAME:-"syslibrary"}"

PACKDIR="${ECODEF:?}/${FOLDER:?}"

if [[ -n "${OVERWRITE_EXISTING_SIMONS_OBSERVATORY_CODE:-}" ]]; then

  rm -rf "${PACKDIR:?}"

fi

if [[ ! -d "${PACKDIR:?}" ]]; then

  cdfolder "${ECODEF}" || return 1;

  "${GIT:?}" clone "${URL:?}" --depth ${GIT_CLONE_MAXIMUM_DEPTH:-1000} \
    --recursive --no-single-branch "${FOLDER:?}" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }

  cdfolder "${PACKDIR}" || return 1;

  if [[ -n "${SO_SYSLIB_GIT_COMMIT:-}" ||
        -n "${SO_SYSLIB_GIT_BRANCH:-}" ||
        -n "${SO_SYSLIB_GIT_TAG:-}" ]]; then
    if [ "$("${GIT:?}" rev-parse --is-shallow-repository)" = "true" ]; then
      "${GIT:?}" fetch --unshallow --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    else
      "${GIT:?}" fetch --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi
  fi

  if [ -n "${SO_SYSLIB_GIT_COMMIT:-}" ]; then
    "${GIT:?}" checkout "${SO_SYSLIB_GIT_COMMIT:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  elif [ -n "${SO_SYSLIB_GIT_BRANCH:-}" ]; then
    "${GIT:?}" checkout -b "${SO_SYSLIB_GIT_BRANCH:?}" "origin/${SO_SYSLIB_GIT_BRANCH:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  elif [ -n "${SO_SYSLIB_GIT_TAG:-}" ]; then
    "${GIT:?}" checkout "tags/${SO_SYSLIB_GIT_TAG:?}" -b "${SO_SYSLIB_GIT_TAG:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  fi

fi

pbottom "SETUP SIMONS OBSERVATORY SYSLIBRARY" || { unset_all; return 1; }

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

ptop "SETUP SIMONS OBSERVATORY MFLIKE" || { unset_all; return 1; }

URL="${SO_MFLIKE_URL:-"https://github.com/simonsobs/LAT_MFLike.git"}"

FOLDER=${SO_MFLIKE_NAME:-"mflike"}

PACKDIR="${ECODEF:?}/${FOLDER:?}"

if [[ -n "${OVERWRITE_EXISTING_SIMONS_OBSERVATORY_CODE:-}" ]]; then    

  rm -rf "${PACKDIR:?}"

fi

if [[ ! -d "${PACKDIR:?}" ]]; then

  cdfolder "${ECODEF}" || return 1;

  "${GIT:?}" clone "${URL:?}" --depth ${GIT_CLONE_MAXIMUM_DEPTH:-1000} \
    --recursive --no-single-branch "${FOLDER:?}" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }

  cdfolder "${PACKDIR}" || return 1;

  if [ -n "${SO_MFLIKE_GIT_COMMIT:-}" ]; then
    if [ "$("${GIT:?}" rev-parse --is-shallow-repository)" = "true" ]; then
      "${GIT:?}" fetch --unshallow --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    else
      "${GIT:?}" fetch --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi
    
    "${GIT:?}" checkout "${SO_MFLIKE_GIT_COMMIT:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  fi

  # PATCH MKFLIKE
  cdfolder "${PACKDIR:?}/mflike"|| return 1;

  cp "${CCCOB:?}/${COBLIKE:?}/mflike/mflike.patch" "${PACKDIR:?}/mflike" \
    2>>"/dev/null" || { error "CP FILE mflike.patch"; return 1; }

  patch -u "mflike.py" -i "mflike.patch" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC17:?} (mflike.patch)"; return 1; }

fi

pbottom "SETUP SIMONS OBSERVATORY MFLIKE" || { unset_all; return 1; }

cdfolder "${ROOTDIR}" || return 1;

#-------------------------------------------------------------------------------

unset_all || return 1;

#-------------------------------------------------------------------------------

return 55; # why this odd number? Setup_cocoa will cache this installation only
           #   if this script runs entirely. What if the user close the terminal 
           #   or the system shuts down in the middle of a git clone?  
           #   In this case, PACKDIR would exists, but it is corrupted


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
