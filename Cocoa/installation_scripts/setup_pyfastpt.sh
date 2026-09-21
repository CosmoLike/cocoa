#!/bin/bash
# ------------------------------------------------------------------------------
# setup_pyfastpt.sh: download FAST-PT (perturbation-theory integrals in
# python).
#
# Sourced by setup_cocoa.sh, never executed directly. Setup scripts MAY use
# the internet (git clone, wget, pip index); the matching compile step never
# does, so anything that must be downloaded happens here.
#
# Skipped (return 99) when IGNORE_FASTPT_CODE is set. Reruns are safe:
# existing results are kept unless OVERWRITE_EXISTING_FASTPT_CODE is set. On
# failure error() names the failing step and unsets everything defined here;
# on success the script returns 55, the value the runner caches so a finished
# step is not repeated.
# ------------------------------------------------------------------------------
if [ -n "${IGNORE_FASTPT_CODE:-}" ]; then
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
  unset -v URL CCIL ECODEF FOLDER PACKDIR CHANGES TFOLDER TFILE TFILEP AL PRINTNAME
  cdroot || return 1;
}

# unset_env_funcs: forget every helper function defined below (functions leak
# into the user's shell exactly like variables).
unset_env_funcs () {
  unset -f cdfolder cpfolder error
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

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

unset_env_vars || return 1

CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

# ---------------------------------------------------------------------------
# E = EXTERNAL, CODE, F=FODLER
ECODEF="${ROOTDIR:?}/external_modules/code"

URL="${FASTPT_URL:-"https://github.com/jablazek/FAST-PT.git"}"

FOLDER="${FASTPT_NAME:-"FAST-PT"}"

PACKDIR="${ECODEF:?}/${FOLDER:?}"

# Name to be printed on this shell script messages
PRINTNAME="PyFAST-PT"

ptop "INSTALLING ${PRINTNAME:?}" || { unset_all; return 1; }

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# FAST-PT
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# In case this script is called twice ----------------------------------------
# ----------------------------------------------------------------------------
if [ -n "${OVERWRITE_EXISTING_FASTPT_CODE:-}" ]; then

  rm -rf "${PACKDIR:?}";

fi

if [[ ! -d "${PACKDIR:?}" ]]; then
  
  # --------------------------------------------------------------------------
  # Clone from original FAST-PT repo -----------------------------------------
  # --------------------------------------------------------------------------
  cdfolder "${ECODEF:?}" || return 1;

  "${GIT:?}" clone "${URL:?}" --depth ${GIT_CLONE_MAXIMUM_DEPTH:-1000} \
    --recursive --no-single-branch "${FOLDER:?}" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }
  
  cdfolder "${PACKDIR}" || return 1;

  if [[ -n "${FASTPT_GIT_COMMIT:-}" ||
        -n "${FASTPT_GIT_BRANCH:-}" ||
        -n "${FASTPT_GIT_TAG:-}" ]]; then
    if [ "$("${GIT:?}" rev-parse --is-shallow-repository)" = "true" ]; then
      "${GIT:?}" fetch --unshallow --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    else
      "${GIT:?}" fetch --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi
  fi

  if [ -n "${FASTPT_GIT_COMMIT:-}" ]; then
    "${GIT:?}" checkout "${FASTPT_GIT_COMMIT:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  elif [ -n "${FASTPT_GIT_BRANCH:-}" ]; then
    "${GIT:?}" checkout -b "${FASTPT_GIT_BRANCH:?}" "origin/${FASTPT_GIT_BRANCH:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  elif [ -n "${FASTPT_GIT_TAG:-}" ]; then
    "${GIT:?}" checkout "tags/${FASTPT_GIT_TAG:?}" -b "${FASTPT_GIT_TAG:?}TMP" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  fi
    
fi

unset -v URL FOLDER PACKDIR

URL="${FASTPT_WRAPPER_URL:-"https://github.com/CosmoLike/fastpt.git"}"

FOLDER="${FASTPT_WRAPPER_NAME:-"PyFAST-PT"}"

PACKDIR="${ECODEF:?}/${FOLDER:?}"

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# FAST-PT cobaya theory wrapper
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -n "${OVERWRITE_EXISTING_FASTPT_CODE:-}" ]; then

  rm -rf "${PACKDIR:?}";

fi

if [[ ! -d "${PACKDIR:?}" ]]; then
  
  # --------------------------------------------------------------------------
  # Clone from original FAST-PT theory wrapper repo --------------------------
  # --------------------------------------------------------------------------
  cdfolder "${ECODEF:?}" || return 1;

  "${GIT:?}" clone "${URL:?}" --depth ${GIT_CLONE_MAXIMUM_DEPTH:-1000} \
    --recursive "${FOLDER:?}" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }
  
  cdfolder "${PACKDIR}" || return 1;

  if [[ -n "${FASTPT_WRAPPER_GIT_COMMIT:-}" ||
        -n "${FASTPT_WRAPPER_GIT_BRANCH:-}" ||
        -n "${FASTPT_WRAPPER_GIT_TAG:-}" ]]; then
    if [ "$("${GIT:?}" rev-parse --is-shallow-repository)" = "true" ]; then
      "${GIT:?}" fetch --unshallow --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    else
      "${GIT:?}" fetch --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi
  fi

  if [ -n "${FASTPT_WRAPPER_GIT_COMMIT:-}" ]; then
    "${GIT:?}" checkout "${FASTPT_WRAPPER_GIT_COMMIT:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  elif [ -n "${FASTPT_WRAPPER_GIT_BRANCH:-}" ]; then
    "${GIT:?}" checkout -b "${FASTPT_WRAPPER_GIT_BRANCH:?}" "origin/${FASTPT_WRAPPER_GIT_BRANCH:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  elif [ -n "${FASTPT_WRAPPER_GIT_TAG:-}" ]; then
    "${GIT:?}" checkout "tags/${FASTPT_WRAPPER_GIT_TAG:?}" -b "${FASTPT_WRAPPER_GIT_TAG:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  fi

fi

cdfolder "${ROOTDIR:?}" || return 1;

pbottom "INSTALLING ${PRINTNAME:?}" || { unset_all; return 1; }

# ----------------------------------------------------------------------------

unset_all || return 1

#-------------------------------------------------------------------------------

return 55; # why this odd number? Setup_cocoa will cache this installation only
           #   if this script runs entirely. What if the user close the terminal 
           #   or the system shuts down in the middle of a git clone?  
           #   In this case, PACKDIR would exists, but it is corrupted

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
