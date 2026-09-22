#!/bin/bash
# ------------------------------------------------------------------------------
# unxv_bao.sh: unpack the BAO likelihood data files into
# external_modules/data.
#
# Sourced, never executed (only `source` keeps its environment changes).
# setup_cocoa.sh runs it during a full installation, and sourcing it on its
# own redoes just this package (see the README FAQ on installing a single
# package). Datasets are downloaded or unpacked at setup time only, never at
# compile or run time.
#
# Skipped (return 99) when IGNORE_BAO_DATA is set. Reruns are safe: existing
# results are kept unless OVERWRITE_EXISTING_BAO_DATA is set. On failure
# error() names the failing step and unsets everything defined here; on
# success the script returns 55, the value the runner caches so a finished
# step is not repeated.
# ------------------------------------------------------------------------------
if [ -n "${IGNORE_BAO_DATA:-}" ]; then
  return 99
fi

if [ -z "${ROOTDIR}" ]; then
  source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
fi

# parenthesis = run in a subshell 
( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

# unset_env_vars: forget every variable the body defines. Sourced scripts
# share the caller's shell, so anything not unset here would leak into the
# user's environment. cdroot returns to the directory the user launched from.
unset_env_vars () {
  unset -v EDATAF FOLDER FILE PACKDIR PRINTNAME
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

# --------------------------------------------------------------------------- 
# --------------------------------------------------------------------------- 
# ---------------------------------------------------------------------------

unset_env_vars || return 1

# E = EXTERNAL, DATA, F=FODLER
EDATAF="${ROOTDIR:?}/external_modules/data"

FOLDER="bao_data"

# PACK = PACKAGE, DIR = DIRECTORY
PACKDIR="${EDATAF:?}/${FOLDER:?}"

FILE="bao_data.xz"

# Name to be printed on this shell script messages
PRINTNAME=BAO

# ---------------------------------------------------------------------------

ptop "SETUP/UNXV ${PRINTNAME:?} DATA" || { unset_all; return 1; }

# ---------------------------------------------------------------------------
# note: in case script run >1x w/ previous run stoped prematurely b/c error
# ---------------------------------------------------------------------------
if [ -n "${OVERWRITE_EXISTING_BAO_DATA:-}" ]; then
  
  rm -rf "${PACKDIR:?}"

fi

if [ ! -d "${PACKDIR:?}" ]; then

  cdfolder "${EDATAF:?}" || return 1;

  tar xf "${FILE:?}" >>${OUT1:?} 2>>${OUT2:?} || { error "${EC25:?}"; return 1; }

fi

cdfolder "${ROOTDIR}" || return 1;

pbottom "SETUP/UNXV ${PRINTNAME:?} DATA" || { unset_all; return 1; }

#-------------------------------------------------------------------------------

unset_all || return 1

#-------------------------------------------------------------------------------

return 55; # why this odd number? Setup_cocoa will cache this installation only
           #   if this script runs entirely. What if the user close the terminal 
           #   or the system shuts down in the middle of a git clone?  
           #   In this case, PACKDIR would exists, but it is corrupted

 
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
