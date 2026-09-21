#!/bin/bash
# ------------------------------------------------------------------------------
# compile_camb.sh: build and install CAMB (the Boltzmann code computing CMB
# and matter power spectra).
#
# Sourced by compile_cocoa.sh, never executed directly. Compile scripts NEVER
# use the internet (they may run on compute nodes without network access):
# they only build what the matching setup script already downloaded.
#
# Skipped (return 99) when IGNORE_CAMB_CODE is set. On failure error() names
# the failing step and unsets everything defined here; on success the script
# returns 55, the value the runner caches so a finished step is not repeated.
# ------------------------------------------------------------------------------
if [ -n "${IGNORE_CAMB_CODE:-}" ]; then
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
  unset -v ECODEF FOLDER PACKDIR PRINTNAME PLIB
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
  cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER ${1}"; return 1; }
}

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

unset_env_vars || return 1

# ----------------------------------------------------------------------------

# E = EXTERNAL, CODE, F=FODLER
ECODEF="${ROOTDIR:?}/external_modules/code"

FOLDER="${CAMB_NAME:-"CAMB"}"

PACKDIR="${ECODEF:?}/${FOLDER:?}"

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

ptop "COMPILING CAMB" || { unset_all; return 1; }

cdfolder "${PACKDIR:?}" || return 1;

# ---------------------------------------------------------------------------- 
# cleaning any previous compilation
# ----------------------------------------------------------------------------
rm -rf "${PACKDIR:?}/build/"
rm -rf "${PACKDIR:?}/camb/__pycache__/"
rm -f  "${PACKDIR:?}/camb/camblib.so"
rm -rf "${PACKDIR:?}/forutils/Releaselib/"
PLIB="${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages"
rm -rf  "${PLIB:?}/camb"
rm -rf  "${PLIB:?}/camb"-*
# ---------------------------------------------------------------------------

if [ -z "${IGNORE_COSMOREC_CODE:-}" ] && [ -n "${IGNORE_HYREC_CODE:-}" ]; then

  (
    RECOMBINATION_FILES="recfast cosmorec" \
    COMPILER="${FORTRAN_COMPILER:?}" F90C="${FORTRAN_COMPILER:?}" \
    "${PYTHON3:?}" setup.py build \
  ) >>${OUT1:?} 2>>${OUT2:?} || { error "${EC4:?}"; return 1; }

elif [ -n "${IGNORE_COSMOREC_CODE:-}" ] && [ -z "${IGNORE_HYREC_CODE:-}" ]; then

  (
    RECOMBINATION_FILES="recfast hyrec" \
    COMPILER="${FORTRAN_COMPILER:?}" F90C="${FORTRAN_COMPILER:?}" \
    "${PYTHON3:?}" setup.py build
  ) >>${OUT1:?} 2>>${OUT2:?} || { error "${EC4:?}"; return 1; }

elif [ -z "${IGNORE_COSMOREC_CODE:-}" ] && [ -z "${IGNORE_HYREC_CODE:-}" ]; then

  (
    RECOMBINATION_FILES="recfast cosmorec hyrec" \
    COMPILER="${FORTRAN_COMPILER:?}" F90C="${FORTRAN_COMPILER:?}" \
    "${PYTHON3:?}" setup.py build
  ) >>${OUT1:?} 2>>${OUT2:?} || { error "${EC4:?}"; return 1; }

else

  (
    COMPILER="${FORTRAN_COMPILER:?}" F90C="${FORTRAN_COMPILER:?}" \
    "${PYTHON3:?}" setup.py build
  ) >>${OUT1:?} 2>>${OUT2:?} || { error "${EC4:?}"; return 1; }

fi

#prevent all compile_XXX.sh from using the internet (run @compute nodes)
#FROM: https://github.com/pypa/pip/issues/12050
#That is why we use --no-dependencies --no-index --no-build-isolation
(
  env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" \
    ${PIP3:?} install "${PACKDIR:?}" \
      --no-dependencies \
      --prefix="${ROOTDIR:?}/.local" \
      --no-index \
      --no-build-isolation 
) >>${OUT1:?} 2>>${OUT2:?} || { error "${EC3:?}"; return 1; }

pbottom "COMPILING CAMB" || { unset_all; return 1; }

cdfolder "${ROOTDIR}" || return 1;

# ---------------------------------------------------------------------------

unset_all || return 1

#-------------------------------------------------------------------------------

return 55; # why this odd number? compile_cocoa will cache this compilation only
           #   if this script runs entirely.
           
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
