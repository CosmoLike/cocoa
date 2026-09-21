#!/bin/bash
# ------------------------------------------------------------------------------
# compile_simons_observatory.sh: build and install the Simons Observatory
# likelihood code (mflike).
#
# Sourced by compile_cocoa.sh, never executed directly. Compile scripts NEVER
# use the internet (they may run on compute nodes without network access):
# they only build what the matching setup script already downloaded.
#
# Skipped (return 99) when IGNORE_SIMONS_OBSERVATORY_LIKELIHOOD_CODE is set.
# On failure error() names the failing step and unsets everything defined
# here; on success the script returns 55, the value the runner caches so a
# finished step is not repeated.
# ------------------------------------------------------------------------------
if [ -n "${IGNORE_SIMONS_OBSERVATORY_LIKELIHOOD_CODE:-}" ]; then
  return 99
fi

if [ -z "${ROOTDIR:-}" ]; then
  source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
fi
  
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

PLIB="${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages"

# ----------------------------------------------------------------------------  
# ----------------------------------------------------------------------------

ptop "COMPILING SIMONS OBSERVATORY SYSLIBRARY" || { unset_all; return 1; }

PACKDIR="${ECODEF:?}/${SO_SYSLIB_NAME:-"syslibrary"}"

cdfolder "${PACKDIR}" || return 1;

# ---------------------------------------------------------------------------- 
# cleaning any previous compilation
# ----------------------------------------------------------------------------
rm -rf "${PACKDIR:?}/build/"
rm -rf "${PACKDIR:?}/syslibrary.egg-info/"  
rm -rf  "${PLIB:?}/syslibrary"
rm -rf  "${PLIB:?}/syslibrary"-*
# ----------------------------------------------------------------------------

#prevent all compile_XXX.sh from using the internet (run @compute nodes)
#FROM: https://github.com/pypa/pip/issues/12050
#That is why we use --no-dependencies --no-index --no-build-isolation
(
  env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" \
    ${PIP3:?} install . \
      --prefix="${ROOTDIR:?}/.local" \
      --no-index \
      --no-deps \
      --no-build-isolation \
) >>${OUT1:?} 2>>${OUT2:?} || { error "${EC13:?}"; return 1; }
  
pbottom "COMPILING SIMONS OBSERVATORY SYSLIBRARY" || { unset_all; return 1; }

cdfolder "${ROOTDIR}" || return 1;

# ----------------------------------------------------------------------------  
# ----------------------------------------------------------------------------

ptop "COMPILING SIMONS OBSERVATORY MKLIKE" || { unset_all; return 1; }

PACKDIR="${ECODEF:?}/${SO_MFLIKE_NAME:-"mflike"}"

cdfolder "${PACKDIR}" || return 1;

# ----------------------------------------------------------------------------
# cleaning any previous compilation
# ----------------------------------------------------------------------------
rm -rf "${PACKDIR:?}/build/"
rm -rf "${PACKDIR:?}/mflike.egg-info/"
rm -rf  "${PLIB:?}/mflike"
rm -rf  "${PLIB:?}/mflike"-*
# ----------------------------------------------------------------------------

(
  env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" \
    ${PIP3:?} install "${PACKDIR:?}" \
      --prefix="${ROOTDIR:?}/.local" \
      --no-index \
      --no-dependencies \
      --no-build-isolation
) >>${OUT1:?} 2>>${OUT2:?} || { error "${EC13:?}"; return 1; }

pbottom "COMPILING SIMONS OBSERVATORY MKLIKE" || { unset_all; return 1; }

cdfolder "${ROOTDIR}" || return 1;

# ----------------------------------------------------------------------------

unset_all || return 1  

#-------------------------------------------------------------------------------

return 55; # why this odd number? compile_cocoa will cache this compilation only
           #   if this script runs entirely.

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
