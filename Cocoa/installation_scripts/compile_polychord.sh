#!/bin/bash
# ------------------------------------------------------------------------------
# compile_polychord.sh: build and install the PolyChord nested sampler.
#
# Sourced by compile_cocoa.sh, never executed directly. Compile scripts NEVER
# use the internet (they may run on compute nodes without network access):
# they only build what the matching setup script already downloaded.
#
# Skipped (return 99) when IGNORE_POLYCHORD_SAMPLER_CODE is set. On failure
# error() names the failing step and unsets everything defined here; on
# success the script returns 55, the value the runner caches so a finished
# step is not repeated.
# ------------------------------------------------------------------------------
if [ -n "${IGNORE_POLYCHORD_SAMPLER_CODE:-}" ]; then
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
  unset -v ECODEF FOLDER PACKDIR PLIB
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

# ----------------------------------------------------------------------------

# E = EXTERNAL, CODE, F=FODLER
ECODEF="${ROOTDIR:?}/external_modules/code"

FOLDER=${POLY_NAME:-"PolyChordLite"}

PACKDIR="${ECODEF:?}/${FOLDER:?}"

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
ptop 'COMPILING POLYCHORD' || { unset_all; return 1; }

cdfolder "${PACKDIR}" || return 1;

# ----------------------------------------------------------------------------
# cleaning any previous compilation
# ----------------------------------------------------------------------------
make clean \
  >>${OUT1:?} 2>>${OUT2:?} || { error "${EC2:?}"; return 1; }
PLIB="${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages"
rm -rf "${PLIB:?}"/pypolychord-*
rm -rf "${PACKDIR:?}/lib"/*.a
rm -rf "${PACKDIR:?}/lib"/*.so
rm -rf "${PACKDIR:?}/build"
rm -rf "${PACKDIR:?}/dist"
rm -rf "${PACKDIR:?}/pypolychord.egg-info"
# ---------------------------------------------------------------------------

make \
  >>${OUT1:?} 2>>${OUT2:?} || { error "${EC7:?}"; return 1; }

make -j "${MNT:-1}" all \
  >>${OUT1:?} 2>>${OUT2:?} || { error "${EC7:?}"; return 1; }

make -j "${MNT:-1}" pypolychord \
  >>${OUT1:?} 2>>${OUT2:?} || { error "${EC8:?}"; return 1; }

(
  CC="${MPI_CC_COMPILER:?}" CXX="${MPI_CXX_COMPILER:?}" \
  "${PYTHON3:?}" setup.py install --prefix "${ROOTDIR:?}/.local" \
) >>${OUT1:?} 2>> ${OUT2:?} || { error "${EC9:?}"; return 1; }

pbottom 'COMPILING POLYCHORD' || { unset_all; return 1; }

cdfolder "${ROOTDIR}" || return 1;

# ----------------------------------------------------------------------------

unset_all || return 1;

#-------------------------------------------------------------------------------

return 55; # why this odd number? compile_cocoa will cache this compilation only
           #   if this script runs entirely.
           
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
