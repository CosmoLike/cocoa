#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${INSTALL_AXIE_CAMB_V2:-}" ]; then
  return 99
fi

if [ -z "${ROOTDIR:-}" ]; then
  source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
fi

# parenthesis = run in a subshell  
( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

unset_env_vars () {
  unset -v ECODEF FOLDER PACKDIR PRINTNAME PLIB
  cdroot || return 1;
}

unset_env_funcs () {
  unset -f cdfolder cpfolder error
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
  cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER ${1}"; return 1; }
}

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

unset_env_vars || return 1

# ----------------------------------------------------------------------------

# E = EXTERNAL, CODE, F=FODLER
ECODEF="${ROOTDIR:?}/external_modules/code"

FOLDER="${AXIE_CAMB_NAME:-"AXIECAMB"}"

PACKDIR="${ECODEF:?}/${FOLDER:?}"

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

ptop "COMPILING AXIECAMB" || { unset_all; return 1; }

cdfolder "${PACKDIR:?}" || { unset_all; return 1; }

# ---------------------------------------------------------------------------- 
# cleaning any previous compilation
# ----------------------------------------------------------------------------
rm -rf "${PACKDIR:?}/build/"
rm -rf "${PACKDIR:?}/camb/__pycache__/"
rm -f  "${PACKDIR:?}/camb/camblib.so"
rm -rf "${PACKDIR:?}/forutils/Releaselib/"
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

pbottom "COMPILING AXIECAMB" || { unset_all; return 1; }

cdfolder "${ROOTDIR}" || { unset_all; return 1; }

# ---------------------------------------------------------------------------

unset_all || return 1

#-------------------------------------------------------------------------------

return 55; # why this odd number? compile_cocoa will cache this compilation only
           #   if this script runs entirely.
           
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
