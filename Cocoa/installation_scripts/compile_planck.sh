#!/bin/bash
# ------------------------------------------------------------------------------
# compile_planck.sh: build and install the Planck 2018 clik likelihood code.
#
# Sourced by compile_cocoa.sh, never executed directly. Compile scripts NEVER
# use the internet (they may run on compute nodes without network access):
# they only build what the matching setup script already downloaded.
#
# Skipped (return 99) when IGNORE_PLANCK_LIKELIHOOD_CODE is set. On failure
# error() names the failing step and unsets everything defined here; on
# success the script returns 55, the value the runner caches so a finished
# step is not repeated.
# ------------------------------------------------------------------------------
if [ -n "${IGNORE_PLANCK_LIKELIHOOD_CODE:-}" ]; then
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
  unset -v ECPCF CLIK_LAPACK_LIBS CLIK_CFITSIO_LIBS PRINTNAME PACKDIR
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

# ---------------------------------------------------------------------------
# Name to be printed on this shell script messages

PRINTNAME="PLANCK LIKELIHOOD"

ECPCF="external_modules/code/planck/code"

if [ -z "${IGNORE_C_CFITSIO_INSTALLATION:-}" ]; then

  CLIK_CFITSIO_LIBS="${ROOTDIR:?}/.local/lib"

else

  if [ -z "${GLOBAL_PACKAGES_LOCATION:-}" ]; then
    pfail "GLOBAL_PACKAGES_LOCATION"; { unset_all; return 1; }
  fi

  CLIK_CFITSIO_LIBS="${GLOBAL_PACKAGES_LOCATION:?}"

fi

if [ -z "${IGNORE_FORTRAN_INSTALLATION:-}" ]; then

  CLIK_LAPACK_LIBS="${ROOTDIR:?}/.local"

else

  if [ -z "${GLOBAL_PACKAGES_LOCATION:-}" ]; then
    pfail "GLOBAL_PACKAGES_LOCATION"; { unset_all; return 1; }
  fi

  CLIK_LAPACK_LIBS="${GLOBAL_PACKAGES_LOCATION:?}"

fi
  
if [ -z "${USE_SPT_CLIK_PLANCK:-}" ]; then 

  PACKDIR="${ROOTDIR:?}/${ECPCF:?}/plc_3.0/plc-3.1"

else

  PACKDIR="${ROOTDIR:?}/${ECPCF:?}/spt_clik"

fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
ptop "COMPILING ${PRINTNAME:?}" || { unset_all; return 1; }

cdfolder "${PACKDIR:?}" || return 1;

# ----------------------------------------------------------------------------
# cleaning any previous compilation
# ----------------------------------------------------------------------------
rm -f  "${ROOTDIR:?}/.local"/bin/clik*
rm -f  "${ROOTDIR:?}/.local"/lib/libclik_f90.so
rm -f  "${ROOTDIR:?}/.local"/lib/libclik.so
rm -rf "${ROOTDIR:?}/.local"/lib/python/site-packages/clik
rm -rf "${ROOTDIR:?}/.local"/share/clik
rm -f  "${ROOTDIR:?}/.local"/include/clik*
rm -f  "${PACKDIR:?}/".lock-waf_*
"${PYTHON3:?}" waf distclean \
  >>${OUT1:?} 2>>${OUT2:?} || { error "${EC18:?}"; return 1; }
# ---------------------------------------------------------------------------

case "$(uname -s)" in
  Linux)
    (
      export LD_LIBRARY_PATH=${CONDA_PREFIX:?}/lib:$LD_LIBRARY_PATH && \
      export LD_LIBRARY_PATH=${ROOTDIR:?}/.local/lib:$LD_LIBRARY_PATH && \
      FC="${FORTRAN_COMPILER:?}" CC="${C_COMPILER:?}" CXX="${CXX_COMPILER:?}" \
      "${PYTHON3:?}" waf configure \
        --gcc \
        --gfortran \
        --cfitsio_islocal \
        --prefix="${ROOTDIR:?}/.local" \
        --lapack_prefix="${CLIK_LAPACK_LIBS:?}" \
        --cfitsio_lib="${CLIK_CFITSIO_LIBS:?}" \
        --python="${PYTHON3:?}" 
    ) >>${OUT1:?} 2>>${OUT2:?} || { error "${EC5:?}"; return 1; }
    ;;
  Darwin)
    (
      export LD_LIBRARY_PATH=${CONDA_PREFIX:?}/lib:$LD_LIBRARY_PATH && \
      export LD_LIBRARY_PATH=${ROOTDIR:?}/.local/lib:$LD_LIBRARY_PATH && \
      FC="${FORTRAN_COMPILER:?}" CC="${C_COMPILER:?}" CXX="${CXX_COMPILER:?}" \
      AR="/usr/bin/ar" RANLIB="/usr/bin/ranlib" CFLAGS="-O2 -fPIC -std=gnu89" \
      "${PYTHON3:?}" waf configure \
        --cfitsio_islocal \
        --prefix="${ROOTDIR:?}/.local" \
        --lapack_prefix="${CLIK_LAPACK_LIBS:?}" \
        --cfitsio_lib="${CLIK_CFITSIO_LIBS:?}" \
        --python="${PYTHON3:?}" 
    ) >>${OUT1:?} 2>>${OUT2:?} || { error "${EC5:?}"; return 1; }
    ;;
  *)
    error "UNSUPPORTED OS: $(uname -s)"; return 1;
    ;;
esac

(
  export LD_LIBRARY_PATH=${CONDA_PREFIX:?}/lib:$LD_LIBRARY_PATH && \
  export LD_LIBRARY_PATH=${ROOTDIR:?}/.local/lib:$LD_LIBRARY_PATH && \
 "${PYTHON3:?}" waf install -v 
)>>${OUT1:?} 2>>${OUT2:?} || { error "${EC6:?}"; return 1; }


pbottom "COMPILING ${PRINTNAME:?}" || { unset_all; return 1; }

cdfolder "${ROOTDIR}" || return 1;

# ---------------------------------------------------------------------------

unset_all || return 1

#-------------------------------------------------------------------------------

return 55; # why this odd number? compile_cocoa will cache this compilation only
           #   if this script runs entirely.
           
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
