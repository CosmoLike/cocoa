#!/bin/bash
# ------------------------------------------------------------------------------
# unxv_simons_observatory.sh: unpack the Simons Observatory likelihood code
# (mflike) into external_modules/data.
#
# Sourced, never executed (only `source` keeps its environment changes).
# setup_cocoa.sh runs it during a full installation, and sourcing it on its
# own redoes just this package (see the README FAQ on installing a single
# package). Datasets are downloaded or unpacked at setup time only, never at
# compile or run time.
#
# Skipped (return 99) when IGNORE_SIMONS_OBSERVATORY_CMB_DATA is set. Reruns
# are safe: existing results are kept unless
# OVERWRITE_EXISTING_SIMONS_OBSERVATORY_CMB_DATA is set. On failure error()
# names the failing step and unsets everything defined here; on success the
# script returns 55, the value the runner caches so a finished step is not
# repeated.
# ------------------------------------------------------------------------------
if [ -n "${IGNORE_SIMONS_OBSERVATORY_CMB_DATA:-}" ]; then
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
  unset -v TMP EDATAF FOLDER URL VER PACKDIR FILE PRINTNAME
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

unset_env_vars || return 1

# E = EXTERNAL, DATA, F=FODLER
EDATAF="${ROOTDIR:?}/external_modules/data"

FOLDER="simons_observatory"

URL=${SO_DATA_URL:-"https://portal.nersc.gov/cfs/sobs/users/MFLike_data"}

VER="${SO_DATA_VERSION:-"v0.8"}"

# PACK = PACKAGE, DIR = DIRECTORY
PACKDIR="${EDATAF:?}/${FOLDER:?}"

ptop "SETUP/UNXV SIMONS OBSERVATORY DATA" || { unset_all; return 1; }

# ---------------------------------------------------------------------------
# in case this script is called twice
# ---------------------------------------------------------------------------
if [ -n "${OVERWRITE_EXISTING_SIMONS_OBSERVATORY_CMB_DATA:-}" ]; then

  rm -rf "${PACKDIR:?}"
  if [ -n "${REDOWNLOAD_EXISTING_SIMONS_OBSERVATORY_CMB_DATA:-}" ]; then
  
    for x in $(echo "${VER:?}")
    do
      FILE="${x}.tar.gz"
      rm -f "${EDATAF:?}/${FILE:?}"
    done
  
  fi

fi

if [ ! -d "${PACKDIR:?}" ]; then
  
  mkdir -p "${PACKDIR:?}" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC20:?}"; return 1; }
  
  cdfolder "${EDATAF:?}" || return 1;

  # note: users can download multiple versions (reproduce existing work)
  # note: For example, SO_DATA_VERSION="v0.7.1 v0.8"
  for x in $(echo "${VER:?}")
  do
    FILE="${x}.tar.gz"

    if [ ! -e "${FILE:?}" ]; then
    
      "${WGET:?}" "${URL}/${FILE:?}" -q --show-progress --no-check-certificate \
        --progress=bar:force:noscroll --timeout=30 --tries=2 --waitretry=0 \
        --retry-connrefused --read-timeout=30 || { error "${EC24:?}"; return 1; }
    
    fi

    TMP=$(tar -tf "${FILE:?}" | head -1 | cut -f1 -d"/")
    
    tar -zxvf "${FILE:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC25:?}"; return 1; }
    
    mv "${TMP:?}" "${PACKDIR:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "MV SO DATA"; return 1; }
  done

fi

pbottom "SETUP/UNXV SIMONS OBSERVATORY DATA" || { unset_all; return 1; }

cdfolder "${ROOTDIR}" || return 1;

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
