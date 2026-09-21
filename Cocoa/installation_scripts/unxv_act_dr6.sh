#!/bin/bash
# ------------------------------------------------------------------------------
# unxv_act_dr6.sh: unpack the ACT DR6 CMB lensing likelihood code into
# external_modules/data.
#
# Sourced by setup_cocoa.sh, never executed directly. Datasets are downloaded
# or unpacked at setup time only, never at compile or run time.
#
# Skipped (return 99) when IGNORE_ACTDR6_DATA is set. Reruns are safe:
# existing results are kept unless OVERWRITE_EXISTING_ACTDR6_CMB_DATA is set.
# On failure error() names the failing step and unsets everything defined
# here; on success the script returns 55, the value the runner caches so a
# finished step is not repeated.
# ------------------------------------------------------------------------------
if [ -n "${IGNORE_ACTDR6_DATA:-}" ]; then
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
  unset -v EDATAF FOLDER PACKDIR FILE URL_BASE URL PRINTNAME TMP
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

# --------------------------------------------------------------------------- 
# ---------------------------------------------------------------------------

ptop "SETUP/UNXV ACT-DR6 DATA (LENSING)" || { unset_all; return 1; }

FOLDER="act"

# PACK = PACKAGE, DIR = DIRECTORY
PACKDIR="${EDATAF:?}/${FOLDER:?}"

FILE="${ACTDR6_LENSING_DATA_FILE:-"ACT_dr6_likelihood_v1.2.tgz"}"

URL_BASE="https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr6/likelihood/data"

URL="${ACTDR6_LENSING_DATA_URL:-"${URL_BASE:?}"}/${FILE:?}"

# ---------------------------------------------------------------------------
# note: in case script run >1x w/ previous run stoped prematurely b/c error
# ---------------------------------------------------------------------------
if [ -n "${OVERWRITE_EXISTING_ACTDR6_CMB_DATA:-}" ]; then
  
  rm -rf "${PACKDIR:?}"

  if [ -n "${REDOWNLOAD_EXISTING_ACTDR6_CMB_DATA:-}" ]; then
    rm -f "${EDATAF:?}/${FILE:?}"
  fi

fi 

if [[ ! -d "${PACKDIR:?}" || ! -d "${PACKDIR:?}/lensing" ]]; then
  
  mkdir -p "${PACKDIR:?}" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC20:?}"; return 1; }
  
  mkdir -p "${PACKDIR:?}/lensing" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC20:?}"; return 1; }
  
  cdfolder "${EDATAF:?}" || return 1;

  if [ ! -e "${FILE:?}" ]; then
    "${WGET:?}" "${URL:?}" -q --show-progress --no-check-certificate \
      --progress=bar:force:noscroll --timeout=30 --tries=2 --waitretry=0 \
      --retry-connrefused --read-timeout=30 || { error "${EC24:?}"; return 1; }
  fi

  TMP=$(tar -tf "${FILE:?}" | head -1 | cut -f1 -d"/")
  
  tar -zxf "${FILE:?}" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC25:?}"; return 1; }
  
  mv "${TMP:?}" "${PACKDIR:?}/lensing"

fi

cdfolder "${ROOTDIR}" || return 1;

pbottom "SETUP/UNXV ACT-DR6 DATA (LENSING)" || { unset_all; return 1; }

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

ptop "SETUP/UNXV ACT-DR6 DATA (CMBONLY)" || { unset_all; return 1; }

FOLDER="act_dr6_cmbonly"

# PACK = PACKAGE, DIR = DIRECTORY
PACKDIR="${EDATAF:?}/${FOLDER:?}"

FILE="${ACTDR6_CMBONLY_DATA_FILE:-"dr6_data_cmbonly.tar.gz"}"

URL_BASE="https://lambda.gsfc.nasa.gov/data/act/pspipe/sacc_files/"

URL="${ACTDR6_CMBONLY_DATA_URL:-"${URL_BASE:?}"}/${FILE:?}"

# ---------------------------------------------------------------------------
# in case this script is called twice
# ---------------------------------------------------------------------------
if [ -n "${OVERWRITE_EXISTING_ACTDR6_CMB_DATA:-}" ]; then    
  
  rm -rf "${PACKDIR:?}"
  
  if [ -n "${REDOWNLOAD_EXISTING_ACTDR6_CMB_DATA:-}" ]; then
    rm -f "${EDATAF:?}/${FILE:?}"
  fi

fi

if [ ! -d "${PACKDIR:?}" ]; then
  
  mkdir -p "${PACKDIR:?}" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC20:?}"; return 1; }
        
  cdfolder "${EDATAF:?}" || return 1;

  if [ ! -e "${FILE:?}" ]; then
  
    "${WGET:?}" "${URL:?}" -q --show-progress --no-check-certificate \
      --progress=bar:force:noscroll --timeout=30 --tries=2 --waitretry=0 \
      --retry-connrefused --read-timeout=30 || { error "${EC24:?}"; return 1; }
  
  fi

  
  TMP=$(tar -tf "${FILE:?}" | head -1 | cut -f1 -d"/")
  
  tar -zxvf "${FILE:?}" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC25:?}"; return 1; }
  
  mv "${TMP:?}" "${PACKDIR:?}"

fi

cdfolder "${ROOTDIR}" || return 1;

pbottom "SETUP/UNXV ACT-DR6 DATA (CMBONLY)" || { unset_all; return 1; }

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

ptop "SETUP/UNXV ACT-DR6 DATA (MFLIKE)" || { unset_all; return 1; }

FOLDER="act_dr6_mflike"

PACKDIR="${EDATAF:?}/${FOLDER:?}"

FILE="${ACTDR6_MFLIKE_DATA_FILE:-"dr6_data.tar.gz"}"

URL_BASE="https://lambda.gsfc.nasa.gov/data/act/pspipe/sacc_files/"

URL="${ACTDR6_MFLIKE_DATA_URL:-"${URL_BASE:?}"}/${FILE:?}"

# ---------------------------------------------------------------------------
# note: in case script run >1x w/ previous run stoped prematurely b/c error
# ---------------------------------------------------------------------------
if [ -n "${OVERWRITE_EXISTING_ACTDR6_CMB_DATA:-}" ]; then 
  
  rm -rf "${PACKDIR:?}"
  
  if [ -n "${REDOWNLOAD_EXISTING_ACTDR6_CMB_DATA:-}" ]; then
    rm -f "${EDATAF:?}/${FILE:?}"
  fi

fi

if [ ! -d "${PACKDIR:?}" ]; then
  
  mkdir -p "${PACKDIR:?}" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC20:?}"; return 1; }
        
  cdfolder "${EDATAF:?}" || return 1;

  if [ ! -e "${FILE:?}" ]; then
  
    "${WGET:?}" "${URL:?}" -q --show-progress --no-check-certificate \
      --progress=bar:force:noscroll --timeout=30 --tries=2 --waitretry=0 \
      --retry-connrefused --read-timeout=30 || { error "${EC24:?}"; return 1; }
  
  fi

  TMP=$(tar -tf "${FILE:?}" | head -1 | cut -f1 -d"/")
  tar -zxvf "${FILE:?}" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC25:?}"; return 1; }
  
  mv "${TMP:?}" "${PACKDIR:?}"

fi

pbottom "SETUP/UNXV ACT-DR6 DATA (MFLIKE)" || { unset_all; return 1; }

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
