#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${ROOTDIR}" ]; then
  pfail 'ROOTDIR'; return 1
fi

# parenthesis = run in a subshell
( source "${ROOTDIR:?}/installation_scripts/.check_flags.sh" ) || return 1;

if [ -z "${ROOTDIR}" ]; then
  pfail 'ROOTDIR'; return 1
fi

unset_env_vars () {
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
  fail_script_msg "$(basename ${BASH_SOURCE[0]})" "${1}"
  unset_all || return 1
}

cdfolder() {
  cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
}

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

ptop2 'SETUP_DECOMPRESS_FILES' || return 1

# ----------------------------------------------------------------------------
# ----------------- COCOA_INSTALLATION_LIBRARIES -----------------------------
# ----------------------------------------------------------------------------
if [ -z "${NO_UNXZ_COCOA_INSTALLATION_LIBRARIES}" ]; then

  ptop 'DECOMPRESSING FILES ON /COCOA_INSTALLATION_LIBRARIES' || return 1

  cdfolder "${ROOTDIR:?}/../cocoa_installation_libraries/" || return 1;

  # parenthesis = run in a subshell 
  ( sh unxv_all.sh )  || { error "${EC22:?} unxv_all.sh"; return 1; }

  cdfolder "${ROOTDIR:?}" || return 1;

  pbottom 'DECOMPRESSING FILES ON /COCOA_INSTALLATION_LIBRARIES' || return 1

fi

# ----------------------------------------------------------------------------
# ----------------------- /EXTERNAL_MODULES/CODE -----------------------------
# ----------------------------------------------------------------------------
#if [ -z "${NO_UNXZ_EXTERNAL_MODULES_CODE}" ]; then
#    # no xz files under code so far
#fi

# ----------------------------------------------------------------------------
# ----------------------- /EXTERNAL_MODULES/DATA -----------------------------
# ----------------------------------------------------------------------------
if [ -z "${NO_UNXZ_EXTERNAL_MODULES_DATA}" ]; then

  ptop 'DECOMPRESSING FILES ON /EXTERNAL_MODULES/DATA' || return 1

  cdfolder "${ROOTDIR:?}/external_modules/data/" || return 1;

  # parenthesis = run in a subshell
  ( sh unxv_all.sh ) || { error "SCRIPT unxv_all.sh"; return 1; }

  cdfolder "${ROOTDIR:?}" || return 1;

  pbottom 'DECOMPRESSING FILES ON /EXTERNAL_MODULES/DATA' || return 1

fi

# ----------------------------------------------------------------------------

unset_all || return 1; 

pbottom2 'SETUP_DECOMPRESS_FILES DONE' || return 1

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------