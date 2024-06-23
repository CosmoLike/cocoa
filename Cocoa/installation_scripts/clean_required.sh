#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${ROOTDIR}" ]; then
  pfail 'ROOTDIR'; return 1
fi

unset_env_funcs () {
  unset -f cdfolder cpfolder error
  unset -f unset_env_funcs
  cdroot || return 1;
}

unset_all () {
  unset_env_funcs
  unset -f unset_all
  cdroot || return 1;
}

error () {
  fail_script_msg "clean_required.sh" "${1}"
  unset_all
  cdroot || return 1;
}

cdfolder() {
  cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
}

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# ----------------------------------------------------------------------------

ptop 'CLEANING COCOA REQUIRED LIBRARIES' || return 1

cdfolder "${ROOTDIR:?}/../cocoa_installation_libraries" || return 1

sh clean_all || { error "SCRIPT clean_all.sh"; return 1; }

unset -v SETUP_PREREQUISITE_DONE
unset_all || return 1;

pbottom 'CLEANING COCOA REQUIRED LIBRARIES' || return 1

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------