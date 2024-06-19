#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
pfail() {
  echo -e "\033[0;31m ERROR ENV VARIABLE ${1} IS NOT DEFINED \033[0m"
  unset pfail
}
if [ -z "${ROOTDIR}" ]; then
  pfail 'ROOTDIR'; return 1
fi
cdroot() {
  cd "${ROOTDIR}" 2>"/dev/null" || { echo -e \
    "\033[0;31m\t\t CD ROOTDIR (${ROOTDIR}) FAILED \033[0m"; return 1; }
  unset cdroot
}
unset_env_vars_sdf () {
  unset pfail
  unset unset_env_vars_sdf
  cdroot || return 1;
}
fail_sdf () {
  local MSG="\033[0;31m\t\t (setup_decompress_files.sh) WE CANNOT RUN \e[3m"
  local MSG2="\033[0m"
  echo -e "${MSG} ${1} ${MSG2}"  
  unset fail_sdf
  unset_env_vars_sdf
}
cdfolder() {
  cd "${1}" 2>"/dev/null" || { fail_sdf "CD FOLDER: ${1}"; return 1; }
}
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
ptop2 'SETUP_DECOMPRESS_FILES'

# ----------------------------------------------------------------------------
# ----------------- COCOA_INSTALLATION_LIBRARIES -----------------------------
# ----------------------------------------------------------------------------
if [ -z "${NO_UNXZ_COCOA_INSTALLATION_LIBRARIES}" ]; then
  ptop 'DECOMPRESSING FILES ON /COCOA_INSTALLATION_LIBRARIES'

  cdfolder "${ROOTDIR}/../cocoa_installation_libraries/" || return 1;

  sh unxv_all.sh || { fail_sdf "SCRIPT unxv_all.sh"; return 1; }

  cdfolder "${ROOTDIR}" || return 1;

  pbottom 'DECOMPRESSING FILES ON /COCOA_INSTALLATION_LIBRARIES'
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
  ptop 'DECOMPRESSING FILES ON /EXTERNAL_MODULES/DATA'

  cdfolder "${ROOTDIR}/external_modules/data/" || return 1;

  sh unxv_all.sh || { fail_sdf "SCRIPT unxv_all.sh"; return 1; }

  cdfolder "${ROOTDIR}" || return 1;

  pbottom 'DECOMPRESSING FILES ON /EXTERNAL_MODULES/DATA'
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

unset_env_vars_sdf
pbottom2 'SETUP_DECOMPRESS_FILES DONE'
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------