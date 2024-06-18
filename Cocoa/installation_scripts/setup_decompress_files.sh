#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
pfail() {
  echo -e "\033[0;31m ERROR ENV VARIABLE ${1} IS NOT DEFINED \033[0m"
  unset pfail
}
if [ -z "${ROOTDIR}" ]; then
  pfail 'ROOTDIR'
  return 1
fi
unset_env_vars_sdf () {
  cd $ROOTDIR
  unset pfail
  unset unset_env_vars_spoly
}
fail_sdf () {
  export MSG="\033[0;31m (setup_decompress_files.sh) WE CANNOT RUN \e[3m"
  export MSG2="\033[0m"
  echo -e "${MSG} ${1} ${MSG2}"  
  unset_env_vars_sdf
  unset MSG
  unset MSG2
  unset fail_sdf
}

ptop2 'SETUP_DECOMPRESS_FILES'

# ----------------------------------------------------------------------------
# ----------------- COCOA_INSTALLATION_LIBRARIES -----------------------------
# ----------------------------------------------------------------------------
if [ -z "${NO_UNXZ_COCOA_INSTALLATION_LIBRARIES}" ]; then
  ptop 'DECOMPRESSING FILES ON /COCOA_INSTALLATION_LIBRARIES'

  cd $ROOTDIR/../cocoa_installation_libraries/
  if [ $? -ne 0 ]; then
    fail_sdf "CD COCOA_INSTALLATION_LIBRARIES"
    return 1
  fi

  sh unxv_all.sh
  if [ $? -ne 0 ]; then
    fail_sdf "SCRIPT unxv_all.sh"
    return 1
  fi

  cd $ROOTDIR
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

  cd $ROOTDIR/external_modules/data/
  if [ $? -ne 0 ]; then
    fail_sdf "CD EXTERNAL_MODULES/DATA"
    return 1
  fi

  sh unxv_all.sh
  if [ $? -ne 0 ]; then
    fail_sdf "SCRIPT unxv_all.sh"
    return 1
  fi

  cd $ROOTDIR
  pbottom 'DECOMPRESSING FILES ON /EXTERNAL_MODULES/DATA'
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
unset_env_vars_sdf
pbottom2 'SETUP_DECOMPRESS_FILES DONE'
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------