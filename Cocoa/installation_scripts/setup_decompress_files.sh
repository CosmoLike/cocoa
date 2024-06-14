#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${ROOTDIR}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE ROOTDIR IS NOT DEFINED''\033[0m'
    return 1
fi

ptop2 'SETUP_DECOMPRESS_FILES'


# ----------------------------------------------------------------------------
# ----------------- COCOA_INSTALLATION_LIBRARIES -----------------------------
# ----------------------------------------------------------------------------
if [ -z "${NO_UNXZ_COCOA_INSTALLATION_LIBRARIES}" ]; then
  ptop 'DECOMPRESSING FILES ON /COCOA_INSTALLATION_LIBRARIES'

  cd $ROOTDIR/../cocoa_installation_libraries/

  sh unxv_all.sh
  if [ $? -ne 0 ]; then
    cd $ROOTDIR
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

  sh unxv_all.sh
  if [ $? -ne 0 ]; then
    cd $ROOTDIR
    return 1
  fi

  cd $ROOTDIR

  pbottom 'DECOMPRESSING FILES ON /EXTERNAL_MODULES/DATA'
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
pbottom2 'SETUP_DECOMPRESS_FILES DONE'
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------