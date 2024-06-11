#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${SKIP_DECOMM_SPT}" ]; then
  if [ -z "${ROOTDIR}" ]; then
      echo 'ERROR ROOTDIR not defined'
      return
  fi
  if [ -z "${GIT}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE GIT IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${DEBUG_UNXV_CLEAN_ALL}" ]; then
    export OUT_UNXV_1="/dev/null"
    export OUT_UNXV_2="/dev/null"
  else
    export OUT_UNXV_1="/dev/tty"
    export OUT_UNXV_2="/dev/tty"
  fi

  cd $ROOTDIR/external_modules/data

  echo -e '\033[0;32m'"\t\t DECOMPRESSING SPT-3G Y1 DATA"'\033[0m'
  
  rm -rf $ROOTDIR/external_modules/data/spt_3g/
  cd $ROOTDIR/external_modules/data
  
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  export SPTURL='https://github.com/SouthPoleTelescope/spt3g_y1_dist.git'
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------

  $GIT clone $SPTURL spt_3g > ${OUT_UNXV_1} 2> ${OUT_UNXV_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"\t\t DECOMPRESSING SPT-3G Y1 DATA FAILED"'\033[0m'
    cd $ROOTDIR
    unset OUT_UNXV_1
    unset OUT_UNXV_2
    unset SPTURL
    return 1
  fi

  cd $ROOTDIR/external_modules/data/spt_3g
  
  $GIT checkout $SPT3G_DATA_GIT_COMMIT > ${OUT_UNXV_1} 2> ${OUT_UNXV_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"\t\t DECOMPRESSING SPT-3G Y1 DATA FAILED"'\033[0m'
    cd $ROOTDIR
    unset OUT_UNXV_1
    unset OUT_UNXV_2
    unset SPTURL
    return 1
  fi

  cd $ROOTDIR/
  unset SPTUR
  unset OUT_UNXV_1
  unset OUT_UNXV_2
  echo -e '\033[0;32m'"\t\t DECOMPRESSING SPT-3G Y1 DATA DONE"'\033[0m'
fi