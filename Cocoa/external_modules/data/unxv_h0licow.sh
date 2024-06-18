#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${SKIP_DECOMM_STRONG_LENSING}" ]; then
  if [ -z "${ROOTDIR}" ]; then
      echo 'ERROR ROOTDIR not defined'
      return
  fi
  if [ -z "${GIT}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE GIT IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${HOLICOW_DATA_GIT_COMMIT}" ]; then
      echo 'ERROR HOLICOW_DATA_GIT_COMMIT not defined'
      return
  fi
  if [ -z "${DEBUG_UNXV_CLEAN_ALL}" ]; then
    export OUT1="/dev/null"
    export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
  fi

  cd $ROOTDIR/external_modules/data

  echo -e '\033[0;32m'"\t\t DECOMPRESSING H0LICOW DATA"'\033[0m'

  rm -rf $ROOTDIR/external_modules/data/holicow_tmp

  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  export HOLIURL='https://github.com/shsuyu/H0LiCOW-public.git'
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------

  git clone $HOLIURL holicow_tmp > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"\t\t DECOMPRESSING H0LICOW DATA FAILED"'\033[0m'
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset HOLIURL
    return 1
  fi

  cd $ROOTDIR/external_modules/data/holicow_tmp

  git checkout $HOLICOW_DATA_GIT_COMMIT > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"\t\t DECOMPRESSING H0LICOW DATA FAILED"'\033[0m'
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset HOLIURL
    return 1
  fi

  mv h0licow_distance_chains ../
  cd $ROOTDIR/external_modules/data
  rm -rf $ROOTDIR/external_modules/data/holicow_tmp

  cd $ROOTDIR/
  unset OUT1
  unset OUT2
  unset HOLIURL
  echo -e '\033[0;32m'"\t\t DECOMPRESSING H0LICOW DATA DONE"'\033[0m'
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------