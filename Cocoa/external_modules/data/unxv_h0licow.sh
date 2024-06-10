#!/bin/bash

if [ -z "${ROOTDIR}" ]; then
    echo 'ERROR ROOTDIR not defined'
    return
fi

if [ -z "${DEBUG_UNXV_CLEAN_ALL}" ]; then
  export OUTPUT_UNXV_ALL_1="/dev/null"
  export OUTPUT_UNXV_ALL_2="/dev/null"
else
  export OUTPUT_UNXV_ALL_1="/dev/tty"
  export OUTPUT_UNXV_ALL_2="/dev/tty"
fi

cd $ROOTDIR/external_modules/data

if [ -z "${SKIP_DECOMM_STRONG_LENSING}" ]; then
  echo -e '\033[0;32m'"\t\t DECOMPRESSING H0LICOW DATA"'\033[0m'

  rm -rf $ROOTDIR/external_modules/data/holicow_tmp

  git clone https://github.com/shsuyu/H0LiCOW-public.git holicow_tmp > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"\t\t DECOMPRESSING H0LICOW DATA FAILED"'\033[0m'
    cd $ROOTDIR
    return 1
  fi

  cd $ROOTDIR/external_modules/data/holicow_tmp

  git checkout f792647d1fd6c09d9e052fef526669cbd702ab82 > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"\t\t DECOMPRESSING H0LICOW DATA FAILED"'\033[0m'
    cd $ROOTDIR
    return 1
  fi

  mv h0licow_distance_chains ../
  
  cd $ROOTDIR/external_modules/data
  rm -rf $ROOTDIR/external_modules/data/holicow_tmp

  echo -e '\033[0;32m'"\t\t DECOMPRESSING H0LICOW DATA DONE"'\033[0m'
fi

cd $ROOTDIR/