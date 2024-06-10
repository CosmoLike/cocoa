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

if [ -z "${SKIP_DECOMM_SPT}" ]; then
  echo -e '\033[0;32m'"\t\t DECOMPRESSING SPT-3G Y1 DATA"'\033[0m'
  
  rm -rf $ROOTDIR/external_modules/data/spt_3g/
  cd $ROOTDIR/external_modules/data
  
  git clone https://github.com/SouthPoleTelescope/spt3g_y1_dist.git spt_3g > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"\t\t DECOMPRESSING SPT-3G Y1 DATA FAILED"'\033[0m'
    cd $ROOTDIR
    return 1
  fi

  cd $ROOTDIR/external_modules/data/spt_3g
  
  git checkout 66da8e9e2f325024566fe13245788bf8ede897bc > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"\t\t DECOMPRESSING SPT-3G Y1 DATA FAILED"'\033[0m'
    cd $ROOTDIR
    return 1
  fi

  echo -e '\033[0;32m'"\t\t DECOMPRESSING SPT-3G Y1 DATA DONE"'\033[0m'
fi

cd $ROOTDIR/