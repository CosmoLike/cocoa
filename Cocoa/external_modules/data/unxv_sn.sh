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

if [ -z "${SKIP_DECOMM_SN}" ]; then
  echo -e '\033[0;32m'"\t\t DECOMPRESSING SN DATA"'\033[0m'
  
  rm -rf $ROOTDIR/external_modules/data/sn_data

  cd $ROOTDIR/external_modules/data

  tar xf sn_data.xz > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"\t\t DECOMPRESSING SN DATA FAILED"'\033[0m'
    cd $ROOTDIR
    return 1
  fi

  echo -e '\033[0;32m'"\t\t DECOMPRESSING SN DATA DONE"'\033[0m'
fi

cd $ROOTDIR/