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

if [ -z "${SKIP_DECOMM_CAMSPEC}" ]; then
  echo -e '\033[0;32m'"\t\t DECOMPRESSING CAMSPEC DATA"'\033[0m'

  rm -f  $ROOTDIR/external_modules/data/planck/CamSpec/CamSpec2021.zip
  rm -rf $ROOTDIR/external_modules/data/planck/CamSpec/CamSpec2021

  cd  $ROOTDIR/external_modules/data/planck/CamSpec/
  
  wget "https://github.com/CobayaSampler/planck_native_data/releases/download/v1/CamSpec2021.zip" > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
  if [ $? -ne 0 ]; then 
    echo -e '\033[0;31m'"\t\t DECOMPRESSING CAMSPEC FAILED"'\033[0m'
    cd $ROOTDIR
    return 1
  fi
  
  unzip CamSpec2021.zip > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"\t\t DECOMPRESSING CAMSPEC FAILED"'\033[0m'
    cd $ROOTDIR
    return 1
  fi

  rm -f  $ROOTDIR/external_modules/data/planck/CamSpec/CamSpec2021.zip

  echo -e '\033[0;32m'"\t\t DECOMPRESSING CAMSPEC DATA DONE"'\033[0m'
fi

cd $ROOTDIR/