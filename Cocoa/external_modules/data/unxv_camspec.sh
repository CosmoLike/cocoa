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
  
  rm -rf $ROOTDIR/external_modules/data/camspec
  rm -f $ROOTDIR/external_modules/data/camspec2020.tgz

  cd $ROOTDIR/external_modules/data
  
  wget "https://people.ast.cam.ac.uk/~stg20/camspec/camspec2020.tgz" > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
  if [ $? -ne 0 ]; then 
    echo -e '\033[0;31m'"\t\t DECOMPRESSING BAO CAMSPEC FAILED"'\033[0m'
    cd $ROOTDIR
    return 1
  fi
  
  tar xzf camspec2020.tgz > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"\t\t DECOMPRESSING BAO CAMSPEC FAILED"'\033[0m'
    cd $ROOTDIR
    return 1
  fi

  rm -f $ROOTDIR/external_modules/data/camspec2020.tgz

  echo -e '\033[0;32m'"\t\t DECOMPRESSING CAMSPEC DATA DONE"'\033[0m'
fi

cd $ROOTDIR/