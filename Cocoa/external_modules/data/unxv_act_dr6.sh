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

if [ -z "${SKIP_DECOMM_ACT}" ]; then
  echo -e '\033[0;32m'"\t\t DECOMPRESSING ACT-DR6 DATA"'\033[0m'

  rm -rf act/
  mkdir -p act
  cd $ROOTDIR/external_modules/data/act
  mkdir -p lensing
  cd $ROOTDIR/external_modules/data/act/lensing
  wget https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr6/likelihood/data/ACT_dr6_likelihood_v1.2.tgz > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
  tar -zxvf ACT_dr6_likelihood_v1.2.tgz > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
  cd $ROOTDIR/external_modules/data

  echo -e '\033[0;32m'"\t\t DECOMPRESSING ACT-DR6 DATA DONE"'\033[0m'
fi

cd $ROOTDIR/