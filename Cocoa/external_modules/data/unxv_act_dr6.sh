#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${SKIP_DECOMM_ACT}" ]; then
  if [ -z "${ROOTDIR}" ]; then
      echo 'ERROR ROOTDIR not defined'
      return
  fi
  if [ -z "${DEBUG_UNXV_CLEAN_ALL}" ]; then
    export OUT_UNXV_1="/dev/null"
    export OUT_UNXV_2="/dev/null"
  else
    export OUT_UNXV_1="/dev/tty"
    export OUT_UNXV_2="/dev/tty"
  fi

  echo -e '\033[0;32m'"\t\t DECOMPRESSING ACT-DR6 DATA"'\033[0m'

  cd $ROOTDIR/external_modules/data
  rm -rf $ROOTDIR/external_modules/data/act/
  mkdir -p act
  
  cd $ROOTDIR/external_modules/data/act
  mkdir -p lensing
  
  cd $ROOTDIR/external_modules/data/act/lensing

  export ACTURL="https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr6/likelihood/data/ACT_dr6_likelihood_v1.2.tgz"
  
  wget $ACTURL > ${OUT_UNXV_1} 2> ${OUT_UNXV_2}
  if [ $? -ne 0 ]; then 
    echo -e '\033[0;31m'"\t\t DECOMPRESSING ACT-DR6 DATA FAILED"'\033[0m'
    cd $ROOTDIR
    unset OUT_UNXV_1
    unset OUT_UNXV_2
    return 1
  fi

  tar -zxvf ACT_dr6_likelihood_v1.2.tgz > ${OUT_UNXV_1} 2> ${OUT_UNXV_2}
  if [ $? -ne 0 ]; then 
    echo -e '\033[0;31m'"\t\t DECOMPRESSING ACT-DR6 DATA FAILED"'\033[0m'
    cd $ROOTDIR
    unset OUT_UNXV_1
    unset OUT_UNXV_2
    return 1
  fi
  rm -f $ROOTDIR/external_modules/data/act/lensing/ACT_dr6_likelihood_v1.2.tgz

  cd $ROOTDIR/
  unset OUT_UNXV_1
  unset OUT_UNXV_2
  echo -e '\033[0;32m'"\t\t DECOMPRESSING ACT-DR6 DATA DONE"'\033[0m'
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------