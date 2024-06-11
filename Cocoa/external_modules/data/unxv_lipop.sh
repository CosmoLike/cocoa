#!/bin/bash

if [ -z "${ROOTDIR}" ]; then
    echo 'ERROR ROOTDIR not defined'
    return
fi

if [ -z "${LIPOP_DATA_VERSION}" ]; then
    echo 'ERROR LIPOP_DATA_VERSION not defined'
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

if [ -z "${SKIP_DECOMM_LIPOP}" ]; then
  echo -e '\033[0;32m'"\t\t DECOMPRESSING HIL/LOL-LIPOP DATA"'\033[0m'

  rm -f  $ROOTDIR/external_modules/data/planck/planck_2020_hillipop_EE_v$LIPOP_DATA_VERSION.tar.gz
  rm -f  $ROOTDIR/external_modules/data/planck/planck_2020_hillipop_TE_v$LIPOP_DATA_VERSION.tar.gz
  rm -f  $ROOTDIR/external_modules/data/planck/planck_2020_hillipop_TT_v$LIPOP_DATA_VERSION.tar.gz
  rm -f  $ROOTDIR/external_modules/data/planck/planck_2020_hillipop_TTTEEE_v$LIPOP_DATA_VERSION.tar.gz
  rm -f  $ROOTDIR/external_modules/data/planck/planck_2020_lollipop.tar.gz
  rm -rf $ROOTDIR/external_modules/data/planck/planck_2020
  rm -rf $ROOTDIR/external_modules/data/planck/hillipop
  rm -rf $ROOTDIR/external_modules/data/planck/lollipop

  cd  $ROOTDIR/external_modules/data/planck/
  
  export NERSC_WGET="https://portal.nersc.gov/cfs/cmb/planck2020/likelihoods/planck_2020"

  wget "${NERSC_WGET}_hillipop_EE_v${LIPOP_DATA_VERSION}.tar.gz" > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
  if [ $? -ne 0 ]; then 
    echo -e '\033[0;31m'"\t\t DECOMPRESSING EE HILLIPOP FAILED"'\033[0m'
    cd $ROOTDIR
    unset NERSC_WGET
    return 1
  fi

  tar -xvzf planck_2020_hillipop_EE_v$LIPOP_DATA_VERSION.tar.gz > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"\t\t DECOMPRESSING EE HILLIPOP FAILED"'\033[0m'
    cd $ROOTDIR
    unset NERSC_WGET
    return 1
  fi
  
  wget "${NERSC_WGET}_hillipop_TE_v${LIPOP_DATA_VERSION}.tar.gz" > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
  if [ $? -ne 0 ]; then 
    echo -e '\033[0;31m'"\t\t DECOMPRESSING TE HILLIPOP FAILED"'\033[0m'
    cd $ROOTDIR
    unset NERSC_WGET
    return 1
  fi

  tar -xvzf planck_2020_hillipop_TE_v$LIPOP_DATA_VERSION.tar.gz > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"\t\t DECOMPRESSING TE HILLIPOP FAILED"'\033[0m'
    cd $ROOTDIR
    unset NERSC_WGET
    return 1
  fi

  wget "${NERSC_WGET}_hillipop_TTTEEE_v${LIPOP_DATA_VERSION}.tar.gz" > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
  if [ $? -ne 0 ]; then 
    echo -e '\033[0;31m'"\t\t DECOMPRESSING TTTEEE HILLIPOP FAILED"'\033[0m'
    cd $ROOTDIR
    unset NERSC_WGET
    return 1
  fi

  tar -xvzf planck_2020_hillipop_TTTEEE_v$LIPOP_DATA_VERSION.tar.gz > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"\t\t DECOMPRESSING TTTEEE HILLIPOP FAILED"'\033[0m'
    cd $ROOTDIR
    unset NERSC_WGET
    return 1
  fi

  wget "${NERSC_WGET}_hillipop_TT_v${LIPOP_DATA_VERSION}.tar.gz" > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
  if [ $? -ne 0 ]; then 
    echo -e '\033[0;31m'"\t\t DECOMPRESSING TT HILLIPOP FAILED"'\033[0m'
    cd $ROOTDIR
    unset NERSC_WGET
    return 1
  fi

  tar -xvzf planck_2020_hillipop_TT_v$LIPOP_DATA_VERSION.tar.gz > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"\t\t DECOMPRESSING TT HILLIPOP FAILED"'\033[0m'
    cd $ROOTDIR
    return 1
  fi

  wget "${NERSC_WGET}_lollipop.tar.gz" > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
  if [ $? -ne 0 ]; then 
    echo -e '\033[0;31m'"\t\t DECOMPRESSING LOLLIPOP FAILED"'\033[0m'
    cd $ROOTDIR
    unset NERSC_WGET
    return 1
  fi

  tar -xvzf planck_2020_lollipop.tar.gz > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"\t\t DECOMPRESSING LOLLIPOP FAILED"'\033[0m'
    cd $ROOTDIR
    unset NERSC_WGET
    return 1
  fi

  mv $ROOTDIR/external_modules/data/planck/planck_2020/hillipop $ROOTDIR/external_modules/data/planck
  mv $ROOTDIR/external_modules/data/planck/planck_2020/lollipop $ROOTDIR/external_modules/data/planck

  rm -f  $ROOTDIR/external_modules/data/planck/planck_2020_hillipop_EE_v$LIPOP_DATA_VERSION.tar.gz
  rm -f  $ROOTDIR/external_modules/data/planck/planck_2020_hillipop_TE_v$LIPOP_DATA_VERSION.tar.gz
  rm -f  $ROOTDIR/external_modules/data/planck/planck_2020_hillipop_TT_v$LIPOP_DATA_VERSION.tar.gz
  rm -f  $ROOTDIR/external_modules/data/planck/planck_2020_lollipop.tar.gz
  rm -f  $ROOTDIR/external_modules/data/planck/planck_2020_hillipop_TTTEEE_v$LIPOP_DATA_VERSION.tar.gz
  rm -rf $ROOTDIR/external_modules/data/planck/planck_2020

  unset NERSC_WGET

  echo -e '\033[0;32m'"\t\t DECOMPRESSING LIPOP DATA DONE"'\033[0m'
fi

cd $ROOTDIR/