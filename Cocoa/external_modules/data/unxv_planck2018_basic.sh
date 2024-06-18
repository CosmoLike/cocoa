#!/bin/bash

if [ -z "${ROOTDIR}" ]; then
    echo 'ERROR ROOTDIR not defined'
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

if [ -z "${SKIP_DECOMM_PLANCK}" ]; then
  echo -e '\033[0;32m'"\t\t DECOMPRESSING SUPPLEMENTAL DATA AND COVARIANCES"'\033[0m'

  rm -rf $ROOTDIR/external_modules/data/planck/planck_supp_data_and_covmats

  cd $ROOTDIR/external_modules/data/planck
  if [ $? -ne 0 ]; then
    cd $ROOTDIR
    return 1
  fi

  tar xf planck_supp_data_and_covmats.xz
  if [ $? -ne 0 ]; then
    cd $ROOTDIR
    return 1
  fi

  echo -e '\033[0;32m'"\t\t DECOMPRESSING SUPPLEMENTAL DATA AND COVARIANCES DONE"'\033[0m'

  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  
  echo -e '\033[0;32m'"\t\t DECOMPRESSING PLANCK-2018 (PLC-3.0) DATA"'\033[0m'

  rm -rf $ROOTDIR/external_modules/data/planck/plc_3.0/low_l
  rm -rf $ROOTDIR/external_modules/data/planck/plc_3.0/lensing
  rm -rf $ROOTDIR/external_modules/data/planck/plc_3.0/hi_l/plik
  rm -rf $ROOTDIR/external_modules/data/planck/plc_3.0/hi_l/plik_lite

  cd $ROOTDIR/external_modules/data/planck/plc_3.0

  tar xf lensing.xz    > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    cd $ROOTDIR
    return 1
  fi
  tar xf low_l.xz      > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    cd $ROOTDIR
    return 1
  fi

  cd $ROOTDIR/external_modules/data/planck/plc_3.0/hi_l
  rm -rf ./plik/
  tar xf plik.xz       > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    cd $ROOTDIR
    return 1
  fi

  rm -rf ./plik_lite/
  
  tar xf plik_lite.xz  > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    cd $ROOTDIR
    return 1
  fi

  echo -e '\033[0;32m'"\t\t DECOMPRESSING PLANCK-2018 (PLC-3.0) DATA DONE"'\033[0m'
fi

cd $ROOTDIR/