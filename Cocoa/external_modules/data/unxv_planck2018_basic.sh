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


if [ -z "${SKIP_DECOMM_PLANCK}" ]; then
  echo -e '\033[0;32m'"\t\t DECOMPRESSING SUPPLEMENTAL DATA AND COVARIANCES"'\033[0m'

  rm -rf $ROOTDIR/external_modules/data/planck/planck_supp_data_and_covmats

  cd $ROOTDIR/external_modules/data/planck
  tar xf planck_supp_data_and_covmats.xz

  echo -e '\033[0;32m'"\t\t DECOMPRESSING SUPPLEMENTAL DATA AND COVARIANCES DONE"'\033[0m'


  echo -e '\033[0;32m'"\t\t DECOMPRESSING PLANCK-2018 (PLC-3.0) DATA"'\033[0m'

  rm -rf $ROOTDIR/external_modules/data/planck/plc_3.0/low_l
  rm -rf $ROOTDIR/external_modules/data/planck/plc_3.0/lensing
  rm -rf $ROOTDIR/external_modules/data/planck/plc_3.0/hi_l/plik
  rm -rf $ROOTDIR/external_modules/data/planck/plc_3.0/hi_l/plik_lite

  cd $ROOTDIR/external_modules/data/planck/plc_3.0
  tar xf lensing.xz    > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
  tar xf low_l.xz      > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
  cd $ROOTDIR/external_modules/data/planck/plc_3.0/hi_l
  rm -rf ./plik/
  tar xf plik.xz       > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
  rm -rf ./plik_lite/
  tar xf plik_lite.xz  > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}

  echo -e '\033[0;32m'"\t\t DECOMPRESSING PLANCK-2018 (PLC-3.0) DATA DONE"'\033[0m'
fi

cd $ROOTDIR/