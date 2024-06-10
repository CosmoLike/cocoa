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

if [ -z "${SKIP_DECOMM_SIMONS_OBSERVATORY}" ]; then
    echo -e '\033[0;32m'"\t\t DECOMPRESSING SIMONS OBSERVATORY"'\033[0m'

    rm -rf $ROOTDIR/external_modules/data/simons_observatory/

    mkdir -p simons_observatory
    cd $ROOTDIR/external_modules/data/simons_observatory
    
    wget https://portal.nersc.gov/cfs/sobs/users/MFLike_data/v0.7.1.tar.gz > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
    tar -zxvf v0.7.1.tar.gz > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
    
    wget https://portal.nersc.gov/cfs/sobs/users/MFLike_data/v0.8.tar.gz > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
    tar -zxvf v0.8.tar.gz > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}

    echo -e '\033[0;32m'"\t\t DECOMPRESSING SIMONS OBSERVATORY DONE"'\033[0m'
fi

cd $ROOTDIR/