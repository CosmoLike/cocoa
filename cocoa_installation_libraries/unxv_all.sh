#!/bin/bash

if [ -z "${ROOTDIR}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE ROOTDIR IS NOT DEFINED''\033[0m'
    return 1
fi

if [ -z "${DEBUG_UNXV_CLEAN_ALL}" ]; then
  export OUTPUT_UNXV_ALL_1="/dev/null"
  export OUTPUT_UNXV_ALL_2="/dev/null"
else
  export OUTPUT_UNXV_ALL_1="/dev/tty"
  export OUTPUT_UNXV_ALL_2="/dev/tty"
fi

sh $ROOTDIR/../cocoa_installation_libraries/clean_all
cd $ROOTDIR/../cocoa_installation_libraries/

if [ -z "${IGNORE_ALL_PIP_INSTALLATION}" ]; then
  echo -e '\033[0;32m'"\t\tDECOMPRESSING PIP CACHE"'\033[0m' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
  tar xf pip_cache.xz
  if [ -z "${MINICONDA_INSTALLATION}" ]; then
    tar xf expat.xz
  fi
fi

if [ -z "${IGNORE_CPP_INSTALLATION}" ]; then
  if [ -z "${IGNORE_CPP_SPDLOG_INSTALLATION}" ]; then
    echo -e '\033[0;32m'"\t\t DECOMPRESSING CPP SPDLOG LIBRARY"'\033[0m' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
    tar xf spdlog.xz
  fi

  if [ -z "${IGNORE_CPP_CARMA_INSTALLATION}" ]; then
    echo -e '\033[0;32m'"\t\t DECOMPRESSING CPP CARMA LIBRARY"'\033[0m' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
    tar xf carma.xz
  fi
fi

if [ -z "${THREAD_UNXZ}" ]; then
  echo -e '\033[0;32m'"\t\t DECOMPRESSING EE2 LIBRARY"'\033[0m' > ${OUTPUT_UNXV_ALL_1} 2> ${OUTPUT_UNXV_ALL_2}
  tar xf ee2.xz
fi

cd $ROOTDIR/