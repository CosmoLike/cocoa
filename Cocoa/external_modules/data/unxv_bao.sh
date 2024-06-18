#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${SKIP_DECOMM_BAO}" ]; then
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

  echo -e '\033[0;32m'"\t\t DECOMPRESSING BAO DATA"'\033[0m'
  
  rm -rf $ROOTDIR/external_modules/data/bao_data

  cd $ROOTDIR/external_modules/data

  tar xf bao_data.xz > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"\t\t DECOMPRESSING BAO DATA FAILED"'\033[0m'
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    return 1
  fi

  cd $ROOTDIR
  unset OUT1
  unset OUT2
  echo -e '\033[0;32m'"\t\t DECOMPRESSING BAO DATA DONE"'\033[0m'
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------