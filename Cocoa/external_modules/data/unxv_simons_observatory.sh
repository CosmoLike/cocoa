#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${SKIP_DECOMM_SIMONS_OBSERVATORY}" ]; then
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

  echo -e '\033[0;32m'"\t\t DECOMPRESSING SIMONS OBSERVATORY"'\033[0m'

  rm -rf $ROOTDIR/external_modules/data/simons_observatory/

  mkdir -p simons_observatory
  cd $ROOTDIR/external_modules/data/simons_observatory
  
  export SOURL='https://portal.nersc.gov/cfs/sobs/users/MFLike_data/v0.7.1.tar.gz'
  wget $SOURL > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset SOURL
    return 1
  fi

  tar -zxvf v0.7.1.tar.gz > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    cd $ROOTDIR
    unset OUT1
    unset OUT2  
    unset SOURL
    return 1
  fi

  export SOURL='https://portal.nersc.gov/cfs/sobs/users/MFLike_data/v0.8.tar.gz'
  wget $SOURL > ${OUT1} 2> ${OUT2}
   if [ $? -ne 0 ]; then
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset SOURL
    return 1
  fi

  tar -zxvf v0.8.tar.gz > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset SOURL
    return 1
  fi

  cd $ROOTDIR/
  unset OUT1
  unset OUT2
  unset SOURL
  echo -e '\033[0;32m'"\t\t DECOMPRESSING SIMONS OBSERVATORY DONE"'\033[0m'
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------