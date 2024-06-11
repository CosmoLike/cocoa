#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${SKIP_DECOMM_LIPOP}" ]; then
  if [ -z "${ROOTDIR}" ]; then
      echo 'ERROR ROOTDIR not defined'
      return
  fi
  if [ -z "${LIPOP_DATA_VERSION}" ]; then
      echo 'ERROR LIPOP_DATA_VERSION not defined'
      return
  fi
  if [ -z "${DEBUG_UNXV_CLEAN_ALL}" ]; then
    export OUT_UNXV_1="/dev/null"
    export OUT_UNXV_2="/dev/null"
  else
    export OUT_UNXV_1="/dev/tty"
    export OUT_UNXV_2="/dev/tty"
  fi

  cd $ROOTDIR/external_modules/data

  echo -e '\033[0;32m'"\t\t DECOMPRESSING HIL/LOL-LIPOP DATA"'\033[0m'

  export PLK_DIR="${ROOTDIR}/external_modules/data/planck"
  export LPVS=$LIPOP_DATA_VERSION

  rm -f  "${PLK_DIR}/planck_2020_hillipop_EE_v$LPVS.tar.gz"
  rm -f  "${PLK_DIR}/planck_2020_hillipop_TE_v$LPVS.tar.gz"
  rm -f  "${PLK_DIR}/planck_2020_hillipop_TT_v$LPVS.tar.gz"
  rm -f  "${PLK_DIR}/planck_2020_hillipop_TTTEEE_v$LPVS.tar.gz"
  rm -f  "${PLK_DIR}/planck_2020_lollipop.tar.gz"
  rm -rf "${PLK_DIR}/planck_2020"
  rm -rf "${PLK_DIR}/hillipop"
  rm -rf "${PLK_DIR}/lollipop"

  cd  $ROOTDIR/external_modules/data/planck/
  
  export PL2020_URL="https://portal.nersc.gov/cfs/cmb/planck2020/likelihoods/planck_2020"

  wget "${PL2020_URL}_hillipop_EE_v${LPVS}.tar.gz" > ${OUT_UNXV_1} 2> ${OUT_UNXV_2}
  if [ $? -ne 0 ]; then 
    echo -e '\033[0;31m'"\t\t DECOMPRESSING EE HILLIPOP FAILED"'\033[0m'
    cd $ROOTDIR
    unset OUT_UNXV_1
    unset OUT_UNXV_2
    unset PL2020_URL
    unset PLK_DIR
    unset LPVS
    return 1
  fi

  tar -xvzf planck_2020_hillipop_EE_v$LPVS.tar.gz > ${OUT_UNXV_1} 2> ${OUT_UNXV_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"\t\t DECOMPRESSING EE HILLIPOP FAILED"'\033[0m'
    cd $ROOTDIR
    unset OUT_UNXV_1
    unset OUT_UNXV_2
    unset PL2020_URL
    unset PLK_DIR
    unset LPVS
    return 1
  fi
  
  wget "${PL2020_URL}_hillipop_TE_v${LPVS}.tar.gz" > ${OUT_UNXV_1} 2> ${OUT_UNXV_2}
  if [ $? -ne 0 ]; then 
    echo -e '\033[0;31m'"\t\t DECOMPRESSING TE HILLIPOP FAILED"'\033[0m'
    cd $ROOTDIR
    unset OUT_UNXV_1
    unset OUT_UNXV_2
    unset PL2020_URL
    unset PLK_DIR
    unset LPVS
    return 1
  fi

  tar -xvzf planck_2020_hillipop_TE_v$LPVS.tar.gz > ${OUT_UNXV_1} 2> ${OUT_UNXV_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"\t\t DECOMPRESSING TE HILLIPOP FAILED"'\033[0m'
    cd $ROOTDIR
    unset OUT_UNXV_1
    unset OUT_UNXV_2
    unset PL2020_URL
    unset PLK_DIR
    unset LPVS
    return 1
  fi

  wget "${PL2020_URL}_hillipop_TTTEEE_v${LPVS}.tar.gz" > ${OUT_UNXV_1} 2> ${OUT_UNXV_2}
  if [ $? -ne 0 ]; then 
    echo -e '\033[0;31m'"\t\t DECOMPRESSING TTTEEE HILLIPOP FAILED"'\033[0m'
    cd $ROOTDIR
    unset OUT_UNXV_1
    unset OUT_UNXV_2
    unset PL2020_URL
    unset PLK_DIR
    unset LPVS
    return 1
  fi

  tar -xvzf planck_2020_hillipop_TTTEEE_v$LPVS.tar.gz > ${OUT_UNXV_1} 2> ${OUT_UNXV_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"\t\t DECOMPRESSING TTTEEE HILLIPOP FAILED"'\033[0m'
    cd $ROOTDIR
    unset OUT_UNXV_1
    unset OUT_UNXV_2
    unset PL2020_URL
    unset PLK_DIR
    unset LPVS
    return 1
  fi

  wget "${PL2020_URL}_hillipop_TT_v${LPVS}.tar.gz" > ${OUT_UNXV_1} 2> ${OUT_UNXV_2}
  if [ $? -ne 0 ]; then 
    echo -e '\033[0;31m'"\t\t DECOMPRESSING TT HILLIPOP FAILED"'\033[0m'
    cd $ROOTDIR
    unset OUT_UNXV_1
    unset OUT_UNXV_2
    unset PL2020_URL
    unset PLK_DIR
    unset LPVS
    return 1
  fi

  tar -xvzf planck_2020_hillipop_TT_v$LPVS.tar.gz > ${OUT_UNXV_1} 2> ${OUT_UNXV_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"\t\t DECOMPRESSING TT HILLIPOP FAILED"'\033[0m'
    cd $ROOTDIR
    unset OUT_UNXV_1
    unset OUT_UNXV_2
    unset PL2020_URL
    unset PLK_DIR
    unset LPVS
    return 1
  fi

  wget "${PL2020_URL}_lollipop.tar.gz" > ${OUT_UNXV_1} 2> ${OUT_UNXV_2}
  if [ $? -ne 0 ]; then 
    echo -e '\033[0;31m'"\t\t DECOMPRESSING LOLLIPOP FAILED"'\033[0m'
    cd $ROOTDIR
    unset PL2020_URL
    unset PLK_DIR
    unset LPVS
    return 1
  fi

  tar -xvzf planck_2020_lollipop.tar.gz > ${OUT_UNXV_1} 2> ${OUT_UNXV_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"\t\t DECOMPRESSING LOLLIPOP FAILED"'\033[0m'
    cd $ROOTDIR
    unset OUT_UNXV_1
    unset OUT_UNXV_2
    unset PL2020_URL
    unset PLK_DIR
    unset LPVS
    return 1
  fi

  mv "${PLK_DIR}/planck_2020/hillipop" $PLK_DIR
  mv "${PLK_DIR}/planck_2020/lollipop" $PLK_DIR

  rm -f  "${PLK_DIR}/planck_2020_hillipop_EE_v${LIPOP_DATA_VERSION}.tar.gz"
  rm -f  "${PLK_DIR}/planck_2020_hillipop_TE_v${LIPOP_DATA_VERSION}.tar.gz"
  rm -f  "${PLK_DIR}/planck_2020_hillipop_TT_v${LIPOP_DATA_VERSION}.tar.gz"
  rm -f  "${PLK_DIR}/planck_2020_lollipop.tar.gz"
  rm -f  "${PLK_DIR}/planck_2020_hillipop_TTTEEE_v${LIPOP_DATA_VERSION}.tar.gz"
  rm -rf "${PLK_DIR}/planck_2020"

  cd $ROOTDIR/
  unset OUT_UNXV_1
  unset OUT_UNXV_2
  unset PL2020_URL
  unset PLK_DIR
  unset LPVS
  echo -e '\033[0;32m'"\t\t DECOMPRESSING LIPOP DATA DONE"'\033[0m'
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------