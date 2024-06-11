#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_PLANCK_COMPILATION}" ]; then
  echo -e '\033[1;34m''\tCLEANING PLANCK LIKELIHOOD''\033[0m'

  if [ -z "${ROOTDIR}" ]; then
      echo -e '\033[0;31m''ERROR ENV VARIABLE ROOTDIR IS NOT DEFINED''\033[0m'
      return 1
  fi
  if [ -z "${PYTHON3}" ]; then
      echo -e '\033[0;31m''ERROR ENV VARIABLE PYTHON3 IS NOT DEFINED''\033[0m'
      cd $ROOTDIR
      return 1
  fi
  if [ -z "${DEBUG_PLANCK_OUTPUT}" ]; then
    export OUTPUT_PLANCK_1="/dev/null"
    export OUTPUT_PLANCK_2="/dev/null"
  else
    export OUTPUT_PLANCK_1="/dev/tty"
    export OUTPUT_PLANCK_2="/dev/tty"
  fi

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  if [ -z "${USE_SPT_CLIK_PLANCK}" ]; then
    cd $ROOTDIR/external_modules/code/planck/code/plc_3.0/plc-3.1/
  else
    cd $ROOTDIR/external_modules/code/planck/code/spt_clik/
  fi
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------

  rm -f $ROOTDIR/.local/bin/clik*

  rm -f $ROOTDIR/.local/lib/libclik_f90.so

  rm -f $ROOTDIR/.local/lib/libclik.so

  rm -rf $ROOTDIR/.local/lib/python/site-packages/clik

  rm -rf $ROOTDIR/.local/share/clik

  rm -f $ROOTDIR/.local/include/clik*

  rm -f .lock-waf_*

  $PYTHON3 waf distclean > ${OUTPUT_PLANCK_1} 2> ${OUTPUT_PLANCK_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"PLANCK COULD NOT RUN \e[3mWAF DISTCLEAN"'\033[0m'
  else
    echo -e '\033[0;32m'"\t\t PLANCK RUN \e[3mWAF CLEAN\e[0m\e\033[0;32m DONE"'\033[0m'
  fi

  cd $ROOTDIR
  unset OUTPUT_PLANCK_1
  unset OUTPUT_PLANCK_2
  echo -e '\033[1;34m''\t\e[4mCLEANING PLANCK LIKELIHOOD DONE''\033[0m'
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------