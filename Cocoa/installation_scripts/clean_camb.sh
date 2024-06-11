#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_CAMB_COMPILATION}" ]; then
  echo -e '\033[1;34m''\tCLEANING CAMB''\033[0m'

  if [ -z "${ROOTDIR}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE ROOTDIR IS NOT DEFINED''\033[0m'
    return 1
  fi
  if [ -z "${PYTHON3}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE PYTHON3 IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
  fi

  if [ -z "${DEBUG_CAMB_OUTPUT}" ]; then
    export OUTPUT_CAMB_1="/dev/null"
    export OUTPUT_CAMB_2="/dev/null"
  else
    export OUTPUT_CAMB_1="/dev/tty"
    export OUTPUT_CAMB_2="/dev/tty"
  fi

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  cd $ROOTDIR/external_modules/code/CAMB/
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------

  rm -rf ./build/

  rm -rf ./camb/__pycache__/

  rm -f ./camb/camblib.so

  rm -rf ./forutils/Releaselib/

  $PYTHON3 setup.py clean > ${OUTPUT_CAMB_1} 2> ${OUTPUT_CAMB_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"CAMB COULD NOT RUN \e[3mPYTHON3 SETUP.PY CLEAN"'\033[0m'
  else
    echo -e '\033[0;32m'"\t\t CAMB \e[3mRUN PYTHON3 SETUP.PY CLEAN\e[0m\e\033[0;32m DONE"'\033[0m'
  fi

  cd $ROOTDIR
  unset OUTPUT_CAMB_1
  unset OUTPUT_CAMB_2
  echo -e '\033[1;34m''\t\e[4mCLEANING CAMB DONE''\033[0m'
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------