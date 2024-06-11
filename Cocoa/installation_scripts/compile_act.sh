#!/bin/bash
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
if [ -z "${IGNORE_ACT_COMPILATION}" ]; then
  echo -e '\033[1;34m''\tCOMPILING ACT''\033[0m'

  if [ -z "${ROOTDIR}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE ROOTDIR IS NOT DEFINED''\033[0m'
    return 1
  fi
  if [ -z "${CXX_COMPILER}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE CXX_COMPILER IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${C_COMPILER}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE C_COMPILER IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${FORTRAN_COMPILER}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE FORTRAN_COMPILER IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${PYTHON3}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE PYTHON3 IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${MAKE_NUM_THREADS}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE MAKE_NUM_THREADS IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
  fi

  source $ROOTDIR/installation_scripts/clean_act.sh

  if [ -z "${DEBUG_ACT_OUTPUT}" ]; then
    export OUTPUT_ACT_1="/dev/null"
    export OUTPUT_ACT_2="/dev/null"
  else
    export OUTPUT_ACT_1="/dev/tty"
    export OUTPUT_ACT_2="/dev/tty"
  fi

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  cd $ROOTDIR/external_modules/code/pyactlike/
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------

  $PIP3 install . --prefix=$ROOTDIR/.local > ${OUTPUT_ACT_1} 2> ${OUTPUT_ACT_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"ACT COULD NOT RUN \e[3mPIP3 INSTALL ."'\033[0m'
    cd $ROOTDIR
    unset OUTPUT_ACT_1
    unset OUTPUT_ACT_2
    return 1
  else
    echo -e '\033[0;32m'"\t\t ACT RUN \e[3mPIP3 INSTALL .\e[0m\e\033[0;32m DONE"'\033[0m'
  fi

  cd $ROOTDIR
  unset OUTPUT_ACT_1
  unset OUTPUT_ACT_2
  echo -e '\033[1;34m''\t\e[4mCOMPILING ACT DONE''\033[0m'
fi
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------