#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_CLASS_COMPILATION}" ]; then
  echo -e '\033[1;34m''\tCLEANING CLASS''\033[0m'

  if [ -z "${ROOTDIR}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE ROOTDIR IS NOT DEFINED''\033[0m'
    return 1
  fi
  if [ -z "${PYTHON3}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE PYTHON3 IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${DEBUG_CLASS_OUTPUT}" ]; then
    export OUTPUT_CLASS_1="/dev/null"
    export OUTPUT_CLASS_2="/dev/null"
  else
    export OUTPUT_CLASS_1="/dev/tty"
    export OUTPUT_CLASS_2="/dev/tty"
  fi

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  cd $ROOTDIR/external_modules/code/class_public/
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------

  rm $ROOTDIR/.local/lib/python$PYTHON_VERSION/site-packages/classy*

  #Workaround Cocoa .gitignore entry on /include
  rm -rf ./include > ${OUTPUT_CLASS_1} 2> ${OUTPUT_CLASS_2}

  make clean > ${OUTPUT_CLASS_1} 2> ${OUTPUT_CLASS_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"CLASS COULD NOT RUN \e[3mMAKE CLEAN"'\033[0m'
  else
    echo -e '\033[0;32m'"\t\t CLASS RUN \e[3mMAKE CLEAN\e[0m\e\033[0;32m DONE"'\033[0m'
  fi

  rm -f class

  cd ./python

  $PYTHON3 setup.py clean > ${OUTPUT_CLASS_1} 2> ${OUTPUT_CLASS_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"CLASS COULD NOT RUN \e[3mPYTHON3 SETUP.PY CLEAN"'\033[0m'
  else
    echo -e '\033[0;32m'"\t\t CLASS RUN \e[3mPYTHON3 SETUP.PY CLEAN\e[0m\e\033[0;32m DONE"'\033[0m'
  fi

  rm -rf ./build/
  rm -rf ./classy.egg-info

  cd $ROOTDIR
  unset OUTPUT_CLASS_1
  unset OUTPUT_CLASS_2
  echo -e '\033[1;34m''\t\eCLEANING CLASS DONE''\033[0m'
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------