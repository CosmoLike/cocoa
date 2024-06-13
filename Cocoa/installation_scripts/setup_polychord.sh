#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_POLYCHORD_COMPILATION}" ]; then
  echo -e '\033[1;34m''\tINSTALLING POLYCHORD VIA SETUP_POLYCHORD''\033[0m'

  if [ -z "${ROOTDIR}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE ROOTDIR IS NOT DEFINED''\033[0m'
    return 1
  fi
  if [ -z "${GIT}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE GIT IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${POLY_NAME}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE POLY_NAME IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
  fi

  if [ -z "${DEBUG_POLY_OUTPUT}" ]; then
    export OUTPUT_POLY_1="/dev/null"
    export OUTPUT_POLY_2="/dev/null"
    export POLY_MAKE_NUM_THREADS="${MAKE_NUM_THREADS}"
  else
    export OUTPUT_POLY_1="/dev/tty"
    export OUTPUT_POLY_2="/dev/tty"
    export POLY_MAKE_NUM_THREADS=1
  fi

  rm -rf $ROOTDIR/external_modules/code/$POLY_NAME

  cd $ROOTDIR/external_modules/code/

  $GIT clone https://github.com/PolyChord/PolyChordLite.git $POLY_NAME > ${OUTPUT_POLY_1} 2> ${OUTPUT_POLY_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"SETUP POLYCHORD COULD NOT RUN \e[3mGIT CLONE"'\033[0m'
    cd $ROOTDIR
    unset OUTPUT_POLY_1
    unset OUTPUT_POLY_2
    unset POLY_MAKE_NUM_THREADS
    return 1
  else
    echo -e '\033[0;32m'"\t\t SETUP POLYCHORD RUN \e[3mGIT CLONE\e[0m\e\033[0;32m DONE"'\033[0m'
  fi

  cd $ROOTDIR/external_modules/code/$POLY_NAME

  $GIT checkout $POLYCHORD_GIT_COMMIT > ${OUTPUT_POLY_1} 2> ${OUTPUT_POLY_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"SETUP POLYCHORD COULD NOT RUN \e[3mGIT CHECKOUT"'\033[0m'
    cd $ROOTDIR
    unset OUTPUT_POLY_1
    unset OUTPUT_POLY_2
    unset POLY_MAKE_NUM_THREADS
    return 1
  else
    echo -e '\033[0;32m'"\t\t SETUP POLYCHORD RUN \e[3mGIT CHECKOUT\e[0m\e\033[0;32m DONE"'\033[0m'
  fi

  cp $ROOTDIR/../cocoa_installation_libraries/polychord_changes/Makefile.patch .

  cp $ROOTDIR/../cocoa_installation_libraries/polychord_changes/setup.patch .

  patch -u Makefile -i Makefile.patch > ${OUTPUT_POLY_1} 2> ${OUTPUT_POLY_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"SETUP POLYCHORD COULD NOT RUN \e[3mPATCH MAKEFILE"'\033[0m'
    cd $ROOTDIR
    unset OUTPUT_POLY_1
    unset OUTPUT_POLY_2
    unset POLY_MAKE_NUM_THREADS
    return 1
  else
    echo -e '\033[0;32m'"\t\t SETUP POLYCHORD RUN \e[3mPATCH SETUP\e[0m\e\033[0;32m DONE"'\033[0m'
  fi

  patch -u setup.py -i setup.patch > ${OUTPUT_POLY_1} 2> ${OUTPUT_POLY_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"SETUP POLYCHORD COULD NOT RUN \e[3mPATCH"'\033[0m'
    cd $ROOTDIR
    unset OUTPUT_POLY_1
    unset OUTPUT_POLY_2
    unset POLY_MAKE_NUM_THREADS
    return 1
  else
    echo -e '\033[0;32m'"\t\t SETUP POLYCHORD RUN \e[3mPATCH\e[0m\e\033[0;32m DONE"'\033[0m'
  fi

  cd $ROOTDIR
  unset OUTPUT_POLY_1
  unset OUTPUT_POLY_2
  unset POLY_MAKE_NUM_THREADS
  
  echo -e '\033[1;34m''\t\e[4mSETUP POLYCHORD DONE''\033[0m'
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------