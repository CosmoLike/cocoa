#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_FORTRAN_INSTALLATION}" ]; then
  echo -e '\033[1;44m''SETUP_FORTRAN_PACKAGES''\033[0m'

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
  if [ -z "${CMAKE}" ]; then
      echo -e '\033[0;31m''ERROR ENV VARIABLE MAKE IS NOT DEFINED''\033[0m'
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


  if [ -z "${DEBUG_FORTRAN_PACKAGES}" ]; then
    export OUTPUT_FORTRAN_1="/dev/null"
    export OUTPUT_FORTRAN_2="/dev/null"
    export FORTRAN_MAKE_NUM_THREADS="${MAKE_NUM_THREADS}"
  else
    export OUTPUT_FORTRAN_1="/dev/tty"
    export OUTPUT_FORTRAN_2="/dev/tty"
    export FORTRAN_MAKE_NUM_THREADS=1
  fi

  # ----------------------------------------------------------------------------
  # ---------------------------------- FORTRAN LAPACK --------------------------
  # ----------------------------------------------------------------------------
  if [ -z "${IGNORE_FORTRAN_LAPACK_INSTALLATION}" ]; then
    echo -e '\033[1;34m''INSTALLING LAPACK FORTRAN LIBRARY - \e[4mIT WILL TAKE A WHILE''\033[0m'

    cd $ROOTDIR/../cocoa_installation_libraries/
    rm -rf lapack-build
    mkdir lapack-build
    cd ./lapack-build

    $CMAKE -DBUILD_SHARED_LIBS=TRUE -DCMAKE_INSTALL_PREFIX=$ROOTDIR/.local \
      -DCMAKE_C_COMPILER=$C_COMPILER --log-level=ERROR \
      ../${COCOA_LAPACK_DIR} > ${OUTPUT_FORTRAN_1} 2> ${OUTPUT_FORTRAN_2}
    if [ $? -eq 0 ]; then
      echo -e '\033[0;32m'"\t\t  LAPACK RUN \e[3mCMAKE\e[0m\e\033[0;32m DONE"'\033[0m'
    else
      echo -e '\033[0;31m'"LAPACK COULD NOT RUN \e[3mCMAKE"'\033[0m'
      cd $ROOTDIR
      unset OUTPUT_FORTRAN_1
      unset OUTPUT_FORTRAN_2
      unset FORTRAN_MAKE_NUM_THREADS
      return 1
    fi

    make -j $FORTRAN_MAKE_NUM_THREADS all > ${OUTPUT_FORTRAN_1} 2> ${OUTPUT_FORTRAN_2}
    if [ $? -eq 0 ]; then
      echo -e '\033[0;32m'"\t\t  LAPACK RUN \e[3mMAKE\e[0m\e\033[0;32m DONE"'\033[0m'
    else
      echo -e '\033[0;31m'"LAPACK COULD NOT RUN \e[3mMAKE"'\033[0m'
      cd $ROOTDIR
      unset OUTPUT_FORTRAN_1
      unset OUTPUT_FORTRAN_2
      unset FORTRAN_MAKE_NUM_THREADS
      return 1
    fi

    make install > ${OUTPUT_FORTRAN_1} 2> ${OUTPUT_FORTRAN_2}
    if [ $? -eq 0 ]; then
      echo -e '\033[0;32m'"\t\t LAPACK RUN \e[3mMAKE INSTALL\e[0m\e\033[0;32m DONE"'\033[0m'
    else
      echo -e '\033[0;31m'"LAPACK COULD NOT RUN \e[3mMAKE INSTALL"'\033[0m'
      cd $ROOTDIR
      unset OUTPUT_FORTRAN_1
      unset OUTPUT_FORTRAN_2
      unset FORTRAN_MAKE_NUM_THREADS
      return 1
    fi
    
    cd $ROOTDIR    
    echo -e '\033[1;34m''\t\e[4mINSTALLING LAPACK FORTRAN LIBRARY DONE''\033[0m'
  fi
  
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  cd $ROOTDIR
  unset OUTPUT_FORTRAN_1
  unset OUTPUT_FORTRAN_2
  unset FORTRAN_MAKE_NUM_THREADS

  echo -e '\033[1;44m''\e[4mSETUP_FORTRAN_PACKAGES DONE''\033[0m'
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------