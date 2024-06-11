#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
echo -e '\033[1;44m''SETUP_OPENBLAS''\033[0m'

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
if [ -z "${MAKE_NUM_THREADS}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE MAKE_NUM_THREADS IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
fi

# ----------------------------------------------------------------------------
# ----------------------------- OPENBLAS LIBRARY  ----------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_OPENBLAS_INSTALLATION}" ]; then
  echo -e '\033[1;34m''INSTALLING OPENBLAS LIBRARY - \e[4mIT WILL TAKE A WHILE''\033[0m'

  if [ -z "${DEBUG_OPENBLAS_PACKAGE}" ]; then
    export OUTPUT_OPENBLAS_1="/dev/null"
    export OUTPUT_OPENBLAS_2="/dev/null"
    export OPENBLAS_MAKE_NUM_THREADS="${MAKE_NUM_THREADS}"
  else
    export OUTPUT_OPENBLAS_1="/dev/tty"
    export OUTPUT_OPENBLAS_2="/dev/tty"
    export OPENBLAS_MAKE_NUM_THREADS=1
  fi

  cd $ROOTDIR/../cocoa_installation_libraries/$COCOA_OPENBLAS_DIR
  
  export MAKE_NB_JOBS=$OPENBLAS_MAKE_NUM_THREADS
  
  make clean > ${OUTPUT_OPENBLAS_1} 2> ${OUTPUT_OPENBLAS_2}

  make CC=$C_COMPILER FC=$FORTRAN_COMPILER \
    USE_OPENMP=1 > ${OUTPUT_OPENBLAS_1} 2> ${OUTPUT_OPENBLAS_2}
  if [ $? -eq 0 ]; then
    echo -e '\033[0;32m'"OPENBLAS RUN \e[3mMAKE\e[0m\e\033[0;32m DONE"'\033[0m'
  else
    echo -e '\033[0;31m'"OPENBLAS COULD NOT RUN \e[3mMAKE"'\033[0m'
    unset MAKE_NB_JOBS
    cd $ROOTDIR
    return 1
  fi

  make install PREFIX=$ROOTDIR/.local > ${OUTPUT_OPENBLAS_1} 2> ${OUTPUT_OPENBLAS_2}
  if [ $? -eq 0 ]; then
    echo -e '\033[0;32m'"OPENBLAS RUN \e[3mMAKE INSTALL\e[0m\e\033[0;32m DONE"'\033[0m'
  else
    echo -e '\033[0;31m'"OPENBLAS COULD NOT RUN \e[3mMAKE INSTALL"'\033[0m'
    unset MAKE_NB_JOBS
    cd $ROOTDIR
    return 1
  fi

  unset MAKE_NB_JOBS

  cd $ROOTDIR/
  echo -e '\033[1;34m''\e[4mINSTALLING OPENBLAS LIBRARY DONE''\033[0m'
fi

echo -e '\033[1;44m''\e[4mSETUP_OPENBLAS DONE''\033[0m'
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------