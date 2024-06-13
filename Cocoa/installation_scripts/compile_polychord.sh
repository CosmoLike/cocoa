#!/bin/bash
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
if [ -z "${IGNORE_POLYCHORD_COMPILATION}" ]; then
  echo -e '\033[1;34m''\tCOMPILING POLYCHORD''\033[0m'

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
  if [ -z "${POLY_NAME}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE POLY_NAME IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
  fi
  
  source $ROOTDIR/installation_scripts/clean_polychord.sh

  if [ -z "${DEBUG_POLY_OUTPUT}" ]; then
    export OUTPUT_POLY_1="/dev/null"
    export OUTPUT_POLY_2="/dev/null"
    export POLY_MAKE_NUM_THREADS="${MAKE_NUM_THREADS}"
  else
    export OUTPUT_POLY_1="/dev/tty"
    export OUTPUT_POLY_2="/dev/tty"
    export POLY_MAKE_NUM_THREADS=1
  fi

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  cd $ROOTDIR/external_modules/code/$POLY_NAME/
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------

  make all > ${OUTPUT_POLY_1} 2> ${OUTPUT_POLY_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"POLYCHORD COULD NOT RUN \e[3mMAKE ALL"'\033[0m'
    cd $ROOTDIR
    unset OUTPUT_POLY_1
    unset OUTPUT_POLY_2
    unset POLY_MAKE_NUM_THREADS
    return 1
  else
    echo -e '\033[0;32m'"\t\t POLYCHORD RUN \e[3mMAKE ALL\e[0m\e\033[0;32m DONE"'\033[0m'
  fi

  make -j $POLY_MAKE_NUM_THREADS pypolychord > ${OUTPUT_POLY_1} 2> ${OUTPUT_POLY_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"POLYCHORD COULD NOT RUN \e[3mMAKE PYPOLYCHORD"'\033[0m'
    cd $ROOTDIR
    unset OUTPUT_POLY_1
    unset OUTPUT_POLY_2
    unset POLY_MAKE_NUM_THREADS
    return 1
  else
    echo -e '\033[0;32m'"\t\t POLYCHORD RUN \e[3mMAKE PYPOLYCHORD\e[0m\e\033[0;32m DONE"'\033[0m'
  fi

  CC=$MPI_CC_COMPILER CXX=$MPI_CXX_COMPILER $PYTHON3 setup.py install \
      --prefix $ROOTDIR/.local  > ${OUTPUT_POLY_1} 2> ${OUTPUT_POLY_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"POLYCHORD COULD NOT RUN \e[3mPYTHON3 SETUP.PY INSTALL"'\033[0m'
    cd $ROOTDIR
    unset OUTPUT_POLY_1
    unset OUTPUT_POLY_2
    unset POLY_MAKE_NUM_THREADS
    return 1
  else
    echo -e '\033[0;32m'"\t\t POLYCHORD RUN \e[3mPYTHON3 SETUP.PY INSTALL\e[0m\e\033[0;32m DONE"'\033[0m'
  fi

  cd $ROOTDIR
  unset OUTPUT_POLY_1
  unset OUTPUT_POLY_2
  unset POLY_MAKE_NUM_THREADS
  echo -e '\033[1;34m''\t\e[4mCOMPILING POLYCHORD DONE''\033[0m'
fi
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------