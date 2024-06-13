#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_HDF5_INSTALLATION}" ]; then
  echo -e '\033[1;44m''SETUP_HDF5''\033[0m'

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
  if [ -z "${CMAKE}" ]; then
      echo -e '\033[0;31m''ERROR ENV VARIABLE MAKE IS NOT DEFINED''\033[0m'
      return 1
  fi

  echo -e '\033[1;34m''\tINSTALLING HFD5 LIBRARY - \e[4mIT WILL TAKE A WHILE''\033[0m'

  unset_env_vars () {
    cd $ROOTDIR
    unset OUTPUT_HDF5_1
    unset OUTPUT_HDF5_2
    unset HDF5_MAKE_NUM_THREADS
  }

  if [ -z "${DEBUG_HDF5_PACKAGES}" ]; then
    export OUTPUT_HDF5_1="/dev/null"
    export OUTPUT_HDF5_2="/dev/null"
    export HDF5_MAKE_NUM_THREADS="${MAKE_NUM_THREADS}"
  else
    export OUTPUT_HDF5_1="/dev/tty"
    export OUTPUT_HDF5_2="/dev/tty"
    export HDF5_MAKE_NUM_THREADS=1
  fi

  cd $ROOTDIR/../cocoa_installation_libraries/$COCOA_HDF5_DIR2
  rm -f  CMakeCache.txt
  rm -rf CMakeFiles/
  rm -rf cocoa_HDF5_build/
  mkdir cocoa_HDF5_build
  cd cocoa_HDF5_build

  $CMAKE -DBUILD_SHARED_LIBS=TRUE -DCMAKE_INSTALL_PREFIX=$ROOTDIR/.local \
    -DCMAKE_C_COMPILER=$C_COMPILER -DCMAKE_CXX_COMPILER=$CXX_COMPILER \
    -DCMAKE_FC_COMPILER=FORTRAN_COMPILER \
    --log-level=ERROR .. > ${OUTPUT_HDF5_1} 2> ${OUTPUT_HDF5_2}
  if [ $? -eq 0 ]; then
    echo -e '\033[0;32m'"\t\t HDF5 RUN \e[3mCMAKE\e[0m\e\033[0;32m DONE"'\033[0m'
  else
    echo -e '\033[0;31m'"HDF5 COULD NOT RUN \e[3mCMAKE"'\033[0m'
    unset_env_vars
    unset unset_env_vars
    return 1
  fi

  make -j $HDF5_MAKE_NUM_THREADS > ${OUTPUT_HDF5_1} 2> ${OUTPUT_HDF5_2}
  if [ $? -eq 0 ]; then
    echo -e '\033[0;32m'"\t\t HDF5 RUN \e[3mMAKE\e[0m\e\033[0;32m DONE"'\033[0m'
  else
    echo -e '\033[0;31m'"HDF5 COULD NOT RUN \e[3mMAKE"'\033[0m'
    unset_env_vars
    unset unset_env_vars
    return 1
  fi

  make install > ${OUTPUT_HDF5_1} 2> ${OUTPUT_HDF5_2}
  if [ $? -eq 0 ]; then
    echo -e '\033[0;32m'"\t\t HDF5 RUN \e[3mMAKE INSTALL\e[0m\e\033[0;32m DONE"'\033[0m'
  else
    echo -e '\033[0;31m'"HDF5 COULD NOT RUN \e[3mMAKE INSTALL"'\033[0m'
    unset_env_vars
    unset unset_env_vars
    return 1
  fi

  unset_env_vars
  unset unset_env_vars
  echo -e '\033[1;34m''\t\e[4mINSTALLING HDF5 LIBRARY DONE''\033[0m'
  echo -e '\033[1;44m''\e[4mSETUP_HDF5 DONE''\033[0m'
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------