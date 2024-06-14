#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_HDF5_INSTALLATION}" ]; then
  pfail() {
    echo -e "\033[0;31m ERROR ENV VARIABLE ${1} IS NOT DEFINED \033[0m"
    unset pfail
  }
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'
    return 1
  fi
  if [ -z "${CXX_COMPILER}" ]; then
    pfail 'CXX_COMPILER'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${C_COMPILER}" ]; then
    pfail 'C_COMPILER'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${FORTRAN_COMPILER}" ]; then
    pfail 'FORTRAN_COMPILER'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${CMAKE}" ]; then
    pfail 'CMAKE'
    cd $ROOTDIR
    return 1
  fi
  unset_env_vars () {
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset HDF5_MNT
    unset pfail
    unset unset_env_vars
  }
  fail () {
    export FAILMSG="\033[0;31m WE CANNOT RUN \e[3m"
    export FAILMSG2="\033[0m"
    echo -e "${FAILMSG} ${1} ${FAILMSG2}"
    unset_env_vars
    unset FAILMSG
    unset FAILMSG2
    unset fail
  }
  if [ -z "${DEBUG_HDF5_PACKAGES}" ]; then
    if [ -z "${MAKE_NUM_THREADS}" ]; then
      pfail 'MAKE_NUM_THREADS'
      cd $ROOTDIR
      return 1
    fi
    export OUT1="/dev/null"
    export OUT2="/dev/null"
    export HDF5_MNT="${MAKE_NUM_THREADS}"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
    export HDF5_MNT=1
  fi

  ptop2 'SETUP_HDF5'
  ptop 'INSTALLING HFD5 LIBRARY'

  if [ -z "${COCOA_HDF5_DIR2}" ]; then
    pfail 'COCOA_HDF5_DIR2'
    cd $ROOTDIR
    return 1
  fi
  
  cd $ROOTDIR/../cocoa_installation_libraries/$COCOA_HDF5_DIR2
  if [ $? -ne 0 ]; then
    fail "CD COCOA_HDF5_DIR2"
    return 1
  fi

  rm -f  CMakeCache.txt
  rm -rf CMakeFiles/
  rm -rf cocoa_HDF5_build/
  
  mkdir cocoa_HDF5_build
  if [ $? -ne 0 ]; then
    fail "MKDIR COCOA_HDF5_BUILD"
    return 1
  fi

  cd cocoa_HDF5_build

  $CMAKE -DBUILD_SHARED_LIBS=TRUE \
    -DCMAKE_INSTALL_PREFIX=$ROOTDIR/.local \
    -DCMAKE_C_COMPILER=$C_COMPILER \
    -DCMAKE_CXX_COMPILER=$CXX_COMPILER \
    -DCMAKE_FC_COMPILER=FORTRAN_COMPILER \
    --log-level=ERROR .. > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "CMAKE"
    return 1
  fi

  make -j $HDF5_MNT > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "MAKE"
    return 1
  fi

  make install > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "MAKE INSTALL"
    return 1
  fi

  unset_env_vars
  pbottom 'INSTALLING HDF5 LIBRARY DONE'
  pbottom2 'SETUP_HDF5'
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------