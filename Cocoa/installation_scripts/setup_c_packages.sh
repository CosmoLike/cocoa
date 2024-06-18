#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_C_INSTALLATION}" ]; then
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
  unset_env_vars_scp () {
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset C_MNT
    unset pfail
    unset unset_env_vars_scp
  }
  fail_scp () {
    export fail_scpMSG="\033[0;31m (setup_c_packages.sh) WE CANNOT RUN \e[3m"
    export fail_scpMSG2="\033[0m"
    echo -e "${fail_scpMSG} ${1} ${fail_scpMSG2}"
    unset_env_vars_scp
    unset fail_scpMSG
    unset fail_scpMSG2
    unset fail_scp
  }
  if [ -z "${DEBUG_C_PACKAGES}" ]; then
    if [ -z "${MAKE_NUM_THREADS}" ]; then
      pfail 'MAKE_NUM_THREADS'
      cd $ROOTDIR
      return 1
    fi
    export OUT1="/dev/null"
    export OUT2="/dev/null"
    export C_MNT="${MAKE_NUM_THREADS}"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
    export C_MNT=1
  fi

  ptop2 'SETUP_C_PACKAGES'

  # ----------------------------------------------------------------------------
  # -------------------------------- FFTW --------------------------------------
  # ----------------------------------------------------------------------------
  if [ -z "${IGNORE_C_FFTW_INSTALLATION}" ]; then
    ptop 'INSTALLING FFTW C LIBRARY'
    
    if [ -z "${COCOA_FFTW_DIR}" ]; then
      pfail 'COCOA_FFTW_DIR'
      return 1
    fi
    
    cd $ROOTDIR/../cocoa_installation_libraries/$COCOA_FFTW_DIR
    if [ $? -ne 0 ]; then
      fail_scp "CD COCOA_FFTW_DIR"
      return 1
    fi

    FC=$FORTRAN_COMPILER CC=$C_COMPILER ./configure \
      --enable-openmp \
      --prefix=$ROOTDIR/.local \
      --enable-shared=yes \
      --enable-static=yes > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_scp "CONFIGURE"
      return 1
    fi

    make -j $C_MNT all > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_scp "MAKE"
      return 1
    fi

    make install > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_scp "MAKE INSTALL"
      return 1
    fi

    cd $ROOTDIR
    pbottom 'INSTALLING FFTW C LIBRARY'
  fi

  # ----------------------------------------------------------------------------
  # ------------------------------- CFITSIO ------------------------------------
  # ----------------------------------------------------------------------------
  if [ -z "${IGNORE_C_CFITSIO_INSTALLATION}" ]; then
    ptop 'INSTALLING CFITSIO C LIBRARY'

    if [ -z "${COCOA_CFITSIO_DIR}" ]; then
      pfail 'COCOA_CFITSIO_DIR'
      cd $ROOTDIR
      return 1
    fi

    cd $ROOTDIR/../cocoa_installation_libraries/$COCOA_CFITSIO_DIR
    if [ $? -ne 0 ]; then
      fail_scp "CD COCOA_CFITSIO_DIR"
      return 1
    fi

    rm -f  CMakeCache.txt
    rm -rf CFITSIOBUILD/

    mkdir ./CFITSIOBUILD
    cd ./CFITSIOBUILD

    $CMAKE -DBUILD_SHARED_LIBS=TRUE \
      -DCMAKE_INSTALL_PREFIX=$ROOTDIR/.local \
      -DCMAKE_C_COMPILER=$C_COMPILER \
      -DCMAKE_CXX_COMPILER=$CXX_COMPILER \
      -DCMAKE_FC_COMPILER=FORTRAN_COMPILER \
      --log-level=ERROR .. > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_scp "CMAKE"
      return 1
    fi

    make -j $C_MNT all > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_scp "MAKE"
      return 1
    fi

    make install > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_scp "MAKE INSTALL"
      return 1
    fi

    cd ../
    rm -rf ./CFITSIOBUILD

    cd $ROOTDIR
    pbottom 'INSTALLING CFITSIO C LIBRARY'
  fi

  # ----------------------------------------------------------------------------
  # ---------------------------------- GSL -------------------------------------
  # ----------------------------------------------------------------------------
  if [ -z "${IGNORE_C_GSL_INSTALLATION}" ]; then
    ptop 'INSTALLING GSL C LIBRARY'

    if [ -z "${COCOA_GSL_DIR}" ]; then
      pfail 'COCOA_GSL_DIR'
      cd $ROOTDIR
      return 1
    fi

    cd $ROOTDIR/../cocoa_installation_libraries/$COCOA_GSL_DIR
    if [ $? -ne 0 ]; then
      fail_scp "CD COCOA_GSL_DIR"
      return 1
    fi

    CC=$C_COMPILER ./configure \
      --prefix=$ROOTDIR/.local \
      --enable-shared=yes \
      --enable-static=yes > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_scp "CONFIGURE"
      return 1
    fi

    make -j $C_MNT all > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_scp "MAKE"
      return 1
    fi

    make install > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_scp "MAKE INSTALL"
      return 1
    fi

    cd $ROOTDIR
    pbottom -e 'INSTALLING GSL C LIBRARY'
  fi

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  unset_env_vars_scp
  pbottom2 'SETUP_C_PACKAGES DONE'
fi

# ------------------------------------------------------------------------------
# ----------------------------------- EUCLID EMU -------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_ALL_PIP_INSTALLATION}" ]; then
  if [ -z "${PIP3}" ]; then
    pfail 'PIP3'
    cd $ROOTDIR
    return 1
  fi
  unset_env_vars_scp_eemul2 () {
    cd $ROOTDIR
    unset OUT1
    unset OUT2
  }
  fail_scp () {
    export MSG="\033[0;31m (setup_c_packages.sh / euclidemul2) WE CANNOT RUN \e[3m"
    export MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"
    unset_env_vars_scp_eemul2
    unset MSG
    unset MSG2
    unset fail_scp
  }
  if [ -z "${DEBUG_PIP_PACKAGES}" ]; then
    export OUT1="/dev/null"
    export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
  fi

  ptop2 'INSTALLING EUCLIDEMU2 DONE'

  # WE MIGRATED euclidemu2 TO setup_c_packages SCRIPT BECAUSE IT DEPENDS ON GSL-GNU LIB
  env CXX=$CXX_COMPILER CC=$C_COMPILER $PIP3 install --global-option=build_ext \
    $ROOTDIR/../cocoa_installation_libraries/euclidemu2-1.2.0 \
    --no-dependencies \
    --prefix=$ROOTDIR/.local --no-index > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_scp "PIP INSTALL EUCLUDEMUL2"
    return 1
  fi

  cd $ROOTDIR
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  unset_env_vars_scp_eemul2
  pbottom2 'INSTALLING EUCLIDEMU2 DONE'
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------