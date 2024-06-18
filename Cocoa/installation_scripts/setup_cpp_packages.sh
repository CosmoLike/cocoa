#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_CPP_INSTALLATION}" ]; then
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
  unset_env_vars_scppp () {
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset CPP_MNT
    unset pfail
    unset unset_env_vars_scppp
  }
  fail_scppp () {
    export MSG="\033[0;31m (setup_cpp_packages.sh) WE CANNOT RUN \e[3m"
    export MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"
    unset_env_vars_scppp
    unset MSG
    unset MSG2
    unset fail_scppp
  }
  if [ -z "${DEBUG_CPP_PACKAGES}" ]; then
    if [ -z "${MAKE_NUM_THREADS}" ]; then
      pfail 'MAKE_NUM_THREADS'
      cd $ROOTDIR
      return 1
    fi
    export OUT1="/dev/null"
    export OUT2="/dev/null"
    export CPP_MNT="${MAKE_NUM_THREADS}"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
    export CPP_MNT=1
  fi

  ptop2 'SETUP_CPP_PACKAGES'

  # ----------------------------------------------------------------------------
  # ---------------------------------- SPDLOG ----------------------------------
  # ----------------------------------------------------------------------------
  if [ -z "${IGNORE_CPP_SPDLOG_INSTALLATION}" ]; then
    ptop 'INSTALLING SPDLOG C++ LIBRARY'

    if [ -z "${COCOA_SPDLOG_DIR}" ]; then
      pfail 'COCOA_SPDLOG_DIR'
      cd $ROOTDIR
      return 1
    fi

    cd $ROOTDIR/../cocoa_installation_libraries/$COCOA_SPDLOG_DIR
    if [ $? -ne 0 ]; then
      fail_scppp "CD COCOA_SPDLOG_DIR"
      return 1
    fi
    
    rm -f CMakeCache.txt

    $CMAKE -DCMAKE_INSTALL_PREFIX=$ROOTDIR/.local \
      -DCMAKE_C_COMPILER=$C_COMPILER \
      -DCMAKE_CXX_COMPILER=$CXX_COMPILER \
      --log-level=ERROR . > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_scppp "CMAKE"
      return 1
    fi

    make -j $CPP_MNT > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_scppp "MAKE"
      return 1
    fi

    make install > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_scppp "MAKE INSTALL"
      return 1
    fi

    cd $ROOTDIR
    pbottom 'INSTALLING SPDLOG C++ LIBRARY'
  fi

  # ----------------------------------------------------------------------------
  # ------------------------------ ARMADILLO -----------------------------------
  # ----------------------------------------------------------------------------
  if [ -z "${IGNORE_CPP_ARMA_INSTALLATION}" ]; then
    ptop 'INSTALLING ARMADILLO C++ LIBRARY'

    if [ -z "${COCOA_ARMADILLO_DIR}" ]; then
      pfail 'COCOA_ARMADILLO_DIR'
      cd $ROOTDIR
      return 1
    fi

    cd $ROOTDIR/../cocoa_installation_libraries/$COCOA_ARMADILLO_DIR
    if [ $? -ne 0 ]; then
      fail_scppp "CD COCOA_ARMADILLO_DIR"
      return 1
    fi

    rm -f CMakeCache.txt

    $CMAKE -DBUILD_SHARED_LIBS=TRUE \
      -DCMAKE_INSTALL_PREFIX=$ROOTDIR/.local \
      -DCMAKE_C_COMPILER=$C_COMPILER \
      -DCMAKE_CXX_COMPILER=$CXX_COMPILER \
      -DLAPACK_FOUND=YES \
      -DLAPACK_LIBRARIES=$ROOTDIR/.local/lib/liblapack.so \
      -DBLAS_FOUND=NO \
      --log-level=ERROR . > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_scppp "CMAKE"
      return 1
    fi

    make clean > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_scppp "MAKE CLEAN"
      return 1
    fi

    make -j $CPP_MNT all -Wno-dev > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_scppp "MAKE"
      return 1
    fi

    make install > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_scppp "MAKE INSTALL"
      return 1
    fi

    cd $ROOTDIR
    pbottom 'INSTALLING ARMADILLO C++ LIBRARY'
  fi

  # ----------------------------------------------------------------------------
  # ---------------------- CARMA (ARMADILLO <-> PYBIND11) ----------------------
  # ----------------------------------------------------------------------------
  if [ -z "${IGNORE_CPP_CARMA_INSTALLATION}" ]; then
    ptop 'INSTALLING CARMA C++ LIBRARY'

    if [ -z "${COCOA_CARMA_DIR}" ]; then
      pfail 'COCOA_CARMA_DIR'
      cd $ROOTDIR
      return 1
    fi

    cd $ROOTDIR/../cocoa_installation_libraries/$COCOA_CARMA_DIR
    if [ $? -ne 0 ]; then
      fail_scppp "CD COCOA_CARMA_DIR"
      return 1
    fi

    rm -rf $ROOTDIR/.local/include/carma/    
    mkdir $ROOTDIR/.local/include/carma/ 
    if [ $? -ne 0 ]; then
      fail_scppp "MKDIR CARMA FOLDER"
      return 1
    fi

    cp ./carma.h $ROOTDIR/.local/include/
    if [ $? -ne 0 ]; then
      fail_scppp "CP CARMA HEADER"
      return 1
    fi

    cp -r ./carma_bits $ROOTDIR/.local/include/
    if [ $? -ne 0 ]; then
      fail_scppp "CP CARMA_BITS"
      return 1
    fi

    cd $ROOTDIR
    pbottom 'INSTALLING CARMA C++ LIBRARY'
  fi

  # ----------------------------------------------------------------------------
  # ---------------------------------- BOOST -----------------------------------
  # ----------------------------------------------------------------------------
  if [ -z "${IGNORE_CPP_BOOST_INSTALLATION}" ]; then
    ptop 'INSTALLING BOOST C++ LIBRARY'

    if [ -z "${COCOA_BOOST_DIR}" ]; then
      pfail 'COCOA_BOOST_DIR'
      cd $ROOTDIR
      return 1
    fi

    cd $ROOTDIR/../cocoa_installation_libraries/$COCOA_BOOST_DIR
    if [ $? -ne 0 ]; then
      fail_scppp "CD COCOA_BOOST_DIR"
      return 1
    fi

    ./bootstrap.sh --prefix=$ROOTDIR/.local > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_scppp "BOOTSTRAP"
      return 1
    fi

    ./b2 --with=regex install --without-python --without-thread \
      --without-timer --without-mpi \
      --without-atomic > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail_scppp "B2 INSTALL"
      return 1
    fi
    
    cd $ROOTDIR
    pbottom 'INSTALLING BOOST C++ LIBRARY'
  fi

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  unset_env_vars_scppp
  pbottom2 'SETUP_CPP_PACKAGES DONE'
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------