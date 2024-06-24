#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_CPP_INSTALLATION}" ]; then
  
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1;
  fi

  source "${ROOTDIR:?}/installation_scripts/.check_flags.sh" || return 1;
    
  unset_env_vars () {
    unset -v CCIL PACKDIR
    cdroot || return 1;
  }
  
  unset_env_funcs () {
    unset -f cdfolder cpfolder error cpfile
    unset -f unset_env_funcs
    cdroot || return 1;
  }

  unset_all () {
    unset_env_vars
    unset_env_funcs
    unset -f unset_all
    cdroot || return 1;
  }

  error () {
    fail_script_msg "setup_cpp_packages.sh" "${1}"
    unset_all || return 1
  }
  
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
  }

  cpfolder() {
    cp -r "${1:?}" "${2:?}" \
      2>"/dev/null" || { error "CP FOLDER ${1} on ${2}"; return 1; }
  }

 cpfile() {
    cp "${1:?}" "${2:?}" \
      2>"/dev/null" || { error "CP FILE ${1} on ${2}"; return 1; }
  }
  
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  ptop2 'SETUP_CPP_PACKAGES' || return 1

  CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

  # ----------------------------------------------------------------------------
  # ---------------------------------- SPDLOG ----------------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_CPP_SPDLOG_INSTALLATION}" ]; then
    
    ptop 'INSTALLING SPDLOG C++ LIBRARY' || return 1

    PACKDIR="${CCIL:?}/${COCOA_SPDLOG_DIR:-"spdlog/"}"

    rm -f "${PACKDIR:?}/CMakeCache.txt"

    cdfolder "${PACKDIR}" || return 1;
    
    ${CMAKE:?} -DCMAKE_INSTALL_PREFIX="${ROOTDIR:?}/.local" \
      -DCMAKE_C_COMPILER="${C_COMPILER:?}" \
      -DCMAKE_CXX_COMPILER="${CXX_COMPILER:?}" \
      --log-level=ERROR . \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC12:?}"; return 1; }

    make -j $MNT >${OUT1:?} 2>${OUT2:?} || { error "${EC8:?}"; return 1; }

    make install >${OUT1:?} 2>${OUT2:?} || { error "${EC10:?}"; return 1; }

    cdfolder "${ROOTDIR}" || return 1;

    pbottom 'INSTALLING SPDLOG C++ LIBRARY' || return 1

  fi

  # ----------------------------------------------------------------------------
  # ------------------------------ ARMADILLO -----------------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_CPP_ARMA_INSTALLATION}" ]; then
    
    ptop 'INSTALLING ARMADILLO C++ LIBRARY' || return 1

    PACKDIR="${CCIL:?}/${COCOA_ARMADILLO_DIR:-"armadillo-12.8.2/"}"

    rm -f "${PACKDIR:?}/CMakeCache.txt"

    cdfolder "${PACKDIR}" || return 1;

    ${CMAKE:?} -DBUILD_SHARED_LIBS=TRUE \
      -DCMAKE_INSTALL_PREFIX="${ROOTDIR:?}/.local" \
      -DCMAKE_C_COMPILER="${C_COMPILER:?}" \
      -DCMAKE_CXX_COMPILER="${CXX_COMPILER:?}" \
      -DLAPACK_FOUND=YES \
      -DLAPACK_LIBRARIES="${ROOTDIR:?}/.local/lib/liblapack.so" \
      -DBLAS_FOUND=NO \
      --log-level=ERROR . \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC12:?}"; return 1; }
    
    make clean >${OUT1:?} 2>${OUT2:?} || { error "${EC2:?}"; return 1; }

    make -j $MNT all -Wno-dev \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC8:?}"; return 1; }

    make install >${OUT1:?} 2>${OUT2:?} || { error "${EC10:?}"; return 1; }

    cdfolder "${ROOTDIR}" || return 1;

    pbottom 'INSTALLING ARMADILLO C++ LIBRARY' || return 1

  fi

  # ----------------------------------------------------------------------------
  # ---------------------- CARMA (ARMADILLO <-> PYBIND11) ----------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_CPP_CARMA_INSTALLATION}" ]; then
    
    ptop 'INSTALLING CARMA C++ LIBRARY' || return 1

    PACKDIR="${CCIL:?}/${COCOA_CARMA_DIR:-"carma/"}"

    rm -rf "${ROOTDIR:?}/.local/include/carma/"   
    
    mkdir "${ROOTDIR:?}/.local/include/carma/" \
      2>${OUT2:?} || { error "${EC20:?}"; return 1; }
    
    cpfile "${PACKDIR:?}/carma.h" "${ROOTDIR:?}/.local/include/" || return 1

    cpfolder "${PACKDIR:?}/carma_bits" "${ROOTDIR:?}/.local/include/" || return 1

    pbottom 'INSTALLING CARMA C++ LIBRARY' || return 1

  fi

  # ----------------------------------------------------------------------------
  # ---------------------------------- BOOST -----------------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_CPP_BOOST_INSTALLATION}" ]; then

    ptop 'INSTALLING BOOST C++ LIBRARY' || return 1

    PACKDIR="${CCIL:?}/${COCOA_BOOST_DIR:-"boost_1_81_0/"}"

    cdfolder "${PACKDIR}" || return 1;

    ./bootstrap.sh --prefix="${ROOTDIR:?}/.local" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC19:?}"; return 1; }

    ./b2 --with=regex install \
      --without-python \
      --without-thread \
      --without-timer  \
      --without-mpi \
      --without-atomic \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC21:?}"; return 1; }
    
    cdfolder "${ROOTDIR}" || return 1;

    pbottom 'INSTALLING BOOST C++ LIBRARY' || return 1

  fi

  # ----------------------------------------------------------------------------
 
  unset_all || return 1; 
  
  pbottom2 'SETUP_CPP_PACKAGES DONE' || return 1

fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------