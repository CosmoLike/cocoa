#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_CORE_INSTALLATION}" ]; then

  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi

  # Parenthesis = run in a subshell
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;
    
  unset_env_vars () {
    unset -v CCIL BDF PACKDIR FOLDER MAKE_NB_JOBS PACKAGE_VERSION DEFAULT
    cdroot || return 1;
  }

  unset_env_funcs () {
    unset -f cdfolder cpfolder error
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
    fail_script_msg "$(basename "${BASH_SOURCE[0]}")" "${1}"
    unset_all || return 1
  }
    
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
  }
  
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  ptop2 'COMPILE_CORE_PACKAGES' || return 1

  unset_env_vars || return 1; 

  CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

  # ----------------------------------------------------------------------------
  # ----------------------------- CMAKE LIBRARY --------------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_CMAKE_INSTALLATION}" ]; then
    
    ptop 'COMPILING CMAKE LIBRARY (CORE LIBS)' || return 1

    PACKAGE_VERSION="${CMAKE_VERSION:-"3.26.4"}"

    DEFAULT="cmake-${PACKAGE_VERSION:?}/"

    FOLDER=${COCOA_CMAKE_DIR:-"${DEFAULT:?}"}

    PACKDIR="${CCIL:?}/${FOLDER:?}"

    cdfolder "${PACKDIR:?}" || return 1;

    env CC="${C_COMPILER:?}" CXX="${CXX_COMPILER:?}" \
      ./bootstrap --prefix="${ROOTDIR:?}/.local" \
      >${OUT1:?} 2>${OUT2:?} || { error "(CMAKE) ${EC19:?}"; return 1; }

    make -j $MNT \
      >${OUT1:?} 2>${OUT2:?} || { error "(CMAKE) ${EC8:?}"; return 1; }

    make install \
      >${OUT1:?} 2>${OUT2:?} || { error "(CMAKE) ${EC10:?}"; return 1; }

    unset -v PACKDIR FOLDER PACKAGE_VERSION DEFAULT

    cdfolder "${ROOTDIR}" || return 1;

    pbottom 'COMPILING CMAKE LIBRARY (CORE LIBS)' || return 1
  
  fi

  # ----------------------------------------------------------------------------
  # ------------------------------ WGET LIBRARY --------------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_WGET_INSTALLATION}" ]; then

    pbottom "COMPILING WGET LIBRARY (CORE LIBS)" || return 1;
    
    PACKAGE_VERSION="${WGET_VERSION:-"1.24.5"}"

    DEFAULT="wget-${PACKAGE_VERSION:?}"

    FOLDER=${COCOA_WGET_DIR:-"${DEFAULT:?}"}

    PACKDIR="${CCIL:?}/${FOLDER:?}"

    cdfolder "${PACKDIR}" || return 1;

    FC="${FORTRAN_COMPILER:?}" CC="${C_COMPILER:?}" \
      ./configure \
      --prefix="${ROOTDIR:?}/.local" \
      >${OUT1:?} 2>${OUT2:?} || { error "(WGET) ${EC11:?}"; return 1; }

    # --------------------------------------------------------------------------
    # note: in case script run >1x w/ previous run stoped prematurely b/c error
    
    make clean \
      >${OUT1:?} 2>${OUT2:?} || { error "(WGET) ${EC2:?}"; return 1; }

    # --------------------------------------------------------------------------
    
    make -j $MNT all \
      >${OUT1:?} 2>${OUT2:?} || { error "(WGET) ${EC7:?}"; return 1; }
      
    make install \
      >${OUT1:?} 2>${OUT2:?} || { error "(WGET) ${EC10:?}"; return 1; }

    unset -v PACKDIR FOLDER DEFAULT PACKAGE_VERSION

    cdfolder "${ROOTDIR}" || return 1;

    pbottom "COMPILING WGET LIBRARY (CORE LIBS)" || return 1;

  fi

  # ----------------------------------------------------------------------------
  # ----------------------- TEXINFO LIBRARY (DISTUTILS) ------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_DISTUTILS_INSTALLATION}" ]; then

    ptop 'COMPILING TEXINFO LIBRARY (CORE LIBS)' || return 1;
    
    PACKAGE_VERSION="${TEXTINFO_VERSION:-"7.0.3"}"

    DEFAULT="texinfo-${PACKAGE_VERSION:?}"

    FOLDER=${COCOA_TEXINFO_DIR:-"${DEFAULT:?}"}

    PACKDIR="${CCIL:?}/${FOLDER:?}"

    cdfolder "${PACKDIR}" || return 1;

    FC="${FORTRAN_COMPILER:?}" CC="${C_COMPILER:?}" \
      ./configure \
      --prefix="${ROOTDIR:?}/.local" \
      --disable-perl-xs \
      >${OUT1:?} 2>${OUT2:?} || { error "(TEXINFO) ${EC11:?}"; return 1; }

    # --------------------------------------------------------------------------
    # note: in case script run >1x w/ previous run stoped prematurely b/c error
    
    make clean \
      >${OUT1:?} 2>${OUT2:?} || { error "(TEXINFO) ${EC2:?}"; return 1; }

    # --------------------------------------------------------------------------

    make -j $MNT all \
      >${OUT1:?} 2>${OUT2:?} || { error "(TEXINFO) ${EC7:?}"; return 1; }
      
    make install \
      >${OUT1:?} 2>${OUT2:?} || { error "(TEXINFO) ${EC10:?}"; return 1; }

    unset -v PACKDIR FOLDER DEFAULT PACKAGE_VERSION

    cdfolder "${ROOTDIR}" || return 1;
    
    pbottom 'COMPILING TEXINFO LIBRARY (CORE LIBS)' || return 1;
  
  fi

  # ----------------------------------------------------------------------------
  # ---------------------------- BINUTILS LIBRARY  -----------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_DISTUTILS_INSTALLATION}" ]; then
  
    ptop 'COMPILING BINUTILS LIBRARY (CORE LIBS)' || return 1;
    
    PACKAGE_VERSION="${BINUTILS_VERSION:-"2.37"}"

    DEFAULT="binutils-${PACKAGE_VERSION:?}"

    FOLDER=${COCOA_BINUTILS_DIR:-"${DEFAULT:?}"}

    PACKDIR="${CCIL:?}/${FOLDER:?}"

    cdfolder "${PACKDIR}" || return 1;

    FC="${FORTRAN_COMPILER:?}" CC="${C_COMPILER:?}" \
      ./configure \
      --prefix="${ROOTDIR:?}/.local" \
      >${OUT1:?} 2>${OUT2:?} || { error "(BINUTILS) ${EC11:?}"; return 1; }

    # --------------------------------------------------------------------------
    # note: in case script run >1x w/ previous run stoped prematurely b/c error
    
    make clean \
      >${OUT1:?} 2>${OUT2:?} || { error "(BINUTILS) ${EC2:?}"; return 1; }

    # --------------------------------------------------------------------------

    make -j $MNT \
      >${OUT1:?} 2>${OUT2:?} || { error "(BINUTILS) ${EC7:?}"; return 1; }

    make install \
      >${OUT1:?} 2>${OUT2:?} || { error "(BINUTILS) ${EC10:?}"; return 1; }

    unset -v PACKDIR FOLDER DEFAULT PACKAGE_VERSION

    cdfolder "${ROOTDIR}" || return 1;

    pbottom 'COMPILING BINUTILS LIBRARY (CORE LIBS)' || return 1;
  
  fi

  # ----------------------------------------------------------------------------
  # ------------------------------ HDF5 LIBRARY --------------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_HDF5_INSTALLATION}" ]; then
  
    ptop 'COMPILING HFD5 LIBRARY (CORE LIBS)' || return 1;
  
    FOLDER=${COCOA_HDF5_DIR:-"hdf5-1.12.3/"}

    PACKDIR="${CCIL:?}/${FOLDER:?}"
    
    BDF="${PACKDIR:?}/cocoa_HDF5_build"

    # --------------------------------------------------------------------------
    # note: in case script run >1x w/ previous run stoped prematurely b/c error

    rm -f  "${PACKDIR:?}/CMakeCache.txt"
    rm -rf "${PACKDIR:?}/CMakeFiles/"
    rm -rf "${BDF:?}"
    
    # --------------------------------------------------------------------------

    mkdir "${BDF:?}" 2>${OUT2:?} || { error "${EC20:?}"; return 1; }

    cdfolder "${BDF:?}" || return 1;

    "${CMAKE:?}" -DBUILD_SHARED_LIBS=TRUE \
      -DCMAKE_INSTALL_PREFIX="${ROOTDIR:?}/.local" \
      -DCMAKE_C_COMPILER="${C_COMPILER:?}" \
      -DCMAKE_CXX_COMPILER="${CXX_COMPILER:?}" \
      -DCMAKE_FC_COMPILER="${FORTRAN_COMPILER:?}" \
      --log-level=ERROR .. \
      >${OUT1:?} 2>${OUT2:?} || { error "(HDF5) ${EC12:?}"; return 1; }

    make -j $MNT \
      >${OUT1:?} 2>${OUT2:?} || { error "(HDF5) ${EC8:?}"; return 1; }

    make install \
      >${OUT1:?} 2>${OUT2:?} || { error "(HDF5) ${EC10:?}"; return 1; }

    unset -v PACKDIR BDF
    
    cdfolder "${ROOTDIR}" || return 1;

    pbottom 'COMPILING HDF5 LIBRARY DONE (CORE LIBS)' || return 1;
  
  fi

  # ----------------------------------------------------------------------------
  # -------------------------- OPENBLAS LIBRARY --------------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_OPENBLAS_INSTALLATION}" ]; then
    
    ptop  'COMPILING OPENBLAS LIBRARY (CORE LIBS)' || return 1;

    PACKAGE_VERSION="${OPENBLAS_VERSION:-"0.3.23"}" 

    DEFAULT="OpenBLAS-${PACKAGE_VERSION:?}"

    FOLDER=${COCOA_OPENBLAS_DIR:-"${DEFAULT:?}"}

    PACKDIR="${CCIL:?}/${FOLDER:?}"

    cdfolder "${PACKDIR:?}" || return 1;

    export MAKE_NB_JOBS=$MNT
    
    # --------------------------------------------------------------------------
    # note: in case script run >1x w/ previous run stoped prematurely b/c error

    make clean \
      >${OUT1:?} 2>${OUT2:?} || { error "(OpenBLAS) ${EC2:?}"; return 1; }

    # --------------------------------------------------------------------------

    make CC="${C_COMPILER:?}" FC="${FORTRAN_COMPILER:?}" USE_OPENMP=1 \
      >${OUT1:?} 2>${OUT2:?} || { error "(OpenBLAS) ${EC8:?}"; return 1; }
    
    make install PREFIX="${ROOTDIR:?}/.local" \
      >${OUT1:?} 2>${OUT2:?} || { error "(OpenBLAS) ${EC10:?}"; return 1; }

    unset -v PACKDIR  MAKE_NB_JOBS FOLDER DEFAULT PACKAGE_VERSION

    cdfolder "${ROOTDIR:?}" || return 1;

    pbottom  'COMPILING OPENBLAS LIBRARY (CORE LIBS)' || return 1;
  
  fi
  
  # ----------------------------------------------------------------------------
  # ---------------------------- LAPACK LIBRARY --------------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_FORTRAN_LAPACK_INSTALLATION}" ]; then
    
    ptop 'COMPILING LAPACK FORTRAN LIBRARY (CORE LIBS)' || return 1;

    PACKAGE_VERSION="${LAPACK_VERSION:-"3.11"}" 

    DEFAULT="lapack-${PACKAGE_VERSION:?}"

    FOLDER=${COCOA_LAPACK_DIR:-"${DEFAULT:?}"}

    PACKDIR="${CCIL:?}/${FOLDER:?}"
    
    BDF="${PACKDIR:?}/lapack-build"

    # --------------------------------------------------------------------------
    # note: in case script run >1x w/ previous run stoped prematurely b/c error

    rm -rf "${BDF:?}"

    # --------------------------------------------------------------------------
    
    mkdir "${BDF:?}" 2>${OUT2:?} || { error "${EC20:?}"; return 1; }

    cdfolder "${BDF:?}" || return 1;

    ${CMAKE:?} -DBUILD_SHARED_LIBS=TRUE \
      -DCMAKE_INSTALL_PREFIX="${ROOTDIR:?}/.local" \
      -DCMAKE_C_COMPILER="${C_COMPILER:?}" \
      --log-level=ERROR "../${PACKDIR:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error "(LAPACK) ${EC12:?}"; return 1; }

    make -j $MNT all \
      >${OUT1:?} 2>${OUT2:?} || { error "(LAPACK) ${EC8:?}"; return 1; }

    make install \
      >${OUT1:?} 2>${OUT2:?} || { error "(LAPACK) ${EC10:?}"; return 1; }

    # --------------------------------------------------------------------------
    # clean build

    rm -rf "${BDF:?}"
    
    # --------------------------------------------------------------------------

    unset -v PACKDIR BDF FOLDER PACKAGE_VERSION DEFAULT

    cdfolder "${ROOTDIR}" || return 1;

    pbottom 'COMPILING LAPACK FORTRAN LIBRARY (CORE LIBS)' || return 1;

  fi
  
  # ----------------------------------------------------------------------------
  # -------------------------------- FFTW --------------------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_C_FFTW_INSTALLATION}" ]; then
    
    ptop 'COMPILING FFTW C LIBRARY (CORE LIBS)' || return 1

    PACKAGE_VERSION="${FFTW_VERSION:-"3.3.10"}" 

    DEFAULT="fftw-${PACKAGE_VERSION:?}"

    FOLDER=${COCOA_FFTW_DIR:-"${DEFAULT:?}"}

    PACKDIR="${CCIL:?}/${FOLDER:?}"

    cdfolder "${PACKDIR}" || return 1;

    FC="${FORTRAN_COMPILER:?}" CC="${C_COMPILER:?}" ./configure \
      --enable-openmp \
      --prefix="${ROOTDIR:?}/.local" \
      --enable-shared=yes \
      --enable-static=yes \
      >${OUT1:?} 2>${OUT2:?} || { error "(FFTW) ${EC11:?}"; return 1; }

    make -j $MNT all \
      >${OUT1:?} 2>${OUT2:?} || { error "(FFTW) ${EC8:?}"; return 1; }

    make install \
      >${OUT1:?} 2>${OUT2:?} || { error "(FFTW) ${EC10:?}"; return 1; }

    unset -v PACKDIR FOLDER PACKAGE_VERSION DEFAULT

    cdfolder "${ROOTDIR}" || return 1;

    pbottom 'COMPILING FFTW C LIBRARY (CORE LIBS)' || return 1

  fi

  # ----------------------------------------------------------------------------
  # ------------------------------- CFITSIO ------------------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_C_CFITSIO_INSTALLATION}" ]; then
    
    ptop 'COMPILING CFITSIO C LIBRARY (CORE LIBS)' || return 1

    PACKAGE_VERSION="${CFITSIO_VERSION:-"4.0.0"}" 

    DEFAULT="cfitsio-${PACKAGE_VERSION:?}"

    FOLDER=${COCOA_CFITSIO_DIR:-"${DEFAULT:?}"}

    PACKDIR="${CCIL:?}/${FOLDER:?}"
    
    BDF="${PACKDIR:?}/CFITSIOBUILD"

    # --------------------------------------------------------------------------
    # note: in case script run >1x w/ previous run stoped prematurely b/c error

    rm -f  "${PACKDIR:?}/CMakeCache.txt"
    rm -rf "${BDF:?}"
    
    # --------------------------------------------------------------------------

    mkdir "${BDF:?}/" 2>${OUT2:?} || { error "${EC14:?}"; return 1; }
    
    cdfolder "${BDF:?}" || return 1;

    ${CMAKE:?} -DBUILD_SHARED_LIBS=TRUE \
      -DCMAKE_INSTALL_PREFIX="${ROOTDIR:?}/.local" \
      -DCMAKE_C_COMPILER="${C_COMPILER:?}" \
      -DCMAKE_CXX_COMPILER="${CXX_COMPILER:?}" \
      -DCMAKE_FC_COMPILER="${FORTRAN_COMPILER:?}" \
      --log-level=ERROR .. \
      >${OUT1:?} 2>${OUT2:?} || { error "(CFITSIO) ${EC12:?}"; return 1; }

    make -j $MNT all \
      >${OUT1:?} 2>${OUT2:?} || { error "(CFITSIO) ${EC8:?}"; return 1; }

    make install \
      >${OUT1:?} 2>${OUT2:?} || { error "(CFITSIO) ${EC10:?}"; return 1; }
    
    # --------------------------------------------------------------------------
    # clean build

    rm -rf "${BDF:?}"
    
    # --------------------------------------------------------------------------

    unset -v PACKDIR BDF FOLDER DEFAULT PACKAGE_VERSION

    cdfolder "${ROOTDIR}" || return 1;

    pbottom 'COMPILING CFITSIO C LIBRARY (CORE LIBS)' || return 1

  fi

  # ----------------------------------------------------------------------------
  # ---------------------------------- GSL -------------------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_C_GSL_INSTALLATION}" ]; then
    
    ptop 'COMPILING GSL C LIBRARY (CORE LIBS)' || return 1

    PACKAGE_VERSION="${GSL_VERSION:-"2.7"}" 

    DEFAULT="gsl-${PACKAGE_VERSION:?}"

    FOLDER=${COCOA_GSL_DIR:-"${DEFAULT:?}"}

    PACKDIR="${CCIL:?}/${FOLDER:?}" 

    cdfolder "${PACKDIR}" || return 1;

    CC="${C_COMPILER:?}" ./configure \
      --prefix="${ROOTDIR:?}/.local" \
      --enable-shared=yes \
      --enable-static=yes \
      >${OUT1:?} 2>${OUT2:?} || { error "(GSL) ${EC11:?}"; return 1; }
 
    make -j $MNT all \
      >${OUT1:?} 2>${OUT2:?} || { error "(GSL) ${EC8:?}"; return 1; }

    make install \
      >${OUT1:?} 2>${OUT2:?} || { error "(GSL) ${EC10:?}"; return 1; }

    unset -v PACKDIR FOLDER DEFAULT PACKAGE_VERSION

    cdfolder "${ROOTDIR}" || return 1;

    pbottom 'COMPILING GSL C LIBRARY (CORE LIBS)' || return 1

  fi

  # ----------------------------------------------------------------------------
  # ----------------------------------- EUCLID EMU -----------------------------
  # ----------------------------------------------------------------------------
  #note: we migrated euclidemu2 to setup_c_packages. 
  #note: Why? it depends on GSL-GNU library
  
  if [ -z "${IGNORE_ALL_PIP_INSTALLATION}" ]; then
            
    ptop 'COMPILING EUCLIDEMU2' || return 1

    env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" ${PIP3:?} install \
      --global-option=build_ext "${CCIL:?}/euclidemu2-1.2.0" \
      --no-dependencies \
      --prefix="${ROOTDIR:?}/.local" \
      --no-index 
      >${OUT1:?} 2>${OUT2:?} || { error "(EUCLIDEMU2) ${EC13:?}"; return 1; }

    cdfolder "${ROOTDIR}" || return 1;  
    
    pbottom 'COMPILING EUCLIDEMU2' || return 1
  
  fi

  # ----------------------------------------------------------------------------
  # ---------------------------------- SPDLOG ----------------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_CPP_SPDLOG_INSTALLATION}" ]; then
    
    ptop 'COMPILING SPDLOG CPP LIBRARY (CORE LIBS)' || return 1

    FOLDER=${COCOA_SPDLOG_DIR:-"spdlog/"}

    PACKDIR="${CCIL:?}/${FOLDER:?}"

    # --------------------------------------------------------------------------
    # note: in case script run >1x w/ previous run stoped prematurely b/c error

    rm -f "${PACKDIR:?}/CMakeCache.txt"

    # --------------------------------------------------------------------------

    cdfolder "${PACKDIR}" || return 1;
    
    ${CMAKE:?} -DCMAKE_INSTALL_PREFIX="${ROOTDIR:?}/.local" \
      -DCMAKE_C_COMPILER="${C_COMPILER:?}" \
      -DCMAKE_CXX_COMPILER="${CXX_COMPILER:?}" \
      --log-level=ERROR . \
      >${OUT1:?} 2>${OUT2:?} || { error "(SPDLOG) ${EC12:?}"; return 1; }

    make -j $MNT \
      >${OUT1:?} 2>${OUT2:?} || { error "(SPDLOG) ${EC8:?}"; return 1; }

    make install \
      >${OUT1:?} 2>${OUT2:?} || { error "(SPDLOG) ${EC10:?}"; return 1; }

    unset -v PACKDIR FOLDER

    cdfolder "${ROOTDIR}" || return 1;

    pbottom 'COMPILING SPDLOG CPP LIBRARY (CORE LIBS)' || return 1

  fi

  # ----------------------------------------------------------------------------
  # ------------------------------ ARMADILLO -----------------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_CPP_ARMA_INSTALLATION}" ]; then
    
    ptop 'COMPILING ARMADILLO CPP LIBRARY (CORE LIBS)' || return 1

    PACKAGE_VERSION="${ARMA_VERSION:-"12.8.2"}" 

    DEFAULT="armadillo-${PACKAGE_VERSION:?}"

    FOLDER=${COCOA_ARMADILLO_DIR:-"${DEFAULT:?}"}

    PACKDIR="${CCIL:?}/${FOLDER:?}"

    # --------------------------------------------------------------------------
    # note: in case script run >1x w/ previous run stoped prematurely b/c error

    rm -f "${PACKDIR:?}/CMakeCache.txt"

    # --------------------------------------------------------------------------

    cdfolder "${PACKDIR}" || return 1;

    ${CMAKE:?} -DBUILD_SHARED_LIBS=TRUE \
      -DCMAKE_INSTALL_PREFIX="${ROOTDIR:?}/.local" \
      -DCMAKE_C_COMPILER="${C_COMPILER:?}" \
      -DCMAKE_CXX_COMPILER="${CXX_COMPILER:?}" \
      -DLAPACK_FOUND=YES \
      -DLAPACK_LIBRARIES="${ROOTDIR:?}/.local/lib/liblapack.so" \
      -DBLAS_FOUND=NO \
      --log-level=ERROR . \
      >${OUT1:?} 2>${OUT2:?} || { error "(ARMA) ${EC12:?}"; return 1; }
    
    make clean \
      >${OUT1:?} 2>${OUT2:?} || { error "(ARMA) ${EC2:?}"; return 1; }

    make -j $MNT all -Wno-dev \
      >${OUT1:?} 2>${OUT2:?} || { error "(ARMA) ${EC8:?}"; return 1; }

    make install \
      >${OUT1:?} 2>${OUT2:?} || { error "(ARMA) ${EC10:?}"; return 1; }

    unset -v PACKDIR FOLDER DEFAULT PACKAGE_VERSION

    cdfolder "${ROOTDIR}" || return 1;

    pbottom 'COMPILING ARMADILLO CPP LIBRARY (CORE LIBS)' || return 1

  fi

  # ----------------------------------------------------------------------------
  # ---------------------- CARMA (ARMADILLO <-> PYBIND11) ----------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_CPP_CARMA_INSTALLATION}" ]; then
    
    ptop 'COMPILING CARMA CPP LIBRARY (CORE LIBS)' || return 1

    FOLDER=${COCOA_CARMA_DIR:-"carma/"}

    PACKDIR="${CCIL:?}/${FOLDER:?}"

    # --------------------------------------------------------------------------
    # note: in case script run >1x w/ previous run stoped prematurely b/c error

    rm -rf "${ROOTDIR:?}/.local/include/carma/"   

    # --------------------------------------------------------------------------

    mkdir "${ROOTDIR:?}/.local/include/carma/" \
      2>${OUT2:?} || { error "${EC20:?}"; return 1; }
    
    cpfile "${PACKDIR:?}/carma.h" "${ROOTDIR:?}/.local/include/" || return 1

    cpfolder "${PACKDIR:?}/carma_bits" "${ROOTDIR:?}/.local/include/" || return 1

    unset -v PACKDIR FOLDER
    
    cdfolder "${ROOTDIR}" || return 1;

    pbottom 'COMPILING CARMA CPP LIBRARY (CORE LIBS)' || return 1

  fi

  # ----------------------------------------------------------------------------
  # ---------------------------------- BOOST -----------------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_CPP_BOOST_INSTALLATION}" ]; then

    ptop 'COMPILING BOOST CPP LIBRARY (CORE LIBS)' || return 1

    PACKAGE_VERSION="${CPP_BOOST_VERSION:-"81"}"

    DEFAULT="boost_1_${PACKAGE_VERSION:?}_0"

    FOLDER=${COCOA_BOOST_DIR:-"${DEFAULT:?}"}

    PACKDIR="${CCIL:?}/${FOLDER:?}"

    cdfolder "${PACKDIR}" || return 1;

    ./bootstrap.sh --prefix="${ROOTDIR:?}/.local" \
      >${OUT1:?} 2>${OUT2:?} || { error "(BOOST) ${EC19:?}"; return 1; }

    ./b2 --with=regex install \
      --without-python \
      --without-thread \
      --without-timer  \
      --without-mpi \
      --without-atomic \
      >${OUT1:?} 2>${OUT2:?} || { error "(BOOST) ${EC21:?}"; return 1; }
    
    unset -v PACKDIR FOLDER DEFAULT PACKAGE_VERSION

    cdfolder "${ROOTDIR}" || return 1;

    pbottom 'COMPILING BOOST CPP LIBRARY (CORE LIBS)' || return 1

  fi

  # ----------------------------------------------------------------------------
  # ------------------------------ CUBA LIBRARY --------------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_CPP_CUBA_INSTALLATION}" ]; then

    pbottom "COMPILING CUBA LIBRARY (CORE LIBS)" || return 1;

    PACKAGE_VERSION="${CPP_CUBA_VERSION:-"4.2.2"}"

    DEFAULT="Cuba-${PACKAGE_VERSION:?}"

    FOLDER=${COCOA_WGET_DIR:-"${DEFAULT:?}"}

    PACKDIR="${CCIL:?}/${FOLDER:?}"

    cdfolder "${PACKDIR}" || return 1;

    FC="${FORTRAN_COMPILER:?}" CC="${C_COMPILER:?}" \
      ./configure \
      --prefix="${ROOTDIR:?}/.local" \
      >${OUT1:?} 2>${OUT2:?} || { error "(CUBA) ${EC11:?}"; return 1; }

    # --------------------------------------------------------------------------
    # note: in case script run >1x w/ previous run stoped prematurely b/c error
    
    make clean \
      >${OUT1:?} 2>${OUT2:?} || { error "(CUBA) ${EC2:?}"; return 1; }

    # --------------------------------------------------------------------------
    
    make -j $MNT all \
      >${OUT1:?} 2>${OUT2:?} || { error "(CUBA) ${EC7:?}"; return 1; }
      
    make install \
      >${OUT1:?} 2>${OUT2:?} || { error "(CUBA) ${EC10:?}"; return 1; }

    unset -v PACKDIR FOLDER DEFAULT PACKAGE_VERSION

    cdfolder "${ROOTDIR}" || return 1;

    pbottom "COMPILING CUBA LIBRARY (CORE LIBS)" || return 1;

  fi

  # ----------------------------------------------------------------------------

  unset_all || return 1; 

  pbottom2 'COMPILE_CORE_PACKAGES' || return 1

fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------