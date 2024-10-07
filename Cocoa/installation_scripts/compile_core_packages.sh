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
    unset -v CHANGES TFOLDER TFILE TFILEP AL
    cdroot || return 1;
  }

  unset_env_funcs () {
    unset -f cdfolder cpfolder cpfile error
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

  cpfolder() {
    cp -r "${1:?}" "${2:?}"  \
      2>"/dev/null" || { error "CP FOLDER ${1} on ${2}"; return 1; }
  }
  
  cpfile() {
  cp "${1:?}" "${2:?}" \
    2>"/dev/null" || { error "CP FILE ${1} on ${2}"; return 1; }
  }
  
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  #ptop2 'COMPILE_CORE_PACKAGES' || return 1

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

    PACKAGE_VERSION="${ARMA_VERSION:-"12.8.4"}" 

    DEFAULT="armadillo-${PACKAGE_VERSION:?}"

    FOLDER=${COCOA_ARMADILLO_DIR:-"${DEFAULT:?}"}

    PACKDIR="${CCIL:?}/${FOLDER:?}"

    # --------------------------------------------------------------------------
    # note: in case script run >1x w/ previous run stoped prematurely b/c error

    rm -f "${PACKDIR:?}/CMakeCache.txt"
    rm -rf "${ROOTDIR:?}"/.local/share/Armadillo/
    rm -rf "${ROOTDIR:?}"/.local/lib/libarmadillo*
    rm -rf "${ROOTDIR:?}"/.local/include/armadillo*
    rm -rf "${ROOTDIR:?}"/.local/lib64/libarmadillo*

    # --------------------------------------------------------------------------

    cdfolder "${PACKDIR}" || return 1;

    #g++ XXX -std=c++11 -O2 -L${CONDA_PREFIX}/lib -L${ROOTDIR}/.local/lib -larmadillo -lopenblas -larpack -llapack
    
    ${CMAKE:?} -DBUILD_SHARED_LIBS=TRUE \
      -DCMAKE_INSTALL_PREFIX="${ROOTDIR:?}/.local" \
      -DCMAKE_C_COMPILER="${C_COMPILER:?}" \
      -DCMAKE_CXX_COMPILER="${CXX_COMPILER:?}" \
      -DLAPACK_FOUND=YES \
      -DARPACK_FOUND=YES  \
      -DATLAS_FOUND=NO \
      -DOpenBLAS_FOUND=YES \
      -DBLAS_FOUND=NO  . \
      >${OUT1:?} 2>${OUT2:?} || { error "(ARMA) ${EC12:?}"; return 1; }
    
    make clean \
      >${OUT1:?} 2>${OUT2:?} || { error "(ARMA) ${EC2:?}"; return 1; }

    make all \
      >${OUT1:?} 2>${OUT2:?} || { error "(ARMA) ${EC8:?}"; return 1; }

    make install \
      >${OUT1:?} 2>${OUT2:?} || { error "(ARMA) ${EC10:?}"; return 1; }

    unset -v PACKDIR FOLDER DEFAULT PACKAGE_VERSION

    cdfolder "${ROOTDIR}" || return 1;

    pbottom 'COMPILING ARMADILLO CPP LIBRARY (CORE LIBS)' || return 1

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

    FOLDER=${COCOA_CUBA_DIR:-"${DEFAULT:?}"}

    PACKDIR="${CCIL:?}/${FOLDER:?}"

    # --------------------------------------------------------------------------
    # In case this script runs twice (after being killed by CTRL-D)

    rm -f "${PACKDIR:?}/libcuba.a"
    rm -f "${PACKDIR:?}/libcuba.so"
    rm -f "${PACKDIR:?}/makefile"
    rm -f "${PACKDIR:?}/makefile.patch"
    rm -f "${PACKDIR:?}/config.h"
    rm -f "${ROOTDIR:?}/.local/lib/libcuba.so"
    rm -f "${ROOTDIR:?}/.local/lib/libcuba.a"

    # --------------------------------------------------------------------------
    cdfolder "${PACKDIR}" || return 1;

    FC="${FORTRAN_COMPILER:?}" CC="${C_COMPILER:?}" \
      ./configure \
      --prefix="${ROOTDIR:?}/.local" \
      >${OUT1:?} 2>${OUT2:?} || { error "(CUBA) ${EC11:?}"; return 1; }

    # --------------------------------------------------------------------------
    # Patch CUBA so it also compiles an .so dynamic library --------------------
    # --------------------------------------------------------------------------
    CHANGES="${CCIL:?}/cuba_changes"
    
    declare -a TFOLDER=("" 
                       ) # If nonblank, path must include /
    
    # T = TMP
    declare -a TFILE=("makefile" 
                     )

    #T = TMP, P = PATCH
    declare -a TFILEP=("makefile.patch" 
                      )

    # AL = Array Length
    AL=${#TFOLDER[@]}

    for (( i=0; i<${AL}; i++ ));
    do
      cdfolder "${PACKDIR:?}/${TFOLDER[$i]}" || return 1

      cpfolder "${CHANGES:?}/${TFOLDER[$i]}${TFILEP[$i]:?}" . \
        2>${OUT2:?} || return 1;

      patch -u "${TFILE[$i]:?}" -i "${TFILEP[$i]:?}" >${OUT1:?} \
        2>${OUT2:?} || { error "${EC17:?} (${TFILE[$i]:?})"; return 1; }
    done

    # --------------------------------------------------------------------------
    
    cdfolder "${PACKDIR}" || return 1;

    make clean \
      >${OUT1:?} 2>${OUT2:?} || { error "(CUBA) ${EC2:?}"; return 1; }
      
    make install \
      >${OUT1:?} 2>${OUT2:?} || { error "(CUBA) ${EC10:?}"; return 1; }

    unset -v PACKDIR FOLDER DEFAULT PACKAGE_VERSION 
    unset -v CHANGES TFOLDER TFILE TFILEP AL

    cdfolder "${ROOTDIR}" || return 1;

    pbottom "COMPILING CUBA LIBRARY (CORE LIBS)" || return 1;

  fi

  # ----------------------------------------------------------------------------
  # --------------------------------- LIBEXPAT ---------------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_EXPAT_CORE_PACKAGE}" ]; then
    
    ptop "INSTALLING EXPAT" || return 1

    PACKDIR="${CCIL:?}/${COCOA_EXPAT_DIR:-"expat-2.5.0/"}"

    cdfolder "${PACKDIR:?}" || return 1;
    
    FC="${FORTRAN_COMPILER:?}" CC="${C_COMPILER:?}" ./configure \
      --prefix="${ROOTDIR:?}/.local" \
      --enable-shared=yes \
      --enable-static=yes \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC11:?}"; return 1; }

    make -j $MNT >${OUT1:?} 2>${OUT2:?} || { error "${EC8:?}"; return 1; }

    make install >${OUT1:?} 2>${OUT2:?} || { error "${EC10:?}"; return 1; }

    LLIB="${ROOTDIR:?}/.local/lib"

    cpfile "${LLIB:?}/libexpat.so.1" "${LLIB:?}/libexpat.so.0" || return 1;
    
    cdfolder "${ROOTDIR}" || return 1;

    pbottom "INSTALLING EXPAT" || return 1

  fi

  # ----------------------------------------------------------------------------
  # --------------------------- PIP CORE PACKAGES ------------------------------
  # ----------------------------------------------------------------------------
  
  
  if [ -z "${IGNORE_PIP_CORE_PACKAGES}" ]; then

    ptop "INSTALLING PYTHON CORE LIBRARIES VIA PIP" || return 1

    env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" ${PIP3:?} install \
        'alabaster==0.7.13' \
        'anytree==2.8.0' \
        'appdirs==1.4.4' \
        'astropy==5.2.2' \
        'babel==2.12.1' \
        'cachetools==5.3.1' \
        'certifi==2023.5.7' \
        'charset-normalizer==3.1.0' \
        'configparser==5.3.0' \
        'contourpy==1.1.0' \
        'corner==2.2.1' \
        'coverage==7.5.4' \
        'cycler==0.11.0' \
        'cython==3.0.10' \
        'decorator==5.1.1' \
        'deprecated==1.2.14' \
        'dill==0.3.6' \
        'distlib==0.3.8' \
        'docutils==0.20.1' \
        'filelock==3.13.4' \
        'fonttools==4.40.0' \
        'fuzzywuzzy==0.18.0' \
        'getdist==1.4.8' \
        'gpy==1.10.0' \
        'h5py==3.8.0' \
        'idna==3.4' \
        'imageio==2.31.1' \
        'imagesize==1.4.1' \
        'iminuit==2.25.2' \
        'importlib_metadata==6.6.0' \
        'importlib-resources==5.12.0' \
        'jax==0.4.12' \
        'Jinja2==3.1.2' \
        'joblib==1.4.0' \
        'johnnydep==1.20.2' \
        'kiwisolver==1.4.4' \
        'lazy_loader==0.2' \
        'lenstronomy==1.11.2' \
        'llvmlite==0.40.1' \
        'markupsafe==2.1.3' \
        'matplotlib==3.7.5' \
        'ml-dtypes==0.2.0' \
        'mpmath==1.3.0' \
        'multiprocess==0.70.14' \
        'networkx==3.1' \
        'notebook==7.1.0' \
        'numba==0.57.0' \
        'numpy==1.23.5'  \
        'numpydoc==1.5.0' \
        'opt-einsum==3.3.0' \
        'oauthlib==3.2.2' \
        'oyaml==1.0' \
        'packaging==23.1' \
        'pandas==2.0.3' \
        'paramz==0.9.5' \
        'pgen==0.2.1' \
        'Pillow==9.5.0' \
        'platformdirs==2.6.2' \
        'portalocker==2.7.0' \
        'protobuf==4.23.2' \
        'Py-bobyqa==1.4' \
        'pybind11==2.12.0' \
        'pyDOE2==1.3.0' \
        'pyerfa==2.0.0.3' \
        'pyfftw==0.13.1' \
        'pygments==2.17.2' \
        'pyparsing==3.0.9' \
        'python-dateutil==2.8.2' \
        'pytz==2023.3' \
        'PyWavelets==1.4.1' \
        'pyxdg==0.28' \
        'PyYAML==6.0' \
        'qp-prob==0.8.3' \
        'requests==2.31.0' \
        'sacc==0.8.1' \
        'schwimmbad==0.3.2' \
        'scikit-image==0.21.0' \
        'scikit-learn==1.2.2' \
        'scipy==1.10.1' \
        'setuptools==67.7.2' \
        'setuptools-scm==7.1.0' \
        'six==1.16.0' \
        'snowballstemmer==2.2.0' \
        'sphinx==7.1.2' \
        'sphinxcontrib-applehelp==1.0.4' \
        'sphinxcontrib-devhelp==1.0.2' \
        'sphinxcontrib-htmlhelp==2.0.1' \
        'sphinxcontrib-jsmath==1.0.1' \
        'sphinxcontrib-qthelp==1.0.3' \
        'sphinxcontrib-serializinghtml==1.1.5' \
        'structlog==23.1.0' \
        'sympy==1.12' \
        'syslibrary==0.1' \
        'tables-io==0.8.1' \
        'tabulate==0.9.0' \
        'threadpoolctl==3.1.0' \
        'tifffile==2023.4.12' \
        'tokenizers==0.13.3' \
        'toml==0.10.2' \
        'tomli==2.0.1' \
        'tqdm==4.65.0' \
        'typing_extensions==4.6.3' \
        'tzdata==2023.3' \
        'urllib3==1.26.16' \
        'virtualenv==20.17.1' \
        'wget==3.2' \
        'wheel==0.40.0' \
        'wimpy==0.6' \
        'wrapt==1.14.1' \
        'zipfile38==0.0.3' \
        'zipp==3.15.0' \
      --no-cache-dir \
      --prefix="${ROOTDIR:?}/.local" \
      >${OUT1:?} 2>${OUT2:?} || { error "(CORE-PACKAGES) ${EC13:?}"; return 1; }

    # mpi4py has a weird bug when installing from conda on a few machines 
    # (e.g., midway) no-cache-dir is important to fix this bug
    # https://github.com/mpi4py/mpi4py/issues/335
    env MPICC=$MPI_CC_COMPILER ${PIP3:?} install \
        'mpi4py==3.1.4' \
        'ipyparallel==8.8.0' \
      --no-cache-dir \
      --prefix="${ROOTDIR:?}/.local" \
      --force \
      >${OUT1:?} 2>${OUT2:?} || { error "(PIP-CORE-PACKAGES) ${EC13:?}"; return 1; }

    pbottom "INSTALLING PYTHON CORE LIBRARIES VIA PIP" || return 1

  else
    
    ptop "INSTALLING PYTHON CORE LIBRARIES VIA PIP" || return 1

    #PS: --force-reinstall - this helps CARMA to see numpy files
    env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" ${PIP3:?} install \
        'numpy==1.23.5' \
      --prefix="${ROOTDIR:?}/.local" \
      --force-reinstall \
      >${OUT1:?} 2>${OUT2:?} || { error "(PIP-CORE-PACKAGES) ${EC13:?}"; return 1; }

    # mpi4py has a weird bug when installing from conda on a few machines 
    # (e.g., midway) no-cache-dir is important to fix this bug
    # https://github.com/mpi4py/mpi4py/issues/335
    env MPICC=$MPI_CC_COMPILER ${PIP3:?} install \
        'mpi4py==3.1.4' \
        'notebook==7.1.1' \
        'ipyparallel==8.8.0' \
      --no-cache-dir \
      --prefix="${ROOTDIR:?}/.local" \
      --force \
      >${OUT1:?} 2>${OUT2:?} || { error "(PIP-CORE-PACKAGES) ${EC13:?}"; return 1; }

    pbottom "INSTALLING PYTHON CORE LIBRARIES VIA PIP" || return 1
  
  fi

  # ----------------------------------------------------------------------------
  # ----------------------------- PIP ML PACKAGES ------------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_EMULATOR_CPU_PIP_PACKAGES}" ]; then
  
    if [ -n "${IGNORE_EMULATOR_GPU_PIP_PACKAGES}" ]; then
      error "${EC28:?} (GPU AND CPU EMULATOR FLAGS)"; return 1;
    fi

    ptop "PIP INSTALL MACHINE LEARNING CPU-ONLY PACKAGES" || return 1

    env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" ${PIP3:?} install \
        'tensorflow-cpu==2.12.0' \
        'tensorflow_probability-0.21.0' \
        'keras==2.12.0' \
        'keras-preprocessing==1.1.2' \
        'torch==1.13.1+cpu' \
        'torchvision==0.14.1+cpu' \
        'torchaudio==0.13.1' \
        'tensiometer==0.1.2' \
      --extra-index-url "https://download.pytorch.org/whl/cpu" \
      --prefix="${ROOTDIR:?}/.local" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC13:?}"; return 1; }
  
    pbottom "PIP INSTALL MACHINE LEARNING CPU-ONLY PACKAGES"
  
  fi

  if [ -z "${IGNORE_EMULATOR_GPU_PIP_PACKAGES}" ]; then
  
    if [ -n "${IGNORE_EMULATOR_CPU_PIP_PACKAGES}" ]; then
      error "${EC28:?} (GPU AND CPU EMULATOR FLAGS)"; return 1;
    fi

    ptop "PIP INSTALL MACHINE LEARNING GPU PACKAGES"

    env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" ${PIP3:?} install \
        'tensorflow==2.12.0' \
        'tensorflow_probability-0.21.0' \
        'keras==2.12.0' \
        'keras-preprocessing==1.1.2' \
        'torch==1.13.1+cu116' \
        'torchvision==0.14.1+cu116' \
        'torchaudio==0.13.1' \
        'tensiometer==0.1.2' \
      --extra-index-url "https://download.pytorch.org/whl/cu116" \
      --prefix="${ROOTDIR:?}/.local" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC13:?}"; return 1; } 

    pbottom "PIP INSTALL MACHINE LEARNING GPU PACKAGES" || return 1
  
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

  unset_all || return 1; 

  #pbottom2 'COMPILE_CORE_PACKAGES' || return 1

fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
