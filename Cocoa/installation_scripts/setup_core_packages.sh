#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_CORE_INSTALLATION}" ]; then

  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1
  fi

  # parenthesis = run in a subshell
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" )  || return 1;

  unset_env_vars () {
    unset -v URL_BASE URL FOLDER VER XZF CCIL CNAME PACKDIR PACKAGE_VERSION 
    cdroot || return 1;
  }

  unset_env_funcs () {
    unset -f cdfolder cpfolder cpfile error wgetact gitact gitact1 gitact2
    unset -f unset_env_funcs wgetact1 wgetact2
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

  wgetact1() {
    #ARGUMENTS: FOLDER, FILE, URL, NEW_FOLDER, XZFILE
    
    local CCIL="${ROOTDIR}/../cocoa_installation_libraries"

    cdfolder "${CCIL:?}" || return 1;

    local FILE="${1:?}.${2:?}"

    # In case this script runs twice -------------------------------------------
    if [ -n "${OVERWRITE_EXISTING_CORE_PACKAGES}" ]; then
      rm -rf "${CCIL:?}/${4:?}"
      rm -rf "${CCIL:?}/${1:?}"
      if [ -n "${REDOWNLOAD_EXISTING_CORE_PACKAGES}" ]; then
        rm -f  "${CCIL:?}/${FILE:?}"
      fi 
    fi

    if [[ ! -d "${CCIL:?}/${4:?}" && ! -d "${CCIL:?}/${1:?}" ]]; then

      if [[ ! -e "${CCIL:?}/${FILE:?}" ]]; then
        wget "${3:?}" --retry-connrefused --waitretry=1 --tries=3 --read-timeout=20 \
          --timeout=15 --waitretry=0 --show-progress --progress=bar:force \
          >${OUT1:?} 2>${OUT2:?} || { error "${EC24:?}"; return 1; }
      fi

      if [ "${2:?}" == "tar.gz" ]; then
        tar zxvf "${FILE}" >${OUT1:?} 2>${OUT2:?} || { error "${EC25:?} (gz)"; return 1; }
      elif [ "${2:?}" == "tar.xz" ]; then
        tar xf "${FILE}" >${OUT1:?} 2>${OUT2:?} || { error "${EC25:?} (xz)"; return 1; }
      else
        error "UNKNOWN FILE EXTENSION"; return 1;
      fi

      if [[ "${1:?}/" != "${4:?}" && "${1:?}" != "${4:?}" && "${1:?}" != "${4:?}/" ]]; then
        # In case this script runs twice (after being killed by CTRL-D)
        rm -rf "${CCIL:?}/${4:?}"
        mv "${1:?}/" "${4:?}" 2>${OUT2:?} || { error "MV FOLDER"; return 1; }
      fi
    fi

    cdfolder "${ROOTDIR}" || return 1;
  }

  wgetact2() {
    #ARGUMENTS: FOLDER, FILE, URL, NEW_FOLDER, XZFILE
    
    local CCIL="${ROOTDIR}/../cocoa_installation_libraries"

    cdfolder "${CCIL:?}" || return 1;

    # In case this script runs twice -------------------------------------------
    if [ -n "${OVERWRITE_EXISTING_CORE_PACKAGES}" ]; then
    
      rm -f  "${CCIL:?}/${5:?}"
    
    fi

    if [ ! -e "${CCIL:?}/${5:?}" ]; then
    
      # Why this compress? In the old Cocoa we saved the libraries in 
      # the github repo using git lfs. So this way preserves the old scripts
      # we set "-k 1" to be the minimum compression
      tar -cf - "${4:?}" | xz -k -1 --threads=$MNT -c - > "${5}" \
        2>${OUT2:?} || { error "TAR (COMPRESS)"; return 1; }

    fi

    cdfolder "${ROOTDIR}" || return 1;
  }

  wgetact() {
    #ARGUMENTS: FOLDER, FILE, URL, NEW_FOLDER, XZFILE
    
    wgetact1 "${1:?}" "${2:?}" "${3:?}" "${4:?}" "${5:?}"

    wgetact2 "${1:?}" "${2:?}" "${3:?}" "${4:?}" "${5:?}"

  }

  gitact1() {
    #ARGUMENTS: FOLDER, VERSION, URL
    
    local CCIL="${ROOTDIR}/../cocoa_installation_libraries"

    cdfolder "${CCIL:?}" || return 1;

    # In case this script runs twice -------------------------------------------
    if [ -n "${OVERWRITE_EXISTING_CORE_PACKAGES}" ]; then
    
      rm -rf "${CCIL:?}/${1:?}"
    
    fi

    if [ ! -d "${CCIL:?}/${1:?}" ]; then
    
      "${CURL:?}" -fsS "${3:?}" \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC27:?} (URL=${3:?})"; return 1; }

      "${GIT:?}" clone "${3:?}" "${1:?}" \
        >${OUT1:?} 2>${OUT2:?} || { error "GIT CLONE"; return 1; }
      
      cdfolder "${CCIL:?}/${1:?}" || return 1;
      
      "${GIT:?}" checkout "${2:?}" \
        >${OUT1:?} 2>${OUT2:?} || { error "GIT CHECKOUT"; return 1; }
      
      rm -rf "${CCIL:?}/${1:?}/.git/"
    
    fi

    cdfolder "${ROOTDIR}" || return 1;
  }

  gitact2() {
    #ARGUMENTS: FOLDER, NEW_FOLDER, XZFILE
    
    local CCIL="${ROOTDIR}/../cocoa_installation_libraries"

    # In case this script runs twice -------------------------------------------
    if [ -n "${OVERWRITE_EXISTING_CORE_PACKAGES}" ]; then

      rm -f  "${CCIL:?}/${3:?}"

    fi

    if [ ! -d "${CCIL:?}/${3:?}" ]; then
    
      cdfolder "${CCIL:?}" || return 1;

      if [ "${1:?}/" != "${2:?}" ] && [ "${1:?}" != "${2:?}" ] \
        && [ "${1:?}" != "${2:?}/" ]; then

        # In case this script runs twice (after being killed by CTRL-D)
        rm -rf "${CCIL:?}/${2:?}"

        mv "${1:?}/" "${2:?}" 2>${OUT2:?} || { error "MV FOLDER"; return 1; }
      
      fi

      # Why this compress? In the old Cocoa we saved the libraries in 
      # the github repo using git lfs. So this way preserves the old scripts
      # we set "-k 1" to be the minimum compression
      tar -cf - "${2:?}" | xz -k -1 --threads=$MNT -c - > "${3}" \
        2>${OUT2:?} || { error "TAR (COMPRESS)"; return 1; }

    fi

    cdfolder "${ROOTDIR}" || return 1;
  }

  gitact() {
    #ARGUMENTS: FOLDER, VERSION, URL, NEW_FOLDER, XZFILE
    gitact1 "${1:?}" "${2:?}" "${3:?}" || return 1;
    
    gitact2 "${1:?}" "${4:?}" "${5:?}" || return 1;
  }

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
   
  unset_env_vars || return 1;

  CCIL="${ROOTDIR}/../cocoa_installation_libraries"

  # ----------------------------------------------------------------------------
  # ------------------------------ XZ LIBRARY ----------------------------------
  # ----------------------------------------------------------------------------

  if [ -z "${IGNORE_XZ_INSTALLATION}" ]; then
  
    ptop "GETTING AND COMPILING XZ LIBRARY (CORE LIBS)" || return 1;

    cdfolder "${CCIL:?}" || return 1;

    #False xz file: just to trigger GIT LFS
    cp xz-5.2.5.tar.gz.xz xz-5.2.5.tar.gz \
      2>${OUT2} ||  { error "CP XZ TAR"; return 1; }

    tar -xf xz-5.2.5.tar.gz.xz \
      >${OUT1:?} 2>${OUT2:?} ||  { error "TAR XZ TAR"; return 1; }

    cdfolder "${CCIL:?}/xz-5.2.5/" || return 1;

    CC="${C_COMPILER:?}" ./configure --prefix="${ROOTDIR:?}/.local" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC11:?}"; return 1; }

    make -j $MNT all >${OUT1:?} 2>${OUT2:?} || { error "${:EC8?}"; return 1; }

    make install >${OUT1:?} 2>${OUT2:?} || { error "${EC10:?}"; return 1; }

    cdfolder "${ROOTDIR}" || return 1;

    pbottom "GETTING AND COMPILING XZ LIBRARY (CORE LIBS)" || return 1;
  
  fi

  # ----------------------------------------------------------------------------
  # ------------------------------ CMAKE LIBRARY --------------------------------
  # ----------------------------------------------------------------------------

  if [ -z "${IGNORE_CMAKE_INSTALLATION}" ]; then
    
    ptop "GETTING CMAKE LIBRARY (CORE LIBS)" || return 1;

    URL='https://github.com/Kitware/CMake.git'
    
    PACKAGE_VERSION="${CMAKE_VERSION:-"3.26.4"}"

    FOLDER="cmake-${PACKAGE_VERSION:?}"
    
    VER=v${PACKAGE_VERSION:?}

    # note: Do not change XFZ filename. 
    # note: Otherwise the script unxv_core_packages.sh will stop working
    XZF="cmake.xz"

    gitact "${FOLDER:?}" "${VER:?}" "${URL:?}" \
      "${COCOA_CMAKE_DIR:-"${FOLDER:?}"}" "${XZF:?}" || return 1;

    unset -v URL FOLDER VER XZF PACKAGE_VERSION

    # COMPILING CMAKE

    pbottom "GETTING CMAKE LIBRARY (CORE LIBS)" || return 1;

  fi

  # ----------------------------------------------------------------------------
  # ------------------------------ WGET LIBRARY --------------------------------
  # ----------------------------------------------------------------------------

  if [ -z "${IGNORE_WGET_INSTALLATION}" ]; then
    
    ptop "GETTING WGET LIBRARY (CORE LIBS)" || return 1;

    PACKAGE_VERSION="${WGET_VERSION:-"1.24.5"}"

    FOLDER="wget-${PACKAGE_VERSION:?}"

    FILE="tar.gz"

    URL="https://ftp.gnu.org/gnu/wget/${FOLDER:?}.${FILE}"

    # note: Do not change XFZ filename. 
    # note: Otherwise the script unxv_core_packages.sh will stop working
    XZF="wget.xz"

    # note: Circular dependence (wget depends on wget!). Why? 
    # note: The wget commands show progress bar when downloading large likelihoods
    # note: The specific progress bar option only exists on recent wget versions  
    wgetact "${FOLDER:?}" "${FILE:?}" "${URL:?}" \
      "${COCOA_WGET_DIR:-"${FOLDER:?}"}" "${XZF:?}" || return 1;

    unset -v URL FOLDER FILE XZF PACKAGE_VERSION

    pbottom "GETTING WGET LIBRARY (CORE LIBS)" || return 1;

  fi

  # ----------------------------------------------------------------------------
  # ------------------------------ BinUtils LIBRARY ----------------------------
  # ----------------------------------------------------------------------------

  if [ -z "${IGNORE_DISTUTILS_INSTALLATION}" ]; then

    ptop "GETTING BINUTILS LIBRARY (CORE LIBS)" || return 1;

    PACKAGE_VERSION="${BINUTILS_VERSION:-"2.37"}"

    FOLDER="binutils-${PACKAGE_VERSION:?}"
    
    FILE="tar.gz"
    
    URL="https://ftp.gnu.org/gnu/binutils/${FOLDER:?}.${FILE:?}"
    
    # note: Do not change XFZ filename. 
    # note: Otherwise the script unxv_core_packages.sh will stop working
    XZF="binutils.xz"

    wgetact "${FOLDER:?}" "${FILE:?}" "${URL:?}" \
      "${COCOA_BINUTILS_DIR:-"${FOLDER:?}"}" "${XZF:?}" || return 1;

    unset -v URL FOLDER FILE XZF PACKAGE_VERSION

    pbottom "GETTING BINUTILS LIBRARY (CORE LIBS)" || return 1;

  fi

  # ----------------------------------------------------------------------------
  # ------------------------------ TEXTINFO LIBRARY ----------------------------
  # ----------------------------------------------------------------------------

  if [ -z "${IGNORE_DISTUTILS_INSTALLATION}" ]; then
    
    ptop "GETTING TEXINFO LIBRARY (CORE LIBS)"

    PACKAGE_VERSION="${TEXTINFO_VERSION:-"7.0.3"}"

    FOLDER="texinfo-${PACKAGE_VERSION:?}"
    
    FILE="tar.xz"
    
    URL="https://ftp.gnu.org/gnu/texinfo/${FOLDER}.${FILE}"
   
    # note: Do not change XFZ filename. 
    # note: Otherwise the script unxv_core_packages.sh will stop working 
    XZF="texinfo.xz"

    wgetact "${FOLDER:?}" "${FILE:?}" "${URL:?}" \
      "${COCOA_TEXINFO_DIR:-"${FOLDER:?}"}" "${XZF:?}" || return 1;

    unset -v URL FOLDER FILE XZF PACKAGE_VERSION

    pbottom "GETTING TEXINFO LIBRARY (CORE LIBS)" || return 1;

  fi

  # ----------------------------------------------------------------------------
  # ------------------------------ OpenBLAS LIBRARY ----------------------------
  # ----------------------------------------------------------------------------

  if [ -z "${IGNORE_OPENBLAS_INSTALLATION}" ]; then

    ptop "GETTING OPENBLAS LIBRARY (CORE LIBS)" || return 1;

    PACKAGE_VERSION="${OPENBLAS_VERSION:-"0.3.23"}" 

    FOLDER="OpenBLAS-${PACKAGE_VERSION:?}"

    URL='https://github.com/OpenMathLib/OpenBLAS.git'

    VER="v${PACKAGE_VERSION:?}"

    # note: Do not change XFZ filename. 
    # note: Otherwise the script unxv_core_packages.sh will stop working
    XZF="OpenBLAS.xz"

    gitact "${FOLDER:?}" "${VER:?}" "${URL:?}" \
      "${COCOA_OPENBLAS_DIR:-"${FOLDER:?}"}" "${XZF:?}" || return 1;

    unset -v URL FOLDER VER XZF PACKAGE_VERSION

    pbottom "GETTING OPENBLAS LIBRARY (CORE LIBS)" || return 1;

  fi

  # ----------------------------------------------------------------------------
  # ------------------------------ LAPACK LIBRARY ------------------------------
  # ----------------------------------------------------------------------------

  if [ -z "${IGNORE_FORTRAN_LAPACK_INSTALLATION}" ]; then

    ptop "GETTING LAPACK LIBRARY (CORE LIBS)" || return 1;

    PACKAGE_VERSION="${LAPACK_VERSION:-"3.11"}" 

    FOLDER="lapack-${PACKAGE_VERSION:?}"

    URL='https://github.com/Reference-LAPACK/lapack.git'

    VER=v"${PACKAGE_VERSION:?}"

    # note: Do not change XFZ filename. 
    # note: Otherwise the script unxv_core_packages.sh will stop working
    XZF="lapack.xz"

    gitact "${FOLDER:?}" "${VER:?}" "${URL:?}" \
      "${COCOA_LAPACK_DIR:-"${FOLDER:?}"}" "${XZF:?}" || return 1;

    unset -v URL FOLDER VER XZF PACKAGE_VERSION

    pbottom "GETTING LAPACK LIBRARY DONE (CORE LIBS)" || return 1;

  fi

  # ----------------------------------------------------------------------------
  # ------------------------------ HDF5 LIBRARY --------------------------------
  # ----------------------------------------------------------------------------

  if [ -z "${IGNORE_HDF5_INSTALLATION}" ]; then

    ptop "GETTING HDF5 LIBRARY (CORE LIBS)" || return 1;

    FOLDER="hdf5-1.12.3"

    FILE="tar.gz"

    URL_BASE="https://hdf-wordpress-1.s3.amazonaws.com/wp-content/uploads"
    
    URL="${URL_BASE:?}/manual/HDF5/HDF5_1_12_3/src/${FOLDER}.${FILE}"

    # note: Do not change XFZ filename. 
    # note: Otherwise the script unxv_core_packages.sh will stop working
    XZF="hdf5.xz"

    wgetact "${FOLDER:?}" "${FILE:?}" "${URL:?}" \
      "${COCOA_HDF5_DIR:-"${FOLDER:?}"}" "${XZF:?}" || return 1;

    unset -v URL_BASE URL FOLDER FILE XZF

    pbottom "GETTING HDF5 LIBRARY DONE (CORE LIBS)" || return 1;

  fi

  # ----------------------------------------------------------------------------
  # ------------------------------ CFITSIO LIBRARY -----------------------------
  # ----------------------------------------------------------------------------

  if [ -z "${IGNORE_C_CFITSIO_INSTALLATION}" ]; then
    
    ptop "GETTING CFITSIO LIBRARY (CORE LIBS)" || return 1;

    PACKAGE_VERSION="${CFITSIO_VERSION:-"4.0.0"}" 

    FOLDER="cfitsio-${PACKAGE_VERSION:?}"
    
    FILE="tar.gz"
    
    URL_BASE="http://heasarc.gsfc.nasa.gov"
    
    URL="${URL_BASE}/FTP/software/fitsio/c/${FOLDER}.${FILE}"
   
    # note: Do not change XFZ filename. 
    # note: Otherwise the script unxv_core_packages.sh will stop working 
    XZF="cfitsio.xz"

    wgetact "${FOLDER:?}" "${FILE:?}" "${URL:?}" \
      "${COCOA_CFITSIO_DIR:-"${FOLDER:?}"}" "${XZF:?}" || return 1;

    unset -v URL_BASE URL FOLDER FILE XZF PACKAGE_VERSION

    pbottom "GETTING CFITSIO LIBRARY DONE (CORE LIBS)" || return 1;

  fi

  # ----------------------------------------------------------------------------
  # ------------------------------ FFTW LIBRARY --------------------------------
  # ----------------------------------------------------------------------------

  if [ -z "${IGNORE_C_FFTW_INSTALLATION}" ]; then

    ptop "GETTING FFTW LIBRARY (CORE LIBS)" || return 1;

    PACKAGE_VERSION="${FFTW_VERSION:-"3.3.10"}" 

    FOLDER="fftw-${PACKAGE_VERSION:?}"

    FILE="tar.gz"

    URL="http://www.fftw.org/${FOLDER}.${FILE}"

    # note: Do not change XFZ filename. 
    # note: Otherwise the script unxv_core_packages.sh will stop working
    XZF="fftw.xz"

    wgetact "${FOLDER:?}" "${FILE:?}" "${URL:?}" \
      "${COCOA_FFTW_DIR:-"${FOLDER:?}"}" "${XZF:?}" || return 1;

    unset -v URL FOLDER FILE XZF PACKAGE_VERSION

    pbottom "GETTING FFTW LIBRARY DONE (CORE LIBS)" || return 1;

  fi

  # ----------------------------------------------------------------------------
  # ------------------------------ GSL LIBRARY ---------------------------------
  # ----------------------------------------------------------------------------

  if [ -z "${IGNORE_C_GSL_INSTALLATION}" ]; then

    ptop "GETTING GSL LIBRARY (CORE LIBS)" || return 1;
    
    PACKAGE_VERSION="${GSL_VERSION:-"2.7"}" 

    FOLDER="gsl-${PACKAGE_VERSION:?}"

    FILE="tar.gz"

    URL="http://ftp.wayne.edu/gnu/gsl/${FOLDER}.${FILE}"

    # note: Do not change XFZ filename. 
    # note: Otherwise the script unxv_core_packages.sh will stop working
    XZF="gsl.xz"

    wgetact "${FOLDER:?}" "${FILE:?}" "${URL:?}" \
      "${COCOA_GSL_DIR:-"${FOLDER:?}"}" "${XZF:?}" || return 1;

    unset -v URL FOLDER FILE XZF PACKAGE_VERSION

    pbottom "GETTING GSL LIBRARY DONE (CORE LIBS)" || return 1;

  fi

  # ----------------------------------------------------------------------------
  # ------------------------------ SPDLOG LIBRARY ------------------------------
  # ----------------------------------------------------------------------------

  if [ -z "${IGNORE_CPP_SPDLOG_INSTALLATION}" ]; then

    ptop "GETTING SPDLOG LIBRARY (CORE LIBS)" || return 1;

    #PACKAGE_VERSION="${SPDLOG_VERSION:-"v1.14.1"}"  
    PACKAGE_VERSION="${SPDLOG_VERSION:-"v1.15.3"}"  

    URL='https://github.com/gabime/spdlog.git'

    FOLDER='spdlog'

    VER="${PACKAGE_VERSION:?}"

    # note: Do not change XFZ filename. 
    # note: Otherwise the script unxv_core_packages.sh will stop working
    XZF="spdlog.xz"

    gitact "${FOLDER:?}" "${VER:?}" "${URL:?}" \
      "${COCOA_SPDLOG_DIR:-"spdlog"}" "${XZF:?}" || return 1;

    unset -v URL FOLDER VER XZF

    pbottom "GETTING SPDLOG LIBRARY (CORE LIBS)" || return 1;

  fi

  # ----------------------------------------------------------------------------
  # ------------------------------ ARMA LIBRARY --------------------------------
  # ----------------------------------------------------------------------------

  if [ -z "${IGNORE_CPP_ARMA_INSTALLATION}" ]; then
    
    ptop "GETTING ARMA LIBRARY DONE (CORE LIBS)" || return 1;

    PACKAGE_VERSION="${ARMA_VERSION:-"12.8.4"}" 
    FOLDER="armadillo-${PACKAGE_VERSION:?}"
    FILE="tar.xz"
    URL="https://sourceforge.net/projects/arma/files/${FOLDER}.${FILE}"

    #PACKAGE_VERSION="${ARMA_VERSION:-"12.8.4"}" 
    #FOLDER="armadillo-code-${PACKAGE_VERSION}"
    #FILE="tar.gz"
    #URL="https://gitlab.com/conradsnicta/armadillo-code/-/archive/${PACKAGE_VERSION}/${FOLDER}.${FILE}"
    
    # note: Do not change XFZ filename. 
    # note: Otherwise the script unxv_core_packages.sh will stop working
    XZF="armadillo.xz"

    wgetact "${FOLDER:?}" "${FILE:?}" "${URL:?}" \
      "${COCOA_ARMADILLO_DIR:-"${FOLDER:?}"}" "${XZF:?}" || return 1;

    unset -v URL FOLDER FILE XZF

    pbottom "GETTING ARMA LIBRARY DONE (CORE LIBS)" || return 1;

  fi

  # ----------------------------------------------------------------------------
  # ------------------------------ BOOST LIBRARY -------------------------------
  # ----------------------------------------------------------------------------

  if [ -z "${IGNORE_CPP_BOOST_INSTALLATION}" ]; then

    ptop "GETTING BOOST LIBRARY (CORE LIBS)" || return 1;

    PACKAGE_VERSION="${CPP_BOOST_VERSION:-"81"}"

    FOLDER="boost_1_${PACKAGE_VERSION:?}_0"

    FILE="tar.gz"

    URL_BASE="https://boostorg.jfrog.io/artifactory/main/release/"
    
    URL="${URL_BASE}/1.${PACKAGE_VERSION:?}.0/source/${FOLDER}.${FILE}"

    # note: Do not change XFZ filename. 
    # note: Otherwise the script unxv_core_packages.sh will stop working
    XZF="boost.xz"

    wgetact "${FOLDER:?}" "${FILE:?}" "${URL:?}" \
      "${COCOA_BOOST_DIR:-"${FOLDER:?}"}" "${XZF:?}" || return 1;

    unset -v URL_BASE URL FOLDER FILE XZF PACKAGE_VERSION

    pbottom "GETTING BOOST LIBRARY (CORE LIBS)" || return 1;

  fi

  # ----------------------------------------------------------------------------
  # ------------------------------ CUBA LIBRARY --------------------------------
  # ----------------------------------------------------------------------------

  if [ -z "${IGNORE_CPP_CUBA_INSTALLATION}" ]; then

    ptop "GETTING CUBA LIBRARY (CORE LIBS)" || return 1;

    PACKAGE_VERSION="${CPP_CUBA_VERSION:-"4.2.2"}"

    FOLDER="Cuba-${PACKAGE_VERSION:?}"

    FILE="tar.gz"
    
    URL="https://feynarts.de/cuba/${FOLDER}.${FILE}"

    # note: Do not change XFZ filename. 
    # note: Otherwise the script unxv_core_packages.sh will stop working
    XZF="cuba.xz"

    wgetact "${FOLDER:?}" "${FILE:?}" "${URL:?}" \
      "${COCOA_CUBA_DIR:-"${FOLDER:?}"}" "${XZF:?}" || return 1;

    unset -v URL FOLDER FILE XZF PACKAGE_VERSION PACKDIR

    pbottom "GETTING CUBA LIBRARY (CORE LIBS)" || return 1;

  fi

  # ----------------------------------------------------------------------------
  # ------------------------------ CARMA LIBRARY -------------------------------
  # ----------------------------------------------------------------------------

  if [ -z "${IGNORE_CPP_CARMA_INSTALLATION}" ]; then

    ptop "GETTING CARMA LIBRARY DONE (CORE LIBS)" || return 1;

    CNAME="${COCOA_CARMA_DIR:-"carma"}"
    
    PACKDIR="${CCIL:?}/${CNAME:?}"

    FOLDER="carma_tmp"

    URL="https://github.com/RUrlus/carma.git"

    VER=v0.7.0

    XZF="carma.xz"

    # In case this script runs twice -------------------------------------------
    if [ -n "${OVERWRITE_EXISTING_CORE_PACKAGES}" ]; then
      
      rm -rf "${CCIL:?}/${FOLDER:?}"
      rm -rf "${PACKDIR:?}"
      rm -rf "${CCIL:?}/include"
    
    fi

    if [ ! -d "${CCIL:?}/${FOLDER:?}" ]; then
      
      gitact1 "${FOLDER:?}" "${VER:?}" "${URL:?}" || return 1;

      # --------------------------------------------------------------------------
      # move/rename include file and carma.h folder
      # --------------------------------------------------------------------------
      mv "${CCIL:?}/${FOLDER:?}/include" "${CCIL:?}" >${OUT1:?} 2>${OUT2:?} || 
        { error "MV CARMA INCLUDE FOLDER"; return 1; }

      mv "${CCIL:?}/include" "${PACKDIR:?}" 2>${OUT2:?} || 
        { error "RENAME CARMA INCLUDE FOLDER"; return 1; }

      mv "${PACKDIR:?}/carma" "${PACKDIR:?}/carma.h" \
        2>${OUT2:?} || { error "RENANE CARMA HEADER"; return 1; }

      # --------------------------------------------------------------------------
      
      gitact2 "${CNAME:?}" "${CNAME:?}" "${XZF:?}" || return 1; 

      unset -v CNAME PACKDIR URL FOLDER VER XZF

    fi

    pbottom "GETTING CARMA LIBRARY DONE (CORE LIBS)" || return 1;

  fi

  # ----------------------------------------------------------------------------

  unset_all || return 1;

fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------