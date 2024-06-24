#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${ROOTDIR}" ]; then
  pfail 'ROOTDIR'; return 1
fi

source "${ROOTDIR:?}/installation_scripts/.check_flags.sh" || return 1;

unset_env_vars_sil () {
  unset -v URL_BASE URL FOLDER VER XZF CCIL CNAME PACKDIR
  cdroot || return 1;
}

unset_env_funcs () {
  unset -f cdfolder cpfolder error wgetact gitact gitact1 gitact2
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
  fail_script_msg "setup_installation_libraries.sh" "${1}"
  unset_all || return 1
}

wgetact() {
  #ARGUMENTS: FOLDER, FILE, URL, NEW_FOLDER, XZFILE
  
  local CCIL="${ROOTDIR}/../cocoa_installation_libraries"

  cdfolder "${CCIL:?}" || return 1;

  local FILE="${1}.${2}"

  # In case this script runs twice (after being killed by CTRL-D)
  rm -rf "${CCIL:?}/${1:?}"
  rm -f  "${CCIL:?}/${FILE:?}"
  
  rm -f  "${CCIL:?}/${5:?}"

  wget -q "${3:?}" >${OUT1:?} 2>${OUT2:?} || { error "${EC24:?}"; return 1; }

  if [ "${2:?}" == "tar.gz" ]; then
    
    tar zxvf "${FILE}" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC25:?} (gz)"; return 1; }
  
  elif [ "${2:?}" == "tar.xz" ]; then
  
    tar xf "${FILE}" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC25:?} (xz)"; return 1; }
  
  else

    error "UNKNOWN FILE EXTENSION"; return 1;
  
  fi

  if [ "${1:?}/" != "${4:?}" ] && [ "${1:?}" != "${4:?}" ] \
    && [ "${1:?}" != "${4:?}/" ] ; then

    # In case this script runs twice (after being killed by CTRL-D)
    rm -rf "${CCIL:?}/${4:?}"
    
    mv "${1:?}/" "${4:?}" 2>${OUT2:?} || { error "MV FOLDER"; return 1; }
  
  fi

  # Why this compress? In the old Cocoa we saved the libraries in 
  # the github repo using git lfs. So this way preserves the old scripts
  # we set "-k 1" to be the minimum compression
  tar -cf - "${4:?}" | xz -k -1 --threads=$MNT -c - > "${5}" \
    2>${OUT2:?} || { error "TAR (COMPRESS)"; return 1; }

  rm -rf "${CCIL:?}/${1:?}"
  rm -f  "${CCIL:?}/${FILE:?}"
  rm -rf "${CCIL:?}/${4:?}"

  cdfolder "${ROOTDIR}" || return 1;
}

gitact1() {
  #ARGUMENTS: FOLDER, VERSION, URL
  
  local CCIL="${ROOTDIR}/../cocoa_installation_libraries"

  cdfolder "${CCIL:?}" || return 1;

  # ---------------------------------------------------------------------------
  # In case this script runs twice
  rm -rf "${CCIL:?}/${1:?}"
  # ---------------------------------------------------------------------------

  "${GIT:?}" clone "${3:?}" "${1:?}" >${OUT1:?} 2>${OUT2:?} || 
    { error "GIT CLONE"; return 1; }
  
  cdfolder "${CCIL:?}/${1:?}" || return 1;
  
  "${GIT:?}" checkout "${2:?}" >${OUT1:?} 2>${OUT2:?} || 
    { error "GIT CHECKOUT"; return 1; }
  
  rm -rf "${CCIL:?}/${1:?}/.git/"

  cdfolder "${ROOTDIR}" || return 1;
}

gitact2() {
  #ARGUMENTS: FOLDER, NEW_FOLDER, XZFILE
  
  local CCIL="${ROOTDIR}/../cocoa_installation_libraries"

  # In case this script runs twice (after being killed by CTRL-D)
  rm -f  "${CCIL:?}/${3:?}"

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

  rm -rf "${CCIL:?}/${1:?}"
  rm -rf "${CCIL:?}/${2:?}"

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
 
ptop2 'SETUP_CORE_LIBRARIES'

unset_env_vars || return 1;

CCIL="${ROOTDIR}/../cocoa_installation_libraries"

# ----------------------------------------------------------------------------
# --------------------------- CMAKE LIBRARY ----------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_CMAKE_INSTALLATION}" ]; then
  
  ptop "GETTING CMAKE LIBRARY" || return 1;

  URL='https://github.com/Kitware/CMake.git'
  
  FOLDER='cmake-3.26.4'
  
  VER=v3.26.4
  
  XZF="cmake.xz"

  gitact "${FOLDER}" "${VER}" "${URL}" \
    "${COCOA_CMAKE_DIR:-"cmake-3.26.4"}" "${XZF}" || return 1;

  unset -v URL FOLDER VER XZF

  pbottom "GETTING CMAKE LIBRARY" || return 1;

fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_DISTUTILS_INSTALLATION}" ]; then

  ptop "GETTING BINUTILS LIBRARY" || return 1;

  FOLDER="binutils-2.37"
  
  FILE="tar.gz"
  
  URL="https://ftp.gnu.org/gnu/binutils/${FOLDER}.${FILE}"
  
  XZF="binutils.xz"

  wgetact "${FOLDER}" "${FILE}" "${URL}" \
    "${COCOA_BINUTILS_DIR:-"binutils-2.37"}" "${XZF}" || return 1;

  unset -v URL FOLDER FILE XZF

  pbottom "GETTING BINUTILS LIBRARY" || return 1;

fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_DISTUTILS_INSTALLATION}" ]; then || return 1;
  
  ptop "GETTING TEXINFO LIBRARY"

  FOLDER="texinfo-7.0.3"
  
  FILE="tar.xz"
  
  URL="https://ftp.gnu.org/gnu/texinfo/${FOLDER}.${FILE}"
  
  XZF="texinfo.xz"

  wgetact "${FOLDER}" "${FILE}" "${URL}" \
    "${COCOA_TEXINFO_DIR:-"texinfo-7.0.3"}" "${XZF}" || return 1;

  unset -v URL FOLDER FILE XZF

  pbottom "GETTING TEXINFO LIBRARY" || return 1;

fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_OPENBLAS_INSTALLATION}" ]; then

  ptop "GETTING OPENBLAS LIBRARY" || return 1;

  URL='https://github.com/OpenMathLib/OpenBLAS.git'

  FOLDER='OpenBLAS-0.3.23'

  VER=v0.3.23

  XZF="OpenBLAS.xz"

  gitact "${FOLDER}" "${VER}" "${URL}" \
    "${COCOA_OPENBLAS_DIR:-"OpenBLAS-0.3.23"}" "${XZF}" || return 1;

  unset -v URL FOLDER VER XZF

  pbottom "GETTING OPENBLAS LIBRARY" || return 1;

fi

# ----------------------------------------------------------------------------
# ------------------------------- LPACK LIBRARY ------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_FORTRAN_LAPACK_INSTALLATION}" ]; then

  ptop "GETTING LAPACK LIBRARY" || return 1;

  URL='https://github.com/Reference-LAPACK/lapack.git'

  FOLDER='lapack-3.11.0'

  VER=v3.11

  XZF="lapack.xz"

  gitact "${FOLDER}" "${VER}" "${URL}" \
    "${COCOA_LAPACK_DIR:-"lapack-3.11.0"}" "${XZF}" || return 1;

  unset -v URL FOLDER VER XZF

  pbottom "GETTING LAPACK LIBRARY DONE" || return 1;

fi

# ----------------------------------------------------------------------------
# ---------------------------- HDF5 LIBRARY ----------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_HDF5_INSTALLATION}" ]; then

  ptop "GETTING HDF5 LIBRARY" || return 1;

  FOLDER="hdf5-1.12.3"

  FILE="tar.gz"

  URL_BASE="https://hdf-wordpress-1.s3.amazonaws.com/wp-content/uploads"
  URL="${URL_BASE:?}/manual/HDF5/HDF5_1_12_3/src/${FOLDER}.${FILE}"

  XZF="hdf5.xz"

  wgetact "${FOLDER}" "${FILE}" "${URL}" \
    "${COCOA_HDF5_DIR:-"hdf5-1.12.3"}" "${XZF}" || return 1;

  unset -v URL_BASE URL FOLDER FILE XZF

  pbottom "GETTING HDF5 LIBRARY DONE" || return 1;

fi

# ----------------------------------------------------------------------------
# ------------------------------ CFITSIO LIBRARY -----------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_C_CFITSIO_INSTALLATION}" ]; then
  
  ptop "GETTING CFITSIO LIBRARY" || return 1;

  FOLDER="cfitsio-4.0.0"
  
  FILE="tar.gz"
  
  URL_BASE="http://heasarc.gsfc.nasa.gov"
  URL="${URL_BASE}/FTP/software/fitsio/c/${FOLDER}.${FILE}"
  
  XZF="cfitsio.xz"

  wgetact "${FOLDER}" "${FILE}" "${URL}" \ 
    "${COCOA_CFITSIO_DIR:-"cfitsio-4.0.0"}" "${XZF}" || return 1;

  unset -v URL_BASE URL FOLDER FILE XZF

  pbottom "GETTING CFITSIO LIBRARY DONE" || return 1;

fi

# ----------------------------------------------------------------------------
# ---------------------------- FFTW LIBRARY ----------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_C_FFTW_INSTALLATION}" ]; then

  ptop "GETTING FFTW LIBRARY" || return 1;

  FOLDER="fftw-3.3.10"

  FILE="tar.gz"

  URL="http://www.fftw.org/${FOLDER}.${FILE}"

  XZF="fftw.xz"

  wgetact "${FOLDER}" "${FILE}" "${URL}" \
    "${COCOA_FFTW_DIR:-"fftw-3.3.10"}" "${XZF}" || return 1;

  unset -v URL FOLDER FILE XZF

  pbottom "GETTING FFTW LIBRARY DONE" || return 1;

fi

# ----------------------------------------------------------------------------
# ------------------------------ GSL LIBRARY ---------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_C_GSL_INSTALLATION}" ]; then

  ptop "GETTING GSL LIBRARY" || return 1;
  
  FOLDER="gsl-2.7"

  FILE="tar.gz"

  URL="http://ftp.wayne.edu/gnu/gsl/${FOLDER}.${FILE}"

  XZF="gsl-2.7.xz"

  wgetact "${FOLDER}" "${FILE}" "${URL}" \ 
    "${COCOA_GSL_DIR:-"gsl-2.7"}" "${XZF}" || return 1;

  unset -v URL FOLDER FILE XZF

  pbottom "GETTING GSL LIBRARY DONE" || return 1;

fi

# ----------------------------------------------------------------------------
# ---------------------------- SPDLOG LIBRARY --------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_CPP_SPDLOG_INSTALLATION}" ]; then

  ptop "GETTING SPDLOG LIBRARY" || return 1;

  URL='https://github.com/gabime/spdlog.git'

  FOLDER='spdlog'

  VER=v1.13.0

  XZF="spdlog.xz"

  gitact "${FOLDER}" "${VER}" "${URL}" \
    "${COCOA_SPDLOG_DIR:-"spdlog"}" "${XZF}" || return 1;

  unset -v URL FOLDER VER XZF

  pbottom "GETTING SPDLOG LIBRARY" || return 1;

fi

# ----------------------------------------------------------------------------
# -------------------------- ARMA LIBRARY ------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_CPP_ARMA_INSTALLATION}" ]; then
  
  ptop "GETTING ARMA LIBRARY DONE" || return 1;

  FOLDER="armadillo-12.8.2"
  
  FILE="tar.xz"
  
  URL="https://sourceforge.net/projects/arma/files/${FOLDER}.${FILE}"
  
  XZF="armadillo.xz"

  wgetact "${FOLDER}" "${FILE}" "${URL}" \
    "${COCOA_ARMADILLO_DIR:-"armadillo-12.8.2"}" "${XZF}" || return 1;

  unset -v URL FOLDER FILE XZF

  pbottom "GETTING ARMA LIBRARY DONE" || return 1;

fi

# ----------------------------------------------------------------------------
# -------------------------- BOOST LIBRARY -----------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_CPP_BOOST_INSTALLATION}" ]; then

  ptop "GETTING BOOST LIBRARY" || return 1;

  if [ -z "${COCOA_BOOST_DIR}" ]; then
    pfail 'COCOA_BOOST_DIR'; cdroot; return 1;
  fi

  FOLDER="boost_1_81_0"

  FILE="tar.gz"

  URL_BASE="https://boostorg.jfrog.io/artifactory/main/release/"
  URL="${URL_BASE}/1.81.0/source/${FOLDER}.${FILE}"

  XZF="boost.xz"

  wgetact "${FOLDER}" "${FILE}" "${URL}" \
    "${COCOA_BOOST_DIR:-"boost_1_81_0"}" "${XZF}" || return 1;

  unset -v URL_BASE URL FOLDER FILE XZF

  pbottom "GETTING BOOST LIBRARY" || return 1;

fi

# ----------------------------------------------------------------------------
# ------------------------ CARMA LIBRARY -------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_CPP_CARMA_INSTALLATION}" ]; then

  ptop "GETTING CARMA LIBRARY DONE" || return 1;

  CNAME="${COCOA_CARMA_DIR:-"carma":?}"
  
  PACKDIR="${CCIL:?}/${CNAME:?}"

  FOLDER="carma_tmp"

  URL="https://github.com/RUrlus/carma.git"

  VER=v0.7.0

  XZF="carma.xz"

  gitact1 "${FOLDER}" "${VER}" "${URL}" || return 1;

  # -------------------------------------------------------------------------
  # In case this script runs twice (after sudden break w/ CTRL-C)
  rm -rf "${PACKDIR:?}"
  rm -rf "${CCIL:?}/include"

  # -------------------------------------------------------------------------
  # move/rename include file and carma.h folder
  # -------------------------------------------------------------------------
  mv "${CCIL:?}/${FOLDER:?}/include" "${CCIL:?}" >${OUT1:?} 2>${OUT2:?} || 
    { error "MV CARMA INCLUDE FOLDER"; return 1; }

  mv "${CCIL:?}/include" "${PACKDIR:?}" 2>${OUT2:?} || 
    { error "RENAME CARMA INCLUDE FOLDER"; return 1; }

  mv "${PACKDIR:?}/carma" "${PACKDIR:?}/carma.h" \
    2>${OUT2:?} || { error "RENANE CARMA HEADER"; return 1; }

  rm -rf "${CCIL:?}/${FOLDER:?}"
  # -------------------------------------------------------------------------
  
  gitact2 "${CNAME}" "${CNAME}" "${XZF}" || return 1; 

  unset -v CNAME PACKDIR URL FOLDER VER XZF

  pbottom "GETTING CARMA LIBRARY DONE" || return 1;

fi

# ----------------------------------------------------------------------------

unset_all || return 1;

pbottom2 "SETUP_CORE_LIBRARIES" || return 1;

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------