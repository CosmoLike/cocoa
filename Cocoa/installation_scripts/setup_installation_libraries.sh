#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
pfail() {
  echo -e \
  "\033[0;31m\t\t ERROR ENV VARIABLE ${1:-"empty arg"} NOT DEFINED \033[0m"
  unset pfail
}

if [ -z "${ROOTDIR}" ]; then
  pfail 'ROOTDIR'; return 1
fi

cdroot() {
  cd "${ROOTDIR:?}" 2>"/dev/null" || { echo -e \
    "\033[0;31m\t\t CD ROOTDIR (${ROOTDIR}) FAILED \033[0m"; return 1; }
  unset cdroot
}

if [ -z "${GIT}" ]; then
  pfail 'GIT'; cdroot; return 1;
fi

unset_env_vars_sil () {
  unset OUT1
  unset OUT2
  unset pfail
  unset URL_BASE
  unset URL
  unset FOLDER
  unset VER
  unset XZF
  unset wgetact
  unset gitact
  unset gitact1
  unset gitact2
  unset CCIL
  unset CNAME 
  unset PACKDIR
  unset unset_env_vars_sil
  cdroot || return 1;
}

fail_sil () {
  local MSG="\033[0;31m\t\t (setup_installation_libraries.sh) WE CANNOT RUN \e[3m"
  local MSG2="\033[0m"
  echo -e "${MSG}${1:-"empty arg"}${MSG2}"
  unset fail_sil
  unset_env_vars_sil
}

if [ -z "${COSMOLIKE_DEBUG_MODE}" ]; then
  export OUT1="/dev/null"; export OUT2="/dev/null"
  export SILMNT="${MAKE_NUM_THREADS:-1}"
else
  export OUT1="/dev/tty"; export OUT2="/dev/tty"
  export SILMNT=1
fi

cdfolder() {
  cd "${1:?}" 2>"/dev/null" || 
    { fail_sil "CD FOLDER: ${1:-"empty arg"}"; return 1; }
}

export CCIL="${ROOTDIR}/../cocoa_installation_libraries"

wgetact() {
  #ARGUMENTS: FOLDER, FILE, URL, NEW_FOLDER, XZFILE
  
  cdfolder "${CCIL:?}" || return 1;

  local FILE="${1}.${2}"

  # In case this script runs twice (after being killed by CTRL-D)
  rm -rf "${CCIL:?}/${1:?}"
  rm -f  "${CCIL:?}/${FILE:?}"
  
  rm -f  "${CCIL:?}/${5:?}"

  wget -q "${3:?}" >${OUT1:?} 2>${OUT2:?} || { fail_sil "WGET"; return 1; }

  if [ "${2:?}" == "tar.gz" ]; then
    
    tar zxvf "${FILE}"  >${OUT1:?} 2>${OUT2:?} || { fail_sil "TAR (gz)"; return 1; }
  
  elif [ "${2:?}" == "tar.xz" ]; then
  
    tar xf "${FILE}" >${OUT1:?} 2>${OUT2:?} || { fail_sil "TAR (xz)"; return 1; }
  
  else

    fail_sil "UNKNOWN FILE EXTENSION"; return 1;
  
  fi

  if [ "${1:?}/" != "${4:?}" ] && [ "${1:?}" != "${4:?}" ] \
    && [ "${1:?}" != "${4:?}/" ] ; then

    # In case this script runs twice (after being killed by CTRL-D)
    rm -rf "${CCIL:?}/${4:?}"
    
    mv "${1:?}/" "${4:?}" 2>${OUT2:?} || { fail_sil "MV FOLDER"; return 1; }
  
  fi

  # Why this compress? In the old Cocoa we saved the libraries in 
  # the github repo using git lfs. So this way preserves the old scripts
  # we set "-k 1" to be the minimum compression
  tar -cf - "${4:?}" | xz -k -1 --threads=$SILMNT -c - > "${5}" \
    2>${OUT2:?} || { fail_sil "TAR (COMPRESS)"; return 1; }

  rm -rf "${CCIL:?}/${1:?}"
  rm -f  "${CCIL:?}/${FILE:?}"
  rm -rf "${CCIL:?}/${4:?}"

  cdfolder "${ROOTDIR}" || return 1;
}

gitact1() {
  #ARGUMENTS: FOLDER, VERSION, URL
  
  cdfolder "${CCIL:?}" || return 1;

  # ---------------------------------------------------------------------------
  # In case this script runs twice
  rm -rf "${CCIL:?}/${1:?}"
  # ---------------------------------------------------------------------------

  "${GIT:?}" clone "${3:?}" "${1:?}" >${OUT1:?} 2>${OUT2:?} || 
    { fail_sil "GIT CLONE"; return 1; }
  
  cdfolder "${CCIL:?}/${1:?}" || return 1;
  
  "${GIT:?}" checkout "${2:?}" >${OUT1:?} 2>${OUT2:?} || 
    { fail_sil "GIT CHECKOUT"; return 1; }
  
  rm -rf "${CCIL:?}/${1:?}/.git/"

  cdfolder "${ROOTDIR}" || return 1;
}

gitact2() {
  #ARGUMENTS: FOLDER, NEW_FOLDER, XZFILE
  
  # In case this script runs twice (after being killed by CTRL-D)
  rm -f  "${CCIL:?}/${3:?}"

  cdfolder "${CCIL:?}" || return 1;

  if [ "${1:?}/" != "${2:?}" ] && [ "${1:?}" != "${2:?}" ] \
    && [ "${1:?}" != "${2:?}/" ]; then

    # In case this script runs twice (after being killed by CTRL-D)
    rm -rf "${CCIL:?}/${2:?}"

    mv "${1:?}/" "${2:?}" 2>${OUT2:?} || { fail_sil "MV FOLDER"; return 1; }
  
  fi

  # Why this compress? In the old Cocoa we saved the libraries in 
  # the github repo using git lfs. So this way preserves the old scripts
  # we set "-k 1" to be the minimum compression
  tar -cf - "${2:?}" | xz -k -1 --threads=$SILMNT -c - > "${3}" \
    2>${OUT2:?} || { fail_sil "TAR (COMPRESS)"; return 1; }

  rm -rf "${CCIL:?}/${1:?}"
  rm -rf "${CCIL:?}/${2:?}"

  cdfolder "${ROOTDIR}" || return 1;
}

gitact() {
  #ARGUMENTS: FOLDER, VERSION, URL, NEW_FOLDER, XZFILE
  gitact1 "${1:?}" "${2:?}" "${3:?}" || return 1;
  
  gitact2 "${1:?}" "${4:?}" "${5:?}" || return 1;
}
 
ptop2 'SETUP_INSTALLATION_LIBRARIES'

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_CMAKE_INSTALLATION}" ]; then
  ptop "GETTING CMAKE LIBRARY"

  export URL='https://github.com/Kitware/CMake.git'
  export FOLDER='cmake-3.26.4'
  export VER=v3.26.4
  export XZF="cmake.xz"

  gitact "${FOLDER}" "${VER}" "${URL}" \
    "${COCOA_CMAKE_DIR:-"cmake-3.26.4"}" "${XZF}" || return 1;

  pbottom "GETTING CMAKE LIBRARY"
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_DISTUTILS_INSTALLATION}" ]; then

  ptop "GETTING BINUTILS LIBRARY"

  export FOLDER="binutils-2.37"
  
  export FILE="tar.gz"
  
  export URL="https://ftp.gnu.org/gnu/binutils/${FOLDER}.${FILE}"
  
  export XZF="binutils.xz"

  wgetact "${FOLDER}" "${FILE}" "${URL}" \
    "${COCOA_BINUTILS_DIR:-"binutils-2.37"}" "${XZF}" || return 1;

  pbottom "GETTING BINUTILS LIBRARY"

fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_DISTUTILS_INSTALLATION}" ]; then
  
  ptop "GETTING TEXINFO LIBRARY"

  export FOLDER="texinfo-7.0.3"
  
  export FILE="tar.xz"
  
  export URL="https://ftp.gnu.org/gnu/texinfo/${FOLDER}.${FILE}"
  
  export XZF="texinfo.xz"

  wgetact "${FOLDER}" "${FILE}" "${URL}" \
    "${COCOA_TEXINFO_DIR:-"texinfo-7.0.3"}" "${XZF}" || return 1;

  pbottom "GETTING TEXINFO LIBRARY"

fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_OPENBLAS_INSTALLATION}" ]; then

  ptop "GETTING OPENBLAS LIBRARY"

  export URL='https://github.com/OpenMathLib/OpenBLAS.git'

  export FOLDER='OpenBLAS-0.3.23'

  export VER=v0.3.23

  export XZF="OpenBLAS.xz"

  gitact "${FOLDER}" "${VER}" "${URL}" \
    "${COCOA_OPENBLAS_DIR:-"OpenBLAS-0.3.23"}" "${XZF}" || return 1;

  pbottom "GETTING OPENBLAS LIBRARY"

fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_FORTRAN_LAPACK_INSTALLATION}" ]; then

  ptop "GETTING LAPACK LIBRARY"

  export URL='https://github.com/Reference-LAPACK/lapack.git'

  export FOLDER='lapack-3.11.0'

  export VER=v3.11

  export XZF="lapack.xz"

  gitact "${FOLDER}" "${VER}" "${URL}" \
    "${COCOA_LAPACK_DIR:-"lapack-3.11.0"}" "${XZF}" || return 1;

  pbottom "GETTING LAPACK LIBRARY DONE"

fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_HDF5_INSTALLATION}" ]; then

  ptop "GETTING HDF5 LIBRARY"

  export FOLDER="hdf5-1.12.3"

  export FILE="tar.gz"

  export URL_BASE="https://hdf-wordpress-1.s3.amazonaws.com/wp-content/uploads"
  export URL="${URL_BASE:?}/manual/HDF5/HDF5_1_12_3/src/${FOLDER}.${FILE}"

  export XZF="hdf5.xz"

  wgetact "${FOLDER}" "${FILE}" "${URL}" \
    "${COCOA_HDF5_DIR:-"hdf5-1.12.3"}" "${XZF}" || return 1;
  
  pbottom "GETTING HDF5 LIBRARY DONE"

fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_C_CFITSIO_INSTALLATION}" ]; then
  
  ptop "GETTING CFITSIO LIBRARY"

  export FOLDER="cfitsio-4.0.0"
  
  export FILE="tar.gz"
  
  export URL_BASE="http://heasarc.gsfc.nasa.gov"
  export URL="${URL_BASE}/FTP/software/fitsio/c/${FOLDER}.${FILE}"
  
  export XZF="cfitsio.xz"

  wgetact "${FOLDER}" "${FILE}" "${URL}" \ 
    "${COCOA_CFITSIO_DIR:-"cfitsio-4.0.0"}" "${XZF}" || return 1;
  
  pbottom "GETTING CFITSIO LIBRARY DONE"

fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_C_FFTW_INSTALLATION}" ]; then

  ptop "GETTING FFTW LIBRARY"

  export FOLDER="fftw-3.3.10"

  export FILE="tar.gz"

  export URL="http://www.fftw.org/${FOLDER}.${FILE}"

  export XZF="fftw.xz"

  wgetact "${FOLDER}" "${FILE}" "${URL}" \
    "${COCOA_FFTW_DIR:-"fftw-3.3.10"}" "${XZF}" || return 1;

  pbottom "GETTING FFTW LIBRARY DONE"

fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_C_GSL_INSTALLATION}" ]; then

  ptop "GETTING GSL LIBRARY"
  
  export FOLDER="gsl-2.7"

  export FILE="tar.gz"

  export URL="http://ftp.wayne.edu/gnu/gsl/${FOLDER}.${FILE}"

  export XZF="gsl-2.7.xz"

  wgetact "${FOLDER}" "${FILE}" "${URL}" \ 
    "${COCOA_GSL_DIR:-"gsl-2.7"}" "${XZF}" || return 1;

  pbottom "GETTING GSL LIBRARY DONE"

fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_CPP_SPDLOG_INSTALLATION}" ]; then

  ptop "GETTING SPDLOG LIBRARY"

  export URL='https://github.com/gabime/spdlog.git'

  export FOLDER='spdlog'

  export VER=v1.13.0

  export XZF="spdlog.xz"

  gitact "${FOLDER}" "${VER}" "${URL}" \
    "${COCOA_SPDLOG_DIR:-"spdlog"}" "${XZF}" || return 1;

  pbottom "GETTING SPDLOG LIBRARY"

fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_CPP_ARMA_INSTALLATION}" ]; then
  
  ptop "GETTING ARMA LIBRARY DONE"

  export FOLDER="armadillo-12.8.2"
  
  export FILE="tar.xz"
  
  export URL="https://sourceforge.net/projects/arma/files/${FOLDER}.${FILE}"
  
  export XZF="armadillo.xz"

  wgetact "${FOLDER}" "${FILE}" "${URL}" \
    "${COCOA_ARMADILLO_DIR:-"armadillo-12.8.2"}" "${XZF}" || return 1;

  pbottom "GETTING ARMA LIBRARY DONE"

fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_CPP_BOOST_INSTALLATION}" ]; then

  ptop "GETTING BOOST LIBRARY"

  if [ -z "${COCOA_BOOST_DIR}" ]; then
    pfail 'COCOA_BOOST_DIR'; cdroot; return 1;
  fi

  export FOLDER="boost_1_81_0"

  export FILE="tar.gz"

  export URL_BASE="https://boostorg.jfrog.io/artifactory/main/release/"
  export URL="${URL_BASE}/1.81.0/source/${FOLDER}.${FILE}"

  export XZF="boost.xz"

  wgetact "${FOLDER}" "${FILE}" "${URL}" \
    "${COCOA_BOOST_DIR:-"boost_1_81_0"}" "${XZF}" || return 1;

  pbottom "GETTING BOOST LIBRARY"

fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_CPP_CARMA_INSTALLATION}" ]; then

  ptop "GETTING CARMA LIBRARY DONE"

  export CNAME="${COCOA_CARMA_DIR:-"carma":?}"
  
  export PACKDIR="${CCIL:?}/${CNAME:?}"

  export FOLDER="carma_tmp"

  export URL="https://github.com/RUrlus/carma.git"

  export VER=v0.7.0

  export XZF="carma.xz"

  gitact1 "${FOLDER}" "${VER}" "${URL}" || return 1;

  # ---------------------------------------------------------------------------
  # In case this script runs twice (after sudden break w/ CTRL-C)
  rm -rf "${PACKDIR:?}"
  rm -rf "${CCIL:?}/include"

  # -------------------------------------------------------------------------
  # move/rename include file and carma.h folder
  # -------------------------------------------------------------------------
  mv "${CCIL:?}/${FOLDER:?}/include" "${CCIL:?}" >${OUT1:?} 2>${OUT2:?} || 
    { fail_sil "MV CARMA INCLUDE FOLDER"; return 1; }

  mv "${CCIL:?}/include" "${PACKDIR:?}" 2>${OUT2:?} || 
    { fail_sil "RENAME CARMA INCLUDE FOLDER"; return 1; }

  mv "${PACKDIR:?}/carma" "${PACKDIR:?}/carma.h" \
    2>${OUT2:?} || { fail_sil "RENANE CARMA HEADER"; return 1; }

  rm -rf "${CCIL:?}/${FOLDER:?}"
  # -------------------------------------------------------------------------
  
  gitact2 "${CNAME}" "${CNAME}" "${XZF}" || return 1; 

  pbottom "GETTING CARMA LIBRARY DONE"

fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

unset_env_vars_sil || return 1;

pbottom2 "SETUP_INSTALLATION_LIBRARIES"

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------