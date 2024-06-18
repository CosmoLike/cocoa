#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
pfail() {
  echo -e "\033[0;31m ERROR ENV VARIABLE ${1} IS NOT DEFINED \033[0m"
  unset pfail
}
if [ -z "${ROOTDIR}" ]; then
  pfail 'ROOTDIR'
  return 1
fi
if [ -z "${GIT}" ]; then
  pfail 'GIT'
  return 1
fi
unset_env_vars () {
  cd $ROOTDIR
  unset OUT1
  unset OUT2
  unset pfail
  unset URL_BASE
  unset URL
  unset FOLDER
  unset VERSION
  unset XZFILE
  unset FILENAME
  unset wget_action
  unset git_action
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

if [ -z "${DEBUG_SIL_OUTPUT}" ]; then
  if [ -z "${MAKE_NUM_THREADS}" ]; then
    pfail 'MAKE_NUM_THREADS'
    cd $ROOTDIR
    return 1
  fi
  export OUT1="/dev/null"
  export OUT2="/dev/null"
  export SIL_MNT="${MAKE_NUM_THREADS}"
else
  export OUT1="/dev/tty"
  export OUT2="/dev/tty"
  export SIL_MNT=1
fi

wget_action() {
  #ARGUMENTS: FOLDER, FILE, URL, NEW_FOLDER, XZFILE
  cd $ROOTDIR/../cocoa_installation_libraries/

  export FILENAME="${1}.${2}"

  # In case this script runs twice
  rm -rf $ROOTDIR/../cocoa_installation_libraries/$1
  rm -f $ROOTDIR/../cocoa_installation_libraries/$FILENAME
  rm -rf $ROOTDIR/../cocoa_installation_libraries/$4
  rm -f $ROOTDIR/../cocoa_installation_libraries/$5

  wget -q $3 > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "WGET"
    return 1
  fi

  if [ "$2" == "tar.gz" ]; then
    tar zxvf $FILENAME > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail "TAR (UNCOMPRESS tar.gz)"
      return 1
    fi
  elif [ "$2" == "tar.xz" ]; then
    tar xf $FILENAME > ${OUT1} 2> ${OUT2}
    if [ $? -ne 0 ]; then
      fail "TAR (UNCOMPRESS tar.xz)"
      return 1
    fi
  else
    fail "UNKNOWN FILE EXTENSION"
  fi

  if [ "${1}/" != "${4}" ] && [ "${1}" != "${4}" ]; then
    mv "${1}/" $4
    if [ $? -ne 0 ]; then
      fail "MV FOLDER"
      return 1
    fi
  fi

  # Why this compress? In the old Cocoa we saved the libraries in 
  # the github repo using git lfs. So this way preserves the old scripts
  tar -cf - $4 | xz -k -1 --threads=$SIL_MNT -c - > $5
  if [ $? -ne 0 ]; then
    fail "TAR (COMPRESS)"
    return 1
  fi 

  rm -rf $ROOTDIR/../cocoa_installation_libraries/$1
  rm -f $ROOTDIR/../cocoa_installation_libraries/$FILENAME
  rm -rf $ROOTDIR/../cocoa_installation_libraries/$4
}

git_action() {
  #ARGUMENTS: FOLDER, VERSION, URL, NEW_FOLDER, XZFILE
  cd $ROOTDIR/../cocoa_installation_libraries/

  # In case this script runs twice
  rm -rf $ROOTDIR/../cocoa_installation_libraries/$1
  rm -rf $ROOTDIR/../cocoa_installation_libraries/$4
  rm -f $ROOTDIR/../cocoa_installation_libraries/$5

  $GIT clone $3 $1 > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "GIT CLONE"
    return 1
  fi
  
  cd $ROOTDIR/../cocoa_installation_libraries/$1
  
  $GIT checkout $2 > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "GIT CHECKOUT"
    return 1
  fi

  rm -rf ./.git/
  cd $ROOTDIR/../cocoa_installation_libraries/

  if [ "${1}/" != "${4}" ] && [ "${1}" != "${4}" ]; then
    mv "${1}/" $4
    if [ $? -ne 0 ]; then
      fail "MV FOLDER"
      return 1
    fi
  fi

  # Why this compress? In the old Cocoa we saved the libraries in 
  # the github repo using git lfs. So this way preserves the old scripts
  tar -cf - $4 | xz -k -1 --threads=$SIL_MNT -c - > $5
  if [ $? -ne 0 ]; then
    fail "TAR (COMPRESS)"
    return 1
  fi 

  rm -rf $ROOTDIR/../cocoa_installation_libraries/$1
  rm -rf $ROOTDIR/../cocoa_installation_libraries/$4
}
 
ptop2 'SETUP_INSTALLATION_LIBRARIES'

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_CMAKE_INSTALLATION}" ]; then
  ptop "GETTING CMAKE LIBRARY"

  export URL='https://github.com/Kitware/CMake.git'
  export FOLDER='cmake-3.26.4'
  export VERSION=v3.26.4
  export XZFILE="cmake.xz"

  git_action $FOLDER $VERSION $URL $COCOA_CMAKE_DIR $XZFILE

  cd $ROOTDIR
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
  export XZFILE="binutils.xz"

  wget_action $FOLDER, $FILE, $URL, $COCOA_BINUTILS_DIR, $XZFILE

  cd $ROOTDIR
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
  export XZFILE="texinfo.xz"

  wget_action $FOLDER, $FILE, $URL, $COCOA_TEXINFO_DIR $XZFILE

  cd $ROOTDIR
  pbottom "GETTING TEXINFO LIBRARY"
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_OPENBLAS_INSTALLATION}" ]; then
  ptop "GETTING OPENBLAS LIBRARY"

  export URL='https://github.com/OpenMathLib/OpenBLAS.git'
  export FOLDER='OpenBLAS-0.3.23'
  export VERSION=v0.3.23
  export XZFILE="OpenBLAS.xz"

  git_action $FOLDER $VERSION $URL $COCOA_OPENBLAS_DIR $XZFILE

  cd $ROOTDIR
  pbottom "GETTING OPENBLAS LIBRARY"
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_FORTRAN_LAPACK_INSTALLATION}" ]; then
  ptop "GETTING LAPACK LIBRARY"

  export URL='https://github.com/Reference-LAPACK/lapack.git'
  export FOLDER='lapack-3.11.0'
  export VERSION=v3.11
  export XZFILE="lapack.xz"

  git_action $FOLDER $VERSION $URL $COCOA_LAPACK_DIR $XZFILE

  cd $ROOTDIR
  pbottom "GETTING LAPACK LIBRARY DONE"
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_HDF5_INSTALLATION}" ]; then
  ptop "GETTING HDF5 LIBRARY"

  export FOLDER="hdf5-1.12.3"
  export FILE="tar.gz"
  export URL="https://hdf-wordpress-1.s3.amazonaws.com/wp-content/uploads/manual/HDF5/HDF5_1_12_3/src/${FOLDER}.${FILE}"
  export XZFILE="hdf5.xz"

  wget_action $FOLDER, $FILE, $URL, $COCOA_HDF5_DIR, $XZFILE

  cd $ROOTDIR
  pbottom "GETTING HDF5 LIBRARY DONE"
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_C_CFITSIO_INSTALLATION}" ]; then
  ptop "GETTING CFITSIO LIBRARY"

  export FOLDER="cfitsio-4.0.0"
  export FILE="tar.gz"
  export URL="http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/${FOLDER}.${FILE}"
  export XZFILE="cfitsio.xz"

  wget_action $FOLDER, $FILE, $URL, $COCOA_CFITSIO_DIR $XZFILE

  cd $ROOTDIR
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
  export XZFILE="fftw.xz"

  wget_action $FOLDER, $FILE, $URL, $COCOA_FFTW_DIR $XZFILE

  cd $ROOTDIR
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
  export XZFILE="gsl-2.7.xz"

  wget_action $FOLDER, $FILE, $URL, $COCOA_GSL_DIR $XZFILE

  cd $ROOTDIR
  pbottom "GETTING GSL LIBRARY DONE"
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_CPP_SPDLOG_INSTALLATION}" ]; then
  ptop "GETTING SPDLOG LIBRARY"

  export URL='https://github.com/gabime/spdlog.git'
  export FOLDER='spdlog'
  export VERSION=v1.13.0
  export XZFILE="spdlog.xz"

  git_action $FOLDER $VERSION $URL $COCOA_SPDLOG_DIR $XZFILE

  cd $ROOTDIR
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
  export XZFILE="armadillo.xz"

  wget_action $FOLDER $FILE $URL $COCOA_ARMADILLO_DIR $XZFILE

  cd $ROOTDIR
  pbottom "GETTING ARMA LIBRARY DONE"
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_CPP_BOOST_INSTALLATION}" ]; then
  ptop "GETTING BOOST LIBRARY"

  export FOLDER="boost_1_81_0"
  export FILE="tar.gz"
  export URL_BASE="https://boostorg.jfrog.io/artifactory/main/release/"
  export URL="${URL_BASE}/1.81.0/source/${FOLDER}.${FILE}"
  export XZFILE="boost.xz"

  wget_action $FOLDER $FILE $URL $COCOA_BOOST_DIR $XZFILE 

  cd $ROOTDIR
  pbottom "GETTING BOOST LIBRARY"
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_CPP_CARMA_INSTALLATION}" ]; then
  ptop "GETTING CARMA LIBRARY DONE"

  export FOLDER="carma_tmp"
  export URL="https://github.com/RUrlus/carma.git"
  export VERSION=v0.7.0
  export XZFILE="carma.xz"

  cd $ROOTDIR/../cocoa_installation_libraries/

  # In case this script runs twice
  rm -rf $ROOTDIR/../cocoa_installation_libraries/$FOLDER
  rm -rf $ROOTDIR/../cocoa_installation_libraries/$COCOA_CARMA_DIR
  rm -f $ROOTDIR/../cocoa_installation_libraries/$XZFILE

  $GIT clone $URL $FOLDER > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "GIT CLONE"
    return 1
  fi

  cd $ROOTDIR/../cocoa_installation_libraries/$FOLDER

  $GIT checkout $VERSION > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "GIT CHECKOUT"
    return 1
  fi

  mv ./include $ROOTDIR/../cocoa_installation_libraries/
  if [ $? -ne 0 ]; then
    fail "MV CARMA INCLUDE FOLDER"
    return 1
  fi

  cd $ROOTDIR/../cocoa_installation_libraries/

  mv ./include $COCOA_CARMA_DIR
  if [ $? -ne 0 ]; then
    fail "RENANE CARMA INCLUDE FOLDER"
    return 1
  fi

  cd $ROOTDIR/../cocoa_installation_libraries/$COCOA_CARMA_DIR

  mv carma carma.h
  if [ $? -ne 0 ]; then
    fail "RENANE CARMA HEADER"
    return 1
  fi

  cd $ROOTDIR/../cocoa_installation_libraries

  # Why this compress? In the old Cocoa we saved the libraries in 
  # the github repo using git lfs. So this way preserves the old scripts
  tar -cf - $COCOA_CARMA_DIR | xz -k -1 --threads=$SIL_MNT -c - > $XZFILE
  if [ $? -ne 0 ]; then
    fail "TAR (COMPRESS)"
    return 1
  fi 

  rm -rf $ROOTDIR/../cocoa_installation_libraries/$FOLDER
  rm -rf $ROOTDIR/../cocoa_installation_libraries/$COCOA_CARMA_DIR

  cd $ROOTDIR
  pbottom "GETTING CARMA LIBRARY DONE"
fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
unset_env_vars
pbottom2 "SETUP_INSTALLATION_LIBRARIES"
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------