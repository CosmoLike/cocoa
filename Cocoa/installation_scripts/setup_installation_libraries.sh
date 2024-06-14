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
  unset CMAKE_URL
  unset ARMA_FILE
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

ptop2 'SETUP_INSTALLATION_LIBRARIES'

if [ -z "${IGNORE_CMAKE_INSTALLATION}" ]; then
  ptop "GETTING CMAKE LIBRARY"

  cd $ROOTDIR/../cocoa_installation_libraries/
  if [ $? -ne 0 ]; then
    fail "CD COCOA_INSTALLATION_LIBRARIES"
    return 1
  fi

  # In case this script runs twice
  rm -rf $ROOTDIR/../cocoa_installation_libraries/cmake-3.26.4
  rm -f $ROOTDIR/../cocoa_installation_libraries/cmake.xz

  export CMAKE_URL='https://github.com/Kitware/CMake.git'
  export CMAKE_FOLDER='cmake-3.26.4'
  $GIT clone $CMAKE_URL $CMAKE_FOLDER > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "GIT CLONE"
    return 1
  fi
  
  cd $ROOTDIR/../cocoa_installation_libraries/$CMAKE_FOLDER
  
  $GIT checkout v3.26.4 > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "GIT CHECKOUT"
    return 1
  fi

  rm -rf $ROOTDIR/../cocoa_installation_libraries/$CMAKE_FOLDER/.git/

  cd $ROOTDIR/../cocoa_installation_libraries/
  
  # Why this compress? In the old Cocoa we saved the libraries in 
  # the github repo using git lfs. So this way preserves the old scripts
  tar -cf - $CMAKE_FOLDER | xz -k -1 --threads=$SIL_MNT -c - > cmake.xz
  if [ $? -ne 0 ]; then
    fail "TAR"
    return 1
  fi

  rm -rf $ROOTDIR/../cocoa_installation_libraries/$CMAKE_FOLDER
  
  cd $ROOTDIR

  pbottom "GETTING CMAKE LIBRARY"
fi

if [ -z "${IGNORE_DISTUTILS_INSTALLATION}" ]; then
  echo -e '\033[0;32m'"\t\tGETTING BINUTILS LIBRARY"'\033[0m' > ${OUT1} 2> ${OUT2}

  cd $ROOTDIR/../cocoa_installation_libraries/

  # In case this script runs twice
  rm -f $ROOTDIR/../cocoa_installation_libraries/binutils-2.37.tar.gz
  rm -f $ROOTDIR/../cocoa_installation_libraries/binutils.xz

  wget -q https://ftp.gnu.org/gnu/binutils/binutils-2.37.tar.gz
  if [ $? -ne 0 ]; then
    fail "WGET"
    return 1
  fi

  tar zxvf binutils-2.37.tar.gz
  if [ $? -ne 0 ]; then
    fail "TAR (UNCOMPRESS)"
    return 1
  fi

  rm -f $ROOTDIR/../cocoa_installation_libraries/binutils-2.37.tar.gz

  # Why this compress? In the old Cocoa we saved the libraries in 
  # the github repo using git lfs. So this way preserves the old scripts
  tar -cf - binutils-2.37/ | xz -k -1 --threads=$SIL_MNT -c - > binutils.xz
  if [ $? -ne 0 ];then
    echo -e '\033[0;31m'"BINUTILS: COULD NOT COMPRESS \e[3mBINUTILS FOLDER"'\033[0m'
    unset_env_vars
    unset unset_env_vars
    return 1
  fi
  rm -rf $ROOTDIR/../cocoa_installation_libraries/binutils-2.37

  cd $ROOTDIR

  echo -e '\033[0;32m'"\t\tGETTING BINUTILS LIBRARY DONE"'\033[0m' > ${OUT1} 2> ${OUT2}
fi

if [ -z "${IGNORE_DISTUTILS_INSTALLATION}" ]; then
  echo -e '\033[0;32m'"\t\tGETTING TEXINFO LIBRARY"'\033[0m' > ${OUT1} 2> ${OUT2}

  cd $ROOTDIR/../cocoa_installation_libraries/

  # In case this script runs twice
  rm -f $ROOTDIR/../cocoa_installation_libraries/texinfo-7.0.3.tar.xz
  rm -f $ROOTDIR/../cocoa_installation_libraries/texinfo.xz

  wget -q https://ftp.gnu.org/gnu/texinfo/texinfo-7.0.3.tar.xz
  if [ $? -ne 0 ]; then
    fail "WGET"
    return 1
  fi

  tar xf texinfo-7.0.3.tar.xz
  if [ $? -ne 0 ]; then
    fail "TAR (UNCOMPRESS)"
    return 1
  fi

  rm -f $ROOTDIR/../cocoa_installation_libraries/texinfo-7.0.3.tar.xz

  # Why this compress? In the old Cocoa we saved the libraries in 
  # the github repo using git lfs. So this way preserves the old scripts
  tar -cf - texinfo-7.0.3/ | xz -k -1 --threads=$SIL_MNT -c - > texinfo.xz
  if [ $? -ne 0 ];then
    echo -e '\033[0;31m'"TEXINFO: COULD NOT COMPRESS \e[3mTEXINFO FOLDER"'\033[0m'
    unset_env_vars
    unset unset_env_vars
    return 1
  fi
  rm -rf $ROOTDIR/../cocoa_installation_libraries/texinfo-7.0.3

  cd $ROOTDIR

  echo -e '\033[0;32m'"\t\tGETTING TEXINFO LIBRARY DONE"'\033[0m' > ${OUT1} 2> ${OUT2}
fi

if [ -z "${IGNORE_OPENBLAS_INSTALLATION}" ]; then
  echo -e '\033[0;32m'"\t\tGETTING OPENBLAS LIBRARY"'\033[0m' > ${OUT1} 2> ${OUT2}
  
  cd $ROOTDIR/../cocoa_installation_libraries/

  # In case this script runs twice
  rm -rf $ROOTDIR/../cocoa_installation_libraries/OpenBLAS-0.3.23
  rm -f $ROOTDIR/../cocoa_installation_libraries/OpenBLAS.xz

  $GIT clone https://github.com/OpenMathLib/OpenBLAS.git OpenBLAS-0.3.23 > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ];then
    echo -e '\033[0;31m'"OPENBLAS: COULD NOT RUN \e[3mGIT CLONE"'\033[0m'
    unset_env_vars
    unset unset_env_vars
    return 1   
  fi

  cd $ROOTDIR/../cocoa_installation_libraries/OpenBLAS-0.3.23

  $GIT checkout v0.3.23 > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ];then
    echo -e '\033[0;31m'"OPENBLAS: COULD NOT RUN \e[3mGIT CHECKOUT"'\033[0m'
    unset_env_vars
    unset unset_env_vars
    return 1
  fi

  rm -rf $ROOTDIR/../cocoa_installation_libraries/OpenBLAS-0.3.23/.git/

  # Why this compress? In the old Cocoa we saved the libraries in 
  # the github repo using git lfs. So this way preserves the old scripts
  tar -cf - OpenBLAS-0.3.23/ | xz -k -1 --threads=$SIL_MNT -c - > OpenBLAS.xz
  if [ $? -ne 0 ];then
    echo -e '\033[0;31m'"OPENBLAS: COULD NOT COMPRESS \e[3mOPENBLAS FOLDER"'\033[0m'
    unset_env_vars
    unset unset_env_vars
    return 1
  fi
  rm -rf $ROOTDIR/../cocoa_installation_libraries/OpenBLAS-0.3.23

  cd $ROOTDIR

  echo -e '\033[0;32m'"\t\tGETTING OPENBLAS LIBRARY DONE"'\033[0m' > ${OUT1} 2> ${OUT2}
fi

if [ -z "${IGNORE_FORTRAN_LAPACK_INSTALLATION}" ]; then
  echo -e '\033[0;32m'"\t\tGETTING LAPACK LIBRARY"'\033[0m' > ${OUT1} 2> ${OUT2}
  
  cd $ROOTDIR/../cocoa_installation_libraries/

  # In case this script runs twice
  rm -rf $ROOTDIR/../cocoa_installation_libraries/lapack
  rm -f $ROOTDIR/../cocoa_installation_libraries/lapack.xz

  $GIT clone https://github.com/Reference-LAPACK/lapack.git lapack > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ];then
    echo -e '\033[0;31m'"LAPACK: COULD NOT RUN \e[3mGIT CLONE"'\033[0m'
    unset_env_vars
    unset unset_env_vars
    return 1
  fi

  cd $ROOTDIR/../cocoa_installation_libraries/lapack-3.11.0
  
  $GIT checkout lapack-3.11.0 > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ];then
    echo -e '\033[0;31m'"LAPACK: COULD NOT RUN \e[3mGIT CHECKOUT"'\033[0m'
    unset_env_vars
    unset unset_env_vars
    return 1
  fi

  rm -rf $ROOTDIR/../cocoa_installation_libraries/lapack-3.11.0/.git/

  cd $ROOTDIR/../cocoa_installation_libraries/

  # Why this compress? In the old Cocoa we saved the libraries in 
  # the github repo using git lfs. So this way preserves the old scripts
  tar -cf - lapack-3.11.0/ | xz -k -1 --threads=$SIL_MNT -c - > lapack.xz
  if [ $? -ne 0 ];then
    echo -e '\033[0;31m'"LAPACK: COULD NOT COMPRESS \e[3mLAPACK FOLDER"'\033[0m'
    unset_env_vars
    unset unset_env_vars
    return 1
  fi
  rm -rf $ROOTDIR/../cocoa_installation_libraries/lapack-3.11.0

  cd $ROOTDIR

  echo -e '\033[0;32m'"\t\tGETTING LAPACK LIBRARY DONE"'\033[0m' > ${OUT1} 2> ${OUT2}
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_HDF5_INSTALLATION}" ]; then
  echo -e '\033[0;32m'"\t\tGETTING HDF5 LIBRARY"'\033[0m' > ${OUT1} 2> ${OUT2}

  cd $ROOTDIR/../cocoa_installation_libraries/

  # In case this script runs twice
  rm -f $ROOTDIR/../cocoa_installation_libraries/hdf5-1.12.3.tar.gz
  rm -f $ROOTDIR/../cocoa_installation_libraries/hdf5.xz

  wget -q https://hdf-wordpress-1.s3.amazonaws.com/wp-content/uploads/manual/HDF5/HDF5_1_12_3/src/hdf5-1.12.3.tar.gz
  if [ $? -ne 0 ]; then
    fail "WGET"
    return 1
  fi

  tar zxvf hdf5-1.12.3.tar.gz
  if [ $? -ne 0 ]; then
    fail "TAR (UNCOMPRESS)"
    return 1
  fi

  rm -f $ROOTDIR/../cocoa_installation_libraries/hdf5-1.12.3.tar.gz

  # Why this compress? In the old Cocoa we saved the libraries in 
  # the github repo using git lfs. So this way preserves the old scripts
  tar -cf - hdf5-1.12.3/ | xz -k -1 --threads=$SIL_MNT -c - > hdf5.xz
  if [ $? -ne 0 ];then
    echo -e '\033[0;31m'"HDF5: COULD NOT COMPRESS \e[3mHDF5 FOLDER"'\033[0m'
    unset_env_vars
    unset unset_env_vars
    return 1
  fi
  rm -rf $ROOTDIR/../cocoa_installation_libraries/hdf5-1.12.3

  cd $ROOTDIR

  echo -e '\033[0;32m'"\t\tGETTING HDF5 LIBRARY DONE"'\033[0m' > ${OUT1} 2> ${OUT2}
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_C_CFITSIO_INSTALLATION}" ]; then
  echo -e '\033[0;32m'"\t\tGETTING CFITSIO LIBRARY"'\033[0m' > ${OUT1} 2> ${OUT2}

  cd $ROOTDIR/../cocoa_installation_libraries/

  # In case this script runs twice
  rm -f $ROOTDIR/../cocoa_installation_libraries/cfitsio-4.0.0.tar.gz
  rm -f $ROOTDIR/../cocoa_installation_libraries/cfitsio.xz

  wget -q http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-4.0.0.tar.gz
  if [ $? -ne 0 ]; then
    fail "WGET"
    return 1
  fi

  tar zxvf cfitsio-4.0.0.tar.gz
  if [ $? -ne 0 ]; then
    fail "TAR (UNCOMPRESS)"
    return 1
  fi

  rm -f $ROOTDIR/../cocoa_installation_libraries/cfitsio-4.0.0.tar.gz

  # Why this compress? In the old Cocoa we saved the libraries in 
  # the github repo using git lfs. So this way preserves the old scripts
  tar -cf - cfitsio-4.0.0/ | xz -k -1 --threads=$SIL_MNT -c - > cfitsio.xz
  if [ $? -ne 0 ];then
    echo -e '\033[0;31m'"CFITSIO: COULD NOT COMPRESS \e[3mCFITSIO FOLDER"'\033[0m'
    unset_env_vars
    unset unset_env_vars
    return 1
  fi
  rm -rf $ROOTDIR/../cocoa_installation_libraries/cfitsio-4.0.0

  cd $ROOTDIR

  echo -e '\033[0;32m'"\t\tGETTING CFITSIO LIBRARY DONE"'\033[0m' > ${OUT1} 2> ${OUT2}
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_C_FFTW_INSTALLATION}" ]; then
  echo -e '\033[0;32m'"\t\tGETTING FFTW LIBRARY"'\033[0m' > ${OUT1} 2> ${OUT2}
  
  cd $ROOTDIR/../cocoa_installation_libraries/

  # In case this script runs twice
  rm -f $ROOTDIR/../cocoa_installation_libraries/fftw-3.3.10.tar.gz 
  rm -f $ROOTDIR/../cocoa_installation_libraries/fftw.xz

  wget -q http://www.fftw.org/fftw-3.3.10.tar.gz 
  if [ $? -ne 0 ]; then
    fail "WGET"
    return 1
  fi

  tar zxvf fftw-3.3.10.tar.gz
  if [ $? -ne 0 ]; then
    fail "TAR (UNCOMPRESS)"
    return 1
  fi

  rm -f fftw-3.3.10.tar.gz

  # Why this compress? In the old Cocoa we saved the libraries in 
  # the github repo using git lfs. So this way preserves the old scripts
  tar -cf - fftw-3.3.10/ | xz -k -1 --threads=$SIL_MNT -c - > fftw.xz
  if [ $? -ne 0 ];then
    echo -e '\033[0;31m'"FFTW: COULD NOT COMPRESS \e[3mFFTW FOLDER"'\033[0m'
    unset_env_vars
    unset unset_env_vars
    return 1
  fi
  rm -rf $ROOTDIR/../cocoa_installation_libraries/fftw-3.3.10

  cd $ROOTDIR

  echo -e '\033[0;32m'"\t\tGETTING FFTW LIBRARY DONE"'\033[0m' > ${OUT1} 2> ${OUT2}
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_C_GSL_INSTALLATION}" ]; then
  echo -e '\033[0;32m'"\t\tGETTING GSL LIBRARY"'\033[0m' > ${OUT1} 2> ${OUT2}

  cd $ROOTDIR/../cocoa_installation_libraries/

  # In case this script runs twice
  rm -f $ROOTDIR/../cocoa_installation_libraries/gsl-2.7.tar.gz
  rm -f $ROOTDIR/../cocoa_installation_libraries/gsl.xz

  wget -q http://ftp.wayne.edu/gnu/gsl/gsl-2.7.tar.gz
  if [ $? -ne 0 ]; then
    fail "WGET"
    return 1
  fi

  tar zxvf gsl-2.7.tar.gz
  if [ $? -ne 0 ]; then
    fail "TAR (UNCOMPRESS)"
    return 1
  fi

  rm -f ROOTDIR/../cocoa_installation_libraries/gsl-2.7.tar.gz

  # Why this compress? In the old Cocoa we saved the libraries in 
  # the github repo using git lfs. So this way preserves the old scripts
  tar -cf - gsl-2.7/ | xz -k -1 --threads=$SIL_MNT -c - > gsl.xz
  if [ $? -ne 0 ];then
    echo -e '\033[0;31m'"GSL: COULD NOT COMPRESS \e[3mGSL FOLDER"'\033[0m'
    unset_env_vars
    unset unset_env_vars
    return 1
  fi
  rm -rf $ROOTDIR/../cocoa_installation_libraries/gsl-2.7

  cd $ROOTDIR

  echo -e '\033[0;32m'"\t\tGETTING GSL LIBRARY DONE"'\033[0m' > ${OUT1} 2> ${OUT2}
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_CPP_SPDLOG_INSTALLATION}" ]; then
  ptop "GETTING SPDLOG LIBRARY"
  
  cd $ROOTDIR/../cocoa_installation_libraries/

  # In case this script runs twice
  rm -rf $ROOTDIR/../cocoa_installation_libraries/spdlog
  rm -f $ROOTDIR/../cocoa_installation_libraries/spdlog.xz

  $GIT clone https://github.com/gabime/spdlog.git spdlog > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "GIT CLONE"
    return 1
  fi

  cd $ROOTDIR/../cocoa_installation_libraries/spdlog

  $GIT checkout v1.13.0 > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "GIT CHECKOUT"
    return 1
  fi

  rm -rf $ROOTDIR/../cocoa_installation_libraries/spdlog/.git/

  cd $ROOTDIR/../cocoa_installation_libraries/

  # Why this compress? In the old Cocoa we saved the libraries in 
  # the github repo using git lfs. So this way preserves the old scripts
  tar -cf - spdlog/ | xz -k -1 --threads=$SIL_MNT -c - > spdlog.xz
  if [ $? -ne 0 ];then
    echo -e '\033[0;31m'"SPDLOG: COULD NOT COMPRESS \e[3mSPDLOG FOLDER"'\033[0m'
    unset_env_vars
    unset unset_env_vars
    return 1
  fi
  rm -rf $ROOTDIR/../cocoa_installation_libraries/spdlog

  cd $ROOTDIR

  pbottom "GETTING SPDLOG LIBRARY"
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_CPP_ARMA_INSTALLATION}" ]; then
  ptop "GETTING ARMA LIBRARY DONE"

  export ARMA_FILE="armadillo-12.8.2"

  cd $ROOTDIR/../cocoa_installation_libraries/

  # In case this script runs twice
  rm -f $ROOTDIR/../cocoa_installation_libraries/$ARMA_FILE.tar.xz
  rm -f $ROOTDIR/../cocoa_installation_libraries/armadillo.xz

  wget -q https://sourceforge.net/projects/arma/files/$ARMA_FILE.tar.xz
  if [ $? -ne 0 ]; then
    fail "WGET"
    return 1
  fi

  tar xf $ARMA_FILE.tar.xz
  if [ $? -ne 0 ]; then
    fail "TAR (UNCOMPRESS)"
    return 1
  fi

  rm -f ROOTDIR/../cocoa_installation_libraries/$ARMA_FILE.tar.xz

  mv "${ARMA_FILE}/" $COCOA_ARMADILLO_DIR
  if [ $? -ne 0 ]; then
    fail "MV FOLDER"
    return 1
  fi

  # Why this compress? In the old Cocoa we saved the libraries in 
  # the github repo using git lfs. So this way preserves the old scripts
  tar -cf - $COCOA_ARMADILLO_DIR | xz -k -1 --threads=$SIL_MNT -c - > armadillo.xz
  if [ $? -ne 0 ]; then
    fail "TAR (COMPRESS)"
    return 1
  fi

  rm -rf $ROOTDIR/../cocoa_installation_libraries/$COCOA_ARMADILLO_DIR

  cd $ROOTDIR
  pbottom "GETTING ARMA LIBRARY DONE"
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_CPP_BOOST_INSTALLATION}" ]; then
  ptop "GETTING BOOST LIBRARY"

  export BOOST_FILE="boost_1_81_0"
  export BOOST_URL_BASE="https://boostorg.jfrog.io/artifactory/main/release/"
  export BOOST_URL="${BOOST_URL}/1.81.0/source/"

  cd $ROOTDIR/../cocoa_installation_libraries/

  # In case this script runs twice
  rm -f $ROOTDIR/../cocoa_installation_libraries/$BOOST_FILE.tar.gz
  rm -f $ROOTDIR/../cocoa_installation_libraries/boost.xz

  wget -q "${BOOST_URL}/${BOOST_FILE}.tar.gz"
  if [ $? -ne 0 ]; then
    fail "WGET"
    return 1
  fi
  
  tar xvzf $BOOST_FILE.tar.gz
  if [ $? -ne 0 ]; then
    fail "TAR (UNCOMPRESS)"
    return 1
  fi

  rm -f ROOTDIR/../cocoa_installation_libraries/$BOOST_FILE.tar.gz

  mv "${BOOST_FILE}/" $COCOA_BOOST_DIR
  if [ $? -ne 0 ]; then
    fail "MV FOLDER"
    return 1
  fi

  # Why this compress? In the old Cocoa we saved the libraries in 
  # the github repo using git lfs. So this way preserves the old scripts
  tar -cf - $COCOA_BOOST_DIR/ | xz -k -1 --threads=$SIL_MNT -c - > boost.xz
  if [ $? -ne 0 ]; then
    fail "TAR (COMPRESS)"
    return 1
  fi
  
  rm -rf $ROOTDIR/../cocoa_installation_libraries/$COCOA_BOOST_DIR

  cd $ROOTDIR
  pbottom "GETTING BOOST LIBRARY"
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_CPP_CARMA_INSTALLATION}" ]; then
  ptop "GETTING CARMA LIBRARY DONE"

  export CARMA_URL="https://github.com/RUrlus/carma.git"
  export CARMA_VERSION=v0.7.0

  cd $ROOTDIR/../cocoa_installation_libraries/

  # In case this script runs twice
  rm -rf $ROOTDIR/../cocoa_installation_libraries/carma_tmp
  rm -rf $ROOTDIR/../cocoa_installation_libraries/$COCOA_CARMA_DIR
  rm -f $ROOTDIR/../cocoa_installation_libraries/carma.xz

  $GIT clone $CARMA_URL carma_tmp > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail "GIT CLONE"
    return 1
  fi

  cd $ROOTDIR/../cocoa_installation_libraries/carma_tmp

  $GIT checkout $CARMA_VERSION > ${OUT1} 2> ${OUT2}
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

  rm -rf $ROOTDIR/../cocoa_installation_libraries/carma_tmp

  cd $ROOTDIR/../cocoa_installation_libraries/$COCOA_CARMA_DIR

  mv carma carma.h
  if [ $? -ne 0 ]; then
    fail "RENANE CARMA HEADER"
    return 1
  fi

  cd $ROOTDIR/../cocoa_installation_libraries

  # Why this compress? In the old Cocoa we saved the libraries in 
  # the github repo using git lfs. So this way preserves the old scripts
  tar -cf - $COCOA_CARMA_DIR | xz -k -1 --threads=$SIL_MNT -c - > carma.xz
  if [ $? -ne 0 ]; then
    fail "TAR (COMPRESS)"
    return 1
  fi 

  rm -rf $ROOTDIR/../cocoa_installation_libraries/$COCOA_CARMA_DIR

  cd $ROOTDIR
  pbottom "GETTING CARMA LIBRARY DONE"
fi

unset_env_vars
unset unset_env_vars
echo -e '\033[1;34m''\t \e[4mSETUP_INSTALLATION_LIBRARIES DONE''\033[0m'
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------