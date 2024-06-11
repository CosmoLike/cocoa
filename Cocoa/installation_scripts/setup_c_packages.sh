#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_C_INSTALLATION}" ]; then
  echo -e '\033[1;44m''SETUP_C_PACKAGES''\033[0m'

  if [ -z "${ROOTDIR}" ]; then
      echo -e '\033[0;31m''ERROR ENV VARIABLE ROOTDIR IS NOT DEFINED''\033[0m'
      return 1
  fi
  if [ -z "${CXX_COMPILER}" ]; then
      echo -e '\033[0;31m''ERROR ENV VARIABLE CXX_COMPILER IS NOT DEFINED''\033[0m'
      cd $ROOTDIR
      return 1
  fi
  if [ -z "${C_COMPILER}" ]; then
      echo -e '\033[0;31m''ERROR ENV VARIABLE C_COMPILER IS NOT DEFINED''\033[0m'
      cd $ROOTDIR
      return 1
  fi
  if [ -z "${CMAKE}" ]; then
      echo -e '\033[0;31m''ERROR ENV VARIABLE MAKE IS NOT DEFINED''\033[0m'
      return 1
  fi
  if [ -z "${FORTRAN_COMPILER}" ]; then
      echo -e '\033[0;31m''ERROR ENV VARIABLE FORTRAN_COMPILER IS NOT DEFINED''\033[0m'
      cd $ROOTDIR
      return 1
  fi
  if [ -z "${PYTHON3}" ]; then
      echo -e '\033[0;31m''ERROR ENV VARIABLE PYTHON3 IS NOT DEFINED''\033[0m'
      cd $ROOTDIR
      return 1
  fi
  if [ -z "${MAKE_NUM_THREADS}" ]; then
      echo -e '\033[0;31m''ERROR ENV VARIABLE MAKE_NUM_THREADS IS NOT DEFINED''\033[0m'
      cd $ROOTDIR
      return 1
  fi

  if [ -z "${DEBUG_C_PACKAGES}" ]; then
    export OUTPUT_CPACKAGES_1="/dev/null"
    export OUTPUT_CPACKAGES_2="/dev/null"
    export C_MAKE_NUM_THREADS="${MAKE_NUM_THREADS}"
  else
    export OUTPUT_CPACKAGES_1="/dev/tty"
    export OUTPUT_CPACKAGES_2="/dev/tty"
    export C_MAKE_NUM_THREADS=1
  fi

  # ----------------------------------------------------------------------------
  # -------------------------------- FFTW --------------------------------------
  # ----------------------------------------------------------------------------
  if [ -z "${IGNORE_C_FFTW_INSTALLATION}" ]; then
    echo -e '\033[1;34m''\t INSTALLING FFTW C LIBRARY - \e[4mIT WILL TAKE A WHILE''\033[0m'

    cd $ROOTDIR/../cocoa_installation_libraries/$COCOA_FFTW_DIR

    FC=$FORTRAN_COMPILER CC=$C_COMPILER ./configure \
      --enable-openmp \
      --prefix=$ROOTDIR/.local \
      --enable-shared=yes \
      --enable-static=yes > ${OUTPUT_CPACKAGES_1} 2> ${OUTPUT_CPACKAGES_2}
    if [ $? -eq 0 ]; then
      echo -e '\033[0;32m'"\t\t FFTW RUN \e[3mCONFIGURE\e[0m\e\033[0;32m DONE"'\033[0m'
    else
      echo -e '\033[0;31m'"FFTW COULD NOT RUN \e[3mCONFIGURE"'\033[0m'
      cd $ROOTDIR
      unset OUTPUT_CPACKAGES_1
      unset OUTPUT_CPACKAGES_2
      unset C_MAKE_NUM_THREADS
      return 1
    fi

    make -j $C_MAKE_NUM_THREADS all > ${OUTPUT_CPACKAGES_1} 2> ${OUTPUT_CPACKAGES_2}
    if [ $? -eq 0 ]; then
      echo -e '\033[0;32m'"\t\t FFTW RUN \e[3mMAKE\e[0m\e\033[0;32m DONE"'\033[0m'
    else
      echo -e '\033[0;31m'"FFTW COULD NOT RUN \e[3mMAKE"'\033[0m'
      cd $ROOTDIR
      unset OUTPUT_CPACKAGES_1
      unset OUTPUT_CPACKAGES_2
      unset C_MAKE_NUM_THREADS
      return 1
    fi

    make install > ${OUTPUT_CPACKAGES_1} 2> ${OUTPUT_CPACKAGES_2}
    if [ $? -eq 0 ]; then
      echo -e '\033[0;32m'"\t\t FFTW RUN \e[3mMAKE INSTALL\e[0m\e\033[0;32m DONE"'\033[0m'
    else
      echo -e '\033[0;31m'"FFTW COULD NOT RUN \e[3mMAKE INSTALL"'\033[0m'
      cd $ROOTDIR
      unset OUTPUT_CPACKAGES_1
      unset OUTPUT_CPACKAGES_2
      unset C_MAKE_NUM_THREADS
      return 1
    fi

    cd $ROOTDIR
    echo -e '\033[1;34m''\t \e[4mINSTALLING FFTW C LIBRARY DONE''\033[0m'
  fi

  # ----------------------------------------------------------------------------
  # ------------------------------- CFITSIO ------------------------------------
  # ----------------------------------------------------------------------------
  if [ -z "${IGNORE_C_CFITSIO_INSTALLATION}" ]; then
    echo -e '\033[1;34m''\t INSTALLING CFITSIO C LIBRARY - \e[4mIT MIGHT TAKE A WHILE''\033[0m'

    cd $ROOTDIR/../cocoa_installation_libraries/$COCOA_CFITSIO_DIR

    rm -f CMakeCache.txt

    mkdir ./CFITSIOBUILD
    cd ./CFITSIOBUILD

    $CMAKE -DBUILD_SHARED_LIBS=TRUE \
      -DCMAKE_INSTALL_PREFIX=$ROOTDIR/.local \
      -DCMAKE_C_COMPILER=$C_COMPILER \
      -DCMAKE_CXX_COMPILER=$CXX_COMPILER \
      -DCMAKE_FC_COMPILER=FORTRAN_COMPILER \
      --log-level=ERROR .. > ${OUTPUT_CPACKAGES_1} 2> ${OUTPUT_CPACKAGES_2}
    if [ $? -eq 0 ]; then
      echo -e '\033[0;32m' "\t\t CFITSIO RUN \e[3mCMAKE\e[0m\e\033[0;32m DONE"'\033[0m'
    else
      echo -e '\033[0;31m'"CFITSIO COULD NOT RUN \e[3mCMAKE"'\033[0m'
      cd $ROOTDIR
      unset OUTPUT_CPACKAGES_1
      unset OUTPUT_CPACKAGES_2
      unset C_MAKE_NUM_THREADS
      return 1
    fi

    make -j $C_MAKE_NUM_THREADS all > ${OUTPUT_CPACKAGES_1} 2> ${OUTPUT_CPACKAGES_2}
    if [ $? -eq 0 ]; then
      echo -e '\033[0;32m'"\t\t CFITSIO RUN \e[3mMAKE\e[0m\e\033[0;32m DONE"'\033[0m'
    else
      echo -e '\033[0;31m'"CFITSIO COULD NOT RUN \e[3mMAKE"'\033[0m'
      cd $ROOTDIR
      unset OUTPUT_CPACKAGES_1
      unset OUTPUT_CPACKAGES_2
      unset C_MAKE_NUM_THREADS
      return 1
    fi

    make install > ${OUTPUT_CPACKAGES_1} 2> ${OUTPUT_CPACKAGES_2}
    if [ $? -eq 0 ]; then
      echo -e '\033[0;32m' "\t\t CFITSIO RUN \e[3mMAKE INSTALL\e[0m\e\033[0;32m DONE"'\033[0m'
    else
      echo -e '\033[0;31m'"CFITSIO COULD NOT RUN \e[3mMAKE INSTALL"'\033[0m'
      cd $ROOTDIR
      unset OUTPUT_CPACKAGES_1
      unset OUTPUT_CPACKAGES_2
      unset C_MAKE_NUM_THREADS
      return 1
    fi

    cd ../
    rm -rf ./CFITSIOBUILD

    cd $ROOTDIR
    echo -e '\033[1;34m''\t \e[4mINSTALLING CFITSIO C LIBRARY DONE''\033[0m'
  fi

  # ----------------------------------------------------------------------------
  # ----------------------------------------- GSL ------------------------------
  # ----------------------------------------------------------------------------
  if [ -z "${IGNORE_C_GSL_INSTALLATION}" ]; then
    echo -e '\033[1;34m''\t INSTALLING GSL C LIBRARY - \e[4mIT WILL TAKE A WHILE''\033[0m'

    cd $ROOTDIR/../cocoa_installation_libraries/$COCOA_GSL_DIR

    CC=$C_COMPILER ./configure \
      --prefix=$ROOTDIR/.local \
      --enable-shared=yes \
      --enable-static=yes > ${OUTPUT_CPACKAGES_1} 2> ${OUTPUT_CPACKAGES_2}
    if [ $? -eq 0 ]; then
      echo -e '\033[0;32m' "\t\t GSL RUN \e[3mCONFIGURE\e[0m\e\033[0;32m DONE"'\033[0m'
    else
      echo -e '\033[0;31m'"GSL COULD NOT RUN \e[3mCONFIGURE"'\033[0m'
      cd $ROOTDIR
      unset OUTPUT_CPACKAGES_1
      unset OUTPUT_CPACKAGES_2
      unset C_MAKE_NUM_THREADS
      return 1
    fi

    make -j $C_MAKE_NUM_THREADS all > ${OUTPUT_CPACKAGES_1} 2> ${OUTPUT_CPACKAGES_2}
    if [ $? -eq 0 ]; then
      echo -e '\033[0;32m'"\t\t GSL RUN \e[3mMAKE\e[0m\e\033[0;32m DONE"'\033[0m'
    else
      echo -e '\033[0;31m'"GSL COULD NOT RUN \e[3mMAKE"'\033[0m'
      cd $ROOTDIR
      unset OUTPUT_CPACKAGES_1
      unset OUTPUT_CPACKAGES_2
      unset C_MAKE_NUM_THREADS
      return 1
    fi

    make install > ${OUTPUT_CPACKAGES_1} 2> ${OUTPUT_CPACKAGES_2}
    if [ $? -eq 0 ]; then
      echo -e '\033[0;32m'"\t\t GSL RUN \e[3mMAKE INSTALL\e[0m\e\033[0;32m DONE"'\033[0m'
    else
      echo -e '\033[0;31m'"GSL COULD NOT RUN \e[3mMAKE INSTALL"'\033[0m'
      cd $ROOTDIR
      unset OUTPUT_CPACKAGES_1
      unset OUTPUT_CPACKAGES_2
      unset C_MAKE_NUM_THREADS
      return 1
    fi

    cd $ROOTDIR
    echo -e '\033[1;34m''\t \e[4mINSTALLING GSL C LIBRARY DONE''\033[0m'
  fi
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  cd $ROOTDIR
  unset OUTPUT_CPACKAGES_1
  unset OUTPUT_CPACKAGES_2
  unset C_MAKE_NUM_THREADS
  echo -e '\033[1;44m''\e[4mSETUP_C_PACKAGES DONE''\033[0m'
fi

# ------------------------------------------------------------------------------
# ----------------------------------- EUCLID EMU -------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_ALL_PIP_INSTALLATION}" ]; then
  if [ -z "${DEBUG_PIP_OUTPUT}" ]; then
    export OUTPUT_PIP_1="/dev/null"
    export OUTPUT_PIP_2="/dev/null"
  else
    export OUTPUT_PIP_1="/dev/tty"
    export OUTPUT_PIP_2="/dev/tty"
  fi

  echo -e '\033[1;34m''INSTALLING EUCLIDEMU2''\033[0m'

  # WE MIGRATED euclidemu2 TO setup_c_packages SCRIPT BECAUSE IT DEPENDS ON GSL-GNU LIB
  env CXX=$CXX_COMPILER CC=$C_COMPILER  $PIP3 install --global-option=build_ext \
    $ROOTDIR/../cocoa_installation_libraries/euclidemu2-1.2.0 \
    --no-dependencies \
    --prefix=$ROOTDIR/.local --no-index > ${OUTPUT_PIP_1} 2> ${OUTPUT_PIP_2}
  if [ $? -ne 0 ]; then
    echo -e '\033[0;31m'"PIP COULD NOT RUN \e[3mPIP INSTALL EUCLUDEMUL"'\033[0m'
    cd $ROOTDIR
    unset OUTPUT_PIP_1
    unset OUTPUT_PIP_2
    return 1
  else
    echo -e '\033[0;32m'"\t\t PIP RUN \e[3mPIP INSTALL EUCLUDEMUL\e[0m\e\033[0;32m DONE"'\033[0m'
  fi

  cd $ROOTDIR
  unset OUTPUT_PIP_1
  unset OUTPUT_PIP_2
  echo -e '\033[1;34m''\e[4mINSTALLING EUCLIDEMU2 DONE''\033[0m'
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------