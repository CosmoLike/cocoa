# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
echo -e '\033[1;44m''SETUP_BINUTILS''\033[0m'

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
if [ -z "${FORTRAN_COMPILER}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE FORTRAN_COMPILER IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
fi
if [ -z "${MAKE_NUM_THREADS}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE MAKE_NUM_THREADS IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
fi

if [ -z "${DEBUG_DISTUTILS_PACKAGE}" ]; then
  export OUTPUT_DISTUTILS_PACKAGES_1="/dev/null"
  export OUTPUT_DISTUTILS_PACKAGES_2="/dev/null"
  export DISTUTILS_MAKE_NUM_THREADS="${MAKE_NUM_THREADS}"
else
  export OUTPUT_DISTUTILS_PACKAGES_1="/dev/tty"
  export OUTPUT_DISTUTILS_PACKAGES_2="/dev/tty"
  export DISTUTILS_MAKE_NUM_THREADS=1
fi

# ----------------------------------------------------------------------------
# ----------------------------- TEXINFO LIBRARY  -----------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_DISTUTILS_INSTALLATION}" ]; then
  echo -e '\033[1;34m''INSTALLING TEXINFO LIBRARY - \e[4mIT WILL TAKE A WHILE''\033[0m'

  cd $ROOTDIR/../cocoa_installation_libraries/$COCOA_TEXINFO_DIR

  FC=$FORTRAN_COMPILER CC=$C_COMPILER ./configure --prefix=$ROOTDIR/.local \
    --disable-perl-xs > ${OUTPUT_DISTUTILS_PACKAGES_1} 2> ${OUTPUT_DISTUTILS_PACKAGES_2}
  if [ $? -eq 0 ]; then
    echo -e '\033[0;32m'"TEXINFO RUN \e[3mCONFIGURE\e[0m\e\033[0;32m DONE"'\033[0m'
  else
    echo -e '\033[0;31m'"TEXINFO COULD NOT RUN \e[3mCONFIGURE"'\033[0m'
    cd $ROOTDIR
    return 1
  fi

  make -j $DISTUTILS_MAKE_NUM_THREADS all > ${OUTPUT_DISTUTILS_PACKAGES_1} 2> ${OUTPUT_DISTUTILS_PACKAGES_2}
  if [ $? -eq 0 ]; then
    echo -e '\033[0;32m'"TEXINFO RUN \e[3mMAKE\e[0m\e\033[0;32m DONE"'\033[0m'
  else
    echo -e '\033[0;31m'"TEXINFO COULD NOT RUN \e[3mMAKE"'\033[0m'
    cd $ROOTDIR
    return 1
  fi

  make install > ${OUTPUT_DISTUTILS_PACKAGES_1} 2> ${OUTPUT_DISTUTILS_PACKAGES_2}
  if [ $? -eq 0 ]; then
    echo -e '\033[0;32m'"TEXINFO RUN \e[3mMAKE INSTALL\e[0m\e\033[0;32m DONE"'\033[0m'
  else
    echo -e '\033[0;31m'"TEXINFO COULD NOT RUN \e[3mMAKE INSTALL"'\033[0m'
    cd $ROOTDIR
    return 1
  fi

  cd $ROOTDIR/
  echo -e '\033[1;34m''\e[4mINSTALLING TEXINFO LIBRARY DONE''\033[0m'
fi

# ----------------------------------------------------------------------------
# ----------------------------- DISTUTILS LIBRARY  ---------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_DISTUTILS_INSTALLATION}" ]; then
  echo -e '\033[1;34m''INSTALLING BINUTILS LIBRARY - \e[4mIT WILL TAKE A WHILE''\033[0m'

  cd $ROOTDIR/../cocoa_installation_libraries/$COCOA_BINUTILS_DIR

  FC=$FORTRAN_COMPILER CC=$C_COMPILER ./configure \
    --prefix=$ROOTDIR/.local \
    > ${OUTPUT_DISTUTILS_PACKAGES_1} 2> ${OUTPUT_DISTUTILS_PACKAGES_2}
  if [ $? -eq 0 ]; then
    echo -e '\033[0;32m'"BINUTILS RUN \e[3mCONFIGURE\e[0m\e\033[0;32m DONE"'\033[0m'
  else
    echo -e '\033[0;31m'"BINUTILS COULD NOT RUN \e[3mCONFIGURE"'\033[0m'
    cd $ROOTDIR
    return 1
  fi

  make -j $DISTUTILS_MAKE_NUM_THREADS \
    > ${OUTPUT_DISTUTILS_PACKAGES_1} 2> ${OUTPUT_DISTUTILS_PACKAGES_2}
  if [ $? -eq 0 ]; then
    echo -e '\033[0;32m'"BINUTILS RUN \e[3mMAKE\e[0m\e\033[0;32m DONE"'\033[0m'
  else
    echo -e '\033[0;31m'"BINUTILS COULD NOT RUN \e[3mMAKE"'\033[0m'
    cd $ROOTDIR
    return 1
  fi

  make install > ${OUTPUT_DISTUTILS_PACKAGES_1} 2> ${OUTPUT_DISTUTILS_PACKAGES_2}
  if [ $? -eq 0 ]; then
    echo -e '\033[0;32m'"BINUTILS RUN \e[3mMAKE INSTALL\e[0m\e\033[0;32m DONE"'\033[0m'
  else
    echo -e '\033[0;31m'"BINUTILS COULD NOT RUN \e[3mMAKE INSTALL"'\033[0m'
    cd $ROOTDIR
    return 1
  fi

  cd $ROOTDIR/
  echo -e '\033[1;34m''\e[4mINSTALLING BINUTILS LIBRARY DONE''\033[0m'
fi

echo -e '\033[1;44m''\e[4mSETUP_BINUTILS DONE''\033[0m'
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------