# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
echo -e '\033[1;44m''SETUP_XZ''\033[0m'

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
if [ -z "${MAKE_NUM_THREADS}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE MAKE_NUM_THREADS IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
fi

# ----------------------------------------------------------------------------
# -------------------------- XZ COMPRESSION Library --------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_XZ_INSTALLATION}" ]; then
  echo -e '\033[1;34m''INSTALLING XZ LIBRARY - \e[4mIT MAY TAKE A WHILE''\033[0m'

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
  cd ../cocoa_installation_libraries/
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
  
  #False xz file: just to trigger GIT LFS
  cp xz-5.2.5.tar.gz.xz xz-5.2.5.tar.gz

  tar -xf xz-5.2.5.tar.gz.xz

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
  cd ./xz-5.2.5/
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

  if [ -z "${DEBUG_XZ_PACKAGE}" ]; then
    export OUTPUT_XZ_1="/dev/null"
    export OUTPUT_XZ_2="/dev/null"
    export XZ_MAKE_NUM_THREADS="${MAKE_NUM_THREADS}"
  else
    export OUTPUT_XZ_1="/dev/tty"
    export OUTPUT_XZ_2="/dev/tty"
    export XZ_MAKE_NUM_THREADS=1
  fi

  CC=$C_COMPILER ./configure \
    --prefix=$ROOTDIR/.local > ${OUTPUT_XZ_1} 2> ${OUTPUT_XZ_2}
  if [ $? -eq 0 ]; then
    echo -e '\033[0;32m'"XZ RUN \e[3mCONFIGURE\e[0m\e\033[0;32m DONE"'\033[0m'
  else
    echo -e '\033[0;31m'"XZ COULD NOT RUN \e[3mCONFIGURE"'\033[0m'
    cd $ROOTDIR
    return 1
  fi

  make -j $XZ_MAKE_NUM_THREADS all > ${OUTPUT_XZ_1} 2> ${OUTPUT_XZ_2}
  if [ $? -eq 0 ]; then
    echo -e '\033[0;32m'"XZ RUN \e[3mMAKE\e[0m\e\033[0;32m DONE"'\033[0m'
  else
    echo -e '\033[0;31m'"XZ COULD NOT RUN \e[3mMAKE"'\033[0m'
    cd $ROOTDIR
    return 1
  fi

  make install > ${OUTPUT_XZ_1} 2> ${OUTPUT_XZ_2}
  if [ $? -eq 0 ]; then
    echo -e '\033[0;32m'"XZ RUN \e[3mMAKE INSTALL\e[0m\e\033[0;32m DONE"'\033[0m'
  else
    echo -e '\033[0;31m'"XZ COULD NOT RUN \e[3mMAKE INSTALL"'\033[0m'
    cd $ROOTDIR
    return 1
  fi

  cd $ROOTDIR
  echo -e '\033[1;34m''\e[4mINSTALLING XZ LIBRARY DONE''\033[0m'
fi

echo -e '\033[1;44m''\e[4mSETUP_XZ DONE''\033[0m'
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------