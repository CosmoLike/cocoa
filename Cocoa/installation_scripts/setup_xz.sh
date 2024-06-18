#!/bin/bash
# ----------------------------------------------------------------------------
# -------------------------- XZ COMPRESSION Library --------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_XZ_INSTALLATION}" ]; then
  pfail() {
    echo -e "\033[0;31m ERROR ENV VARIABLE ${1} IS NOT DEFINED \033[0m"
    unset pfail
  }
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'
    return 1
  fi
  if [ -z "${C_COMPILER}" ]; then
    pfail 'C_COMPILER'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${MAKE_NUM_THREADS}" ]; then
    pfail 'MAKE_NUM_THREADS'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${DEBUG_XZ_PACKAGE}" ]; then
    export OUT1="/dev/null"
    export OUT2="/dev/null"
    export XZ_MNT="${MAKE_NUM_THREADS}"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
    export XZ_MNT=1
  fi
  unset_env_vars_sxz () {
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset POLY_URL
    unset POLY_CHANGES
    unset pfail
    unset unset_env_vars_spoly
  }
  fail_sxz () {
    export MSG="\033[0;31m (setup_xz.sh) WE CANNOT RUN \e[3m"
    export MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"  
    unset_env_vars_sxz
    unset MSG
    unset MSG2
    unset fail_sxz
  }

  ptop2 "SETUP_XZ"
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  ptop "INSTALLING XZ LIBRARY"

  cd ${ROOTDIR}/../cocoa_installation_libraries/
  if [ $? -ne 0 ]; then
    fail_sxz "CD COCOA_INSTALLATION_LIBRARIES"
    return 1
  fi

  #False xz file: just to trigger GIT LFS
  cp xz-5.2.5.tar.gz.xz xz-5.2.5.tar.gz
  if [ $? -ne 0 ]; then
    fail_sxz "CP XZ TAR"
    return 1
  fi

  tar -xf xz-5.2.5.tar.gz.xz
  if [ $? -ne 0 ]; then
    fail_sxz "TAR XZ TAR"
    return 1
  fi

  cd ./xz-5.2.5/
   if [ $? -ne 0 ]; then
    fail_sxz "CD XZ FOLDEr"
    return 1
  fi 

  CC=$C_COMPILER ./configure --prefix=$ROOTDIR/.local > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_sxz "CONFIGURE"
    return 1
  fi

  make -j $XZ_MNT all > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_sxz "MAKE"
    return 1
  fi

  make install > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_sxz "MAKE INSTALL"
    return 1
  fi

  cd $ROOTDIR
  pbottom "INSTALLING XZ LIBRARY"
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  pbottom2 "SETUP_XZ"
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------