#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_CAMB_COMPILATION}" ]; then
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
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${CAMB_GIT_COMMIT}" ]; then
    pfail 'CAMB_GIT_COMMIT'
    cd $ROOTDIR
    return 1
  fi
  if [ -z "${CAMB_NAME}" ]; then
    pfail 'CAMB_NAME'
    cd $ROOTDIR
    return 1
  fi
  unset_env_vars_scamb () {
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset URL
    unset CHANGES
    unset pfail
    unset unset_env_vars_scamb
  }
  fail_scamb () {
    export MSG="\033[0;31m (setup_camb.sh) WE CANNOT RUN \e[3m"
    export MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"
    unset_env_vars_scamb
    unset MSG
    unset MSG2
    unset fail_scamb
  }
  if [ -z "${DEBUG_CAMB_OUTPUT}" ]; then
    export OUT1="/dev/null"
    export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
  fi
  
  # ---------------------------------------------------------------------------
  ptop2 'SETUP_CAMB'
  ptop 'INSTALLING CAMB'

  export URL="https://github.com/cmbant/CAMB"
  export CHANGES="${ROOTDIR}/../cocoa_installation_libraries/camb_changes"

  # ---------------------------------------------------------------------------
  # in case this script is called twice
  # ---------------------------------------------------------------------------
  rm -rf $ROOTDIR/external_modules/code/$CAMB_NAME

  # ---------------------------------------------------------------------------
  # clone from original repo
  # ---------------------------------------------------------------------------
  cd $ROOTDIR/external_modules/code/

  $GIT clone $URL --recursive $CAMB_NAME > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_scamb "GIT CLONE FROM CAMB REPO"
    return 1
  fi

  cd $ROOTDIR/external_modules/code/$CAMB_NAME

  $GIT checkout $CAMB_GIT_COMMIT > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_scamb "GIT CHECKOUT CAMB"
    return 1
  fi
  
  # ---------------------------------------------------------------------------
  # patch CAMB to be compatible w/ COCOA
  # ---------------------------------------------------------------------------
  cd $ROOTDIR/external_modules/code/$CAMB_NAME/camb/
  if [ $? -ne 0 ]; then
    fail_scamb "CD CAMB CAMB FOLDER"
    return 1
  fi
  cp $CHANGES/camb/_compilers.patch .
  if [ $? -ne 0 ]; then
    fail_scamb "CP FILE PATCH (_compilers)"
    return 1
  fi
  patch -u _compilers.py -i _compilers.patch > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_scamb "SCRIPT FILE PATCH (_compilers)"
    return 1
  fi

  cd $ROOTDIR/external_modules/code/$CAMB_NAME/fortran
  if [ $? -ne 0 ]; then
    fail_scamb "CD CAMB FORTRAN FOLDER"
    return 1
  fi
  cp $CHANGES/fortran/Makefile.patch .
  if [ $? -ne 0 ]; then
    fail_scamb "CP FILE PATCH (Makefile)"
    return 1
  fi
  patch -u Makefile -i Makefile.patch > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_scamb "SCRIPT FILE PATCH (Makefile)"
    return 1
  fi

  cd $ROOTDIR/external_modules/code/$CAMB_NAME/forutils
  if [ $? -ne 0 ]; then
    fail_scamb "CD CAMB FORUTILS FOLDER"
    return 1
  fi 
  cp $CHANGES/forutils/Makefile_compiler.patch .
  if [ $? -ne 0 ]; then
    fail_scamb "CP FILE PATCH (Makefile_compiler)"
    return 1
  fi
  patch -u Makefile_compiler -i Makefile_compiler.patch > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_scamb "SCRIPT FILE PATCH (Makefile)"
    return 1
  fi

  cd $ROOTDIR
  pbottom 'INSTALLING CAMB'
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  unset_env_vars_scamb
  pbottom2 'SETUP_CAMB'
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------