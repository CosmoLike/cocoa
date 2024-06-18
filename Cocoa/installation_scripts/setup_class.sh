#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_CLASS_COMPILATION}" ]; then
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
  if [ -z "${CLASS_NAME}" ]; then
    pfail 'CLASS_NAME'
    cd $ROOTDIR
    return 1
  fi
  unset_env_vars_sclass () {
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset URL
    unset CHANGES
    unset pfail
    unset unset_env_vars_sclass
  }
  fail_sclass () {
    export MSG="\033[0;31m (setup_class.sh) WE CANNOT RUN \e[3m"
    export MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"
    unset_env_vars_sclass
    unset MSG
    unset MSG2
    unset fail_sclass
  }
  if [ -z "${DEBUG_CAMB_OUTPUT}" ]; then
    export OUT1="/dev/null"
    export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
  fi

  # ---------------------------------------------------------------------------
  ptop2 'SETUP_CLASS'
  ptop 'INSTALLING CLASS'

  export URL="https://github.com/lesgourg/class_public.git"
  export CHANGES="${ROOTDIR}/../cocoa_installation_libraries/class_changes"

  # ---------------------------------------------------------------------------
  # in case this script is called twice
  # ---------------------------------------------------------------------------
  rm -rf $ROOTDIR/external_modules/code/$CLASS_NAME

  # ---------------------------------------------------------------------------
  # clone from original repo
  # ---------------------------------------------------------------------------
  cd $ROOTDIR/external_modules/code/

  $GIT clone $URL --recursive $CLASS_NAME > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_sclass "GIT CLONE FROM CLASS REPO"
    return 1
  fi

  cd $ROOTDIR/external_modules/code/$CLASS_NAME

  $GIT checkout $CLASS_GIT_COMMIT > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_sclass "GIT CHECKOUT CLASS"
    return 1
  fi

  # ---------------------------------------------------------------------------
  # historical reasons (we used to save class_python on Cocoa Branch)
  # historical: Workaround Cocoa .gitignore entry on /include
  # ---------------------------------------------------------------------------
  mv ./include ./include2/
  
  # ---------------------------------------------------------------------------
  # patch CLASS to be compatible w/ COCOA
  # ---------------------------------------------------------------------------
  cd $ROOTDIR/external_modules/code/$CLASS_NAME/
  if [ $? -ne 0 ]; then
    fail_sclass "CD CLASS"
    return 1
  fi
  cp $CHANGES/Makefile.patch .
  if [ $? -ne 0 ]; then
    fail_sclass "CP FILE PATCH (Makefile)"
    return 1
  fi
  patch -u Makefile -i Makefile.patch > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_sclass "SCRIPT FILE PATCH (Makefile)"
    return 1
  fi

  cd $ROOTDIR/external_modules/code/$CLASS_NAME/python
  if [ $? -ne 0 ]; then
    fail_sclass "CD CLASS PYTHON"
    return 1
  fi
  cp $CHANGES/python/setup.patch .
  if [ $? -ne 0 ]; then
    fail_sclass "CP FILE PATCH (setup)"
    return 1
  fi
  patch -u setup.py -i setup.patch > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_sclass "SCRIPT FILE PATCH (setup)"
    return 1
  fi

  cd $ROOTDIR
  pbottom 'INSTALLING CLASS'
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  unset_env_vars_sclass
  pbottom2 'SETUP_CLASS'
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
