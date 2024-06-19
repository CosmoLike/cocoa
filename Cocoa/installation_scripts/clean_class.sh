#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_CLASS_COMPILATION}" ]; then
  pfail() {
    echo -e "\033[0;31m ERROR ENV VARIABLE ${1} NOT DEFINED \033[0m"
    unset pfail
  }
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'
    return 1
  fi
  cdroot() {
    cd "${ROOTDIR}" 2>"/dev/null" || { echo -e \
      "\033[0;31m\t\t CD ROOTDIR (${ROOTDIR}) FAILED \033[0m"; return 1; }
    unset cdroot
  }
  if [ -z "${PYTHON3}" ]; then
    pfail "PYTHON3"
    cdroot
    return 1
  fi
  if [ -z "${CLASS_NAME}" ]; then
    pfail 'CAMB_NAME'
    cdroot
    return 1
  fi
  unset_env_vars_clean_class () {
    unset OUT1
    unset OUT2
    unset PLIB
    unset pfail
    unset unset_env_vars_clean_class
    cdroot
  }
  fail_clean_class () {
    local MSG="\033[0;31m\t\t (clean_class.sh) WE CANNOT RUN \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"
    unset fail_clean_class
    unset_env_vars_clean_class
  }
  cdfolder() {
    cd "${1}" 2>"/dev/null" || { fail_clean_class "CD FOLDER: ${1}"; return 1; }
  }
  if [ -z "${DEBUG_CLASS_OUTPUT}" ]; then
    export OUT1="/dev/null"
    export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
  fi

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  ptop 'CLEANING CLASS'

  cdfolder "${ROOTDIR}/external_modules/code/${CLASS_NAME}"/ || return 1

  export PLIB="${ROOTDIR}/.local/lib/python${PYTHON_VERSION}/site-packages"
  rm -rf "${PLIB}"/classy*

  make clean >${OUT1} 2>${OUT2} || { fail_clean_class "MAKE CLEAN"; return 1; }

  rm -f class

  cdfolder ./python || return 1

  $PYTHON3 setup.py clean >${OUT1} 2>${OUT2} ||
    { fail_clean_class "PYTHON SETUP CLEAN"; return 1; }

  # ---------------------------------------------------------------------------
  # Historical Workaround Cocoa .gitignore entry on /include
  # ---------------------------------------------------------------------------
  rm -rf ./include 
  # ---------------------------------------------------------------------------
  rm -rf ./build/
  rm -rf ./classy.egg-info

  unset_env_vars_clean_class || return 1
  pbottom 'CLEANING CLASS'
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------