#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_CLASS_COMPILATION}" ]; then
  
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
  
  if [ -z "${PYTHON3}" ]; then
    pfail "PYTHON3"; cdroot; return 1
  fi
  
  if [ -z "${PYTHON_VERSION}" ]; then
    pfail "PYTHON_VERSION"; cdroot; return 1
  fi
    
  unset_env_vars_clean_class () {
    unset OUT1
    unset OUT2
    unset pfail
    unset PLIB
    unset PACKDIR
    unset unset_env_vars_clean_class
    cdroot || return 1;
  }
  
  fail_clcls () {
    local MSG="\033[0;31m\t\t (clean_class.sh) we cannot run \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG} ${1:-"empty arg"} ${MSG2}"
    unset fail_clcls
    unset_env_vars_clean_class
  }
  
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { fail_clcls "CD FOLDER: ${1}"; return 1; }
  }
  
  if [ -z "${COCOA_OUTPUT_VERBOSE}" ]; then
    export OUT1="/dev/null"; export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"; export OUT2="/dev/tty"
  fi

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------

  ptop 'CLEANING CLASS'

  export PLIB="${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages"
  export PACKDIR="${ROOTDIR:?}/external_modules/code/${CLASS_NAME:-"class_public"}/"

  cdfolder "${PACKDIR}" || return 1

  make clean >${OUT1:?} 2>${OUT2:?} || { fail_clcls "MAKE CLEAN"; return 1; }

  cdfolder "${PACKDIR:?}/python"|| return 1
  
  # ---------------------------------------------------------------------------
  # below we ignore if something goes wrong (related to include/ relocation)
  "${PYTHON3:?}" setup.py clean >${OUT1:?} 2>${OUT2:?}

  rm -rf "${PLIB:?}"/classy*
  rm -rf "${PACKDIR:?}/python/build/"
  rm -rf "${PACKDIR:?}/python/classy.egg-info"  
  rm -rf "${PACKDIR:?}/build/"
  rm -f "${PACKDIR:?}/class"
  rm -f "${PACKDIR:?}/libclass.a"
  # ---------------------------------------------------------------------------
  # Historical Workaround Cocoa .gitignore entry on /include
  # ---------------------------------------------------------------------------
  rm -rf "${PACKDIR:?}/include"

  unset_env_vars_clean_class || return 1
  
  pbottom 'CLEANING CLASS'
  
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------