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
  
  if [ -z "${CLASS_NAME}" ]; then
    pfail 'CAMB_NAME'; cdroot; return 1
  fi
  
  unset_env_vars_clean_class () {
    unset OUT1
    unset OUT2
    unset PLIB
    unset pfail
    unset CLASSDIR
    unset unset_env_vars_clean_class
    cdroot
  }
  
  fail_cl_class () {
    local MSG="\033[0;31m\t\t (clean_class.sh) we cannot run \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG} ${1:-"empty arg"} ${MSG2}"
    unset fail_cl_class
    unset_env_vars_clean_class
  }
  
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { fail_cl_class "CD FOLDER: ${1}"; return 1; }
  }
  
  if [ -z "${DEBUG_CLASS_OUTPUT}" ]; then
    export OUT1="/dev/null"; export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"; export OUT2="/dev/tty"
  fi

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------

  ptop 'CLEANING CLASS'

  extern CLASSDIR="${ROOTDIR:?}/external_modules/code/${CLASS_NAME:?}/"

  cdfolder "${CLASSDIR:?}" || return 1

  make clean >${OUT1:?} 2>${OUT2:?} || { fail_cl_class "MAKE CLEAN"; return 1; }

  cdfolder "${CLASSDIR:?}/python"|| return 1
  
  "${PYTHON3:?}" setup.py clean >${OUT1:?} 2>${OUT2:?} ||
    { fail_cl_class "PYTHON SETUP CLEAN"; return 1; }

  export PLIB="${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages"
  rm -rf "${PLIB:?}"/classy*
  rm -rf "${CLASSDIR:?}/python/build/"
  rm -rf "${CLASSDIR:?}/python/classy.egg-info"  
  rm -rf "${CLASSDIR:?}/build/"
  rm -f "${CLASSDIR:?}/class"
  rm -f "${CLASSDIR:?}/libclass.a"
  # ---------------------------------------------------------------------------
  # Historical Workaround Cocoa .gitignore entry on /include
  # ---------------------------------------------------------------------------
  rm -rf "${CLASSDIR:?}/include"

  unset_env_vars_clean_class || return 1
  
  pbottom 'CLEANING CLASS'
  
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------