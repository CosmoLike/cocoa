#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_CLASS_COMPILATION}" ]; then
    
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1
  fi
  
  # parenthesis = run in a subshell
  ( source "${ROOTDIR:?}/installation_scripts/.check_flags.sh" ) || return 1;
    
  unset_env_vars () {
    unset -v PLIB PACKDIR
    cdroot || return 1;
  }
  
  unset_env_funcs () {
    unset -f cdfolder cpfolder error
    unset -f unset_env_funcs
    cdroot || return 1;
  }

  unset_all () {
    unset_env_vars
    unset_env_funcs
    unset -f unset_all
    cdroot || return 1;
  }

  error () {
    fail_script_msg "clean_class.sh" "${1}"
    unset_all || return 1
  }
  
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER ${1}"; return 1; }
  }
  
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  
  ptop 'CLEANING CLASS' || return 1

  unset_env_vars || return 1

  PLIB="${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages"
  
  PACKDIR="${ROOTDIR:?}/external_modules/code/${CLASS_NAME:-"class_public"}/"

  cdfolder "${PACKDIR}" || return 1

  make clean >${OUT1:?} 2>${OUT2:?} || { error "${EC2:?}"; return 1; }

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

  unset_all || return 1
  
  pbottom 'CLEANING CLASS' || return 1

fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------