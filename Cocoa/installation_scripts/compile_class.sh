#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_CLASS_COMPILATION}" ]; then
  
  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi
    
  # parenthesis = run in a subshell 
  ( source "${ROOTDIR:?}/installation_scripts/.check_flags.sh" ) || return 1;
      
  unset_env_vars () {
    unset -v ECODEF FOLDER PACKDIR PLIB
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
    fail_script_msg "$(basename ${BASH_SOURCE[0]})" "${1}"
    unset_all || return 1
  }
  
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER ${1}"; return 1; }
  }

  cpfolder() {
    cp -r "${1:?}" "${2:?}" \
      2>"/dev/null" || { error "CP FOLDER ${1} on ${2}"; return 1; }
  }

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  
  ptop 'COMPILING CLASS' || return 1

  unset_env_vars || return 1

  # E = EXTERNAL, CODE, F=FODLER
  ECODEF="${ROOTDIR:?}/external_modules/code"

  FOLDER=${CLASS_NAME:-"class_public"}

  PACKDIR="${ECODEF:?}/${FOLDER:?}"

  cdfolder "${PACKDIR}" || return 1

  # ---------------------------------------------------------------------------
  # cleaning any previous compilation
  
  # note: below we ignore if something goes wrong (related to include/ mv)
  "${PYTHON3:?}" setup.py clean >${OUT1:?} 2>${OUT2:?}

  PLIB="${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages"
  rm -rf "${PLIB:?}"/classy*

  rm -rf "${PACKDIR:?}/python/build/"
  rm -rf "${PACKDIR:?}/python/classy.egg-info"  
  rm -rf "${PACKDIR:?}/build/"
  rm -f "${PACKDIR:?}/class"
  rm -f "${PACKDIR:?}/libclass.a"
  
  # ---------------------------------------------------------------------------
  # note: historical motivation when class was inside Cocoa git repo
  # (motivation): Cocoa has /include entry on .gitignore. Why? Class compilation
  # (motivation): changes files on /include. Maintaining original files on
  # (motivation): /include2 avoids GIT changes to be triggered on compilation.
  # Delete folder in case this script is called twice
  # --------------------------------------------------------------------------- 
  rm -rf "${PACKDIR:?}/include"

  # ---------------------------------------------------------------------------
  
  cpfolder "${PACKDIR:?}/include2" "${PACKDIR:?}/include" || return 1;

  CC="${C_COMPILER:?}" PYTHON=${PYTHON3:?} make all \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC7:?}"; return 1; }
   
  cdfolder "${PACKDIR}/python" || return 1

  CC="${C_COMPILER:?}" ${PYTHON3:?} setup.py build \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC4:?}"; return 1; }

  unset_all || return 1

  pbottom 'COMPILING CLASS' || return 1

fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------