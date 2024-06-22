#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_CLASS_COMPILATION}" ]; then
  
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1
  fi
  
  # ---------------------------------------------------------------------------
  # Clean any previous compilation
  source "${ROOTDIR:?}/installation_scripts/clean_class.sh" || return 1
  # ---------------------------------------------------------------------------
  
  source "${ROOTDIR:?}/installation_scripts/.check_flags.sh" || return 1;
      
  unset_env_vars () {
    unset -v PACKDIR
    cdroot || return 1;
  }
  
  error () {
    fail_script_msg "compile_class.sh" "${1}"
    unset -f error
    unset_env_vars || return 1
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
  
  ptop 'COMPILING CLASS' || return 1

  unset_env_vars || return 1

  PACKDIR="${ROOTDIR:?}/external_modules/code/${CLASS_NAME:-"class_public"}"

  cdfolder "${PACKDIR}" || return 1

  # ---------------------------------------------------------------------------
  # (motivation): Cocoa has /include entry on .gitignore. Why? Class compilation
  # (motivation): changes files on /include. Maintaining original files on
  # (motivation): /include2 avoids GIT changes to be triggered on compilation.
  # Delete folder in case this script is called twice
  # ---------------------------------------------------------------------------  
  rm -rf "${PACKDIR:?}/include"

  cpfolder "${PACKDIR:?}/include2" "${PACKDIR:?}/include" || return 1;

  # --------------------------------------------------------------------------- 
  CC="${C_COMPILER:?}" PYTHON=${PYTHON3:?} make all \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC7:?}"; return 1; }
   
  cdfolder "${PACKDIR}/python" || return 1

  CC="${C_COMPILER:?}" ${PYTHON3:?} setup.py build \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC4:?}"; return 1; }

  unset_env_vars || return 1

  pbottom 'COMPILING CLASS' || return 1

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------