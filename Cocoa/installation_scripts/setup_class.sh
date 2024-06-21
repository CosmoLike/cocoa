#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_CLASS_COMPILATION}" ]; then
  
  pfail() {
    echo -e \
    "\033[0;31m\t\t ERROR ENV VARIABLE ${1:-"empty arg"} NOT DEFINED \033[0m"
    unset pfail
  }

  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'
    return 1
  fi
  
  cdroot() {
    cd "${ROOTDIR:?}" 2>"/dev/null" || { echo -e \
      "\033[0;31m\t\t CD ROOTDIR (${ROOTDIR}) FAILED \033[0m"; return 1; }
    unset cdroot
  }
  
  if [ -z "${C_COMPILER}" ]; then
    pfail 'C_COMPILER'; cdroot; return 1
  fi

  if [ -z "${FORTRAN_COMPILER}" ]; then
    pfail 'FORTRAN_COMPILER'; cdroot; return 1
  fi 

  if [ -z "${GIT}" ]; then
    pfail 'GIT'; cdroot; return 1
  fi
    
  unset_env_vars_sclass () {
    unset OUT1
    unset OUT2
    unset URL
    unset CHANGES
    unset pfail
    unset ECODEF
    unset PACKDIR
    unset CLNAME
    unset unset_env_vars_sclass
    cdroot || return 1;
  }
  
  fail_sclass () {
    local MSG="\033[0;31m (setup_class.sh) WE CANNOT RUN \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG} ${1:-"empty arg"} ${MSG2}"
    unset fail_sclass
    unset_env_vars_sclass
  }

  if [ -z "${COCOA_OUTPUT_VERBOSE}" ]; then
    export OUT1="/dev/null"; export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"; export OUT2="/dev/tty"
  fi

  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { fail_sclass "CD FOLDER: ${1}"; return 1; }
  }

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  
  ptop2 'SETUP_CLASS'
  
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------

  ptop 'INSTALLING CLASS'

  export URL="https://github.com/lesgourg/class_public.git"
  
  export CHANGES="${ROOTDIR:?}/../cocoa_installation_libraries/class_changes"

  export ECODEF="${ROOTDIR:?}/external_modules/code"

  export CLNAME=${CLASS_NAME:-"class_public"}
  
  export PACKDIR="${ECODEF:?}/${CLNAME:?}/"

  # ---------------------------------------------------------------------------
  # in case this script is called twice
  # ---------------------------------------------------------------------------
  rm -rf "${PACKDIR:?}"

  # ---------------------------------------------------------------------------
  # clone from original repo
  # ---------------------------------------------------------------------------
  cdfolder "${ECODEF}" || return 1;

  "${GIT:?}" clone "${URL:?}" --recursive "${CLNAME:?}" >${OUT1:?} \
    2>${OUT2:?} || { fail_sclass "GIT CLONE FROM CLASS REPO"; return 1; }
  
  cdfolder "${PACKDIR}" || return 1;

  if [ -n "${CLASS_GIT_COMMIT}" ]; then
    "${GIT:?}" checkout "${CLASS_GIT_COMMIT:?}" >${OUT1:?} 2>${OUT2:?} ||
      { fail_sclass "GIT CHECKOUT CLASS"; return 1; }
  fi
  
  # ---------------------------------------------------------------------------
  # historical reasons (we used to save class_python on Cocoa Branch)
  # historical: Workaround Cocoa .gitignore entry on /include
  # ---------------------------------------------------------------------------
  mv ./include ./include2/ >${OUT1:?} 2>${OUT2:?} || 
    { fail_sclass "MKDIR FROM INCLUDE TO INCLUDE2"; return 1; }
  
  # ---------------------------------------------------------------------------
  # patch CLASS to be compatible w/ COCOA
  # ---------------------------------------------------------------------------
  
  # ---------------------------------------------------------------------------
  cdfolder "${PACKDIR}" || return 1;

  cp "${CHANGES:?}/Makefile.patch" .  2>${OUT2:?} || 
    { fail_sclass "CP FILE PATCH (Makefile.patch)"; return 1; }
  
  patch -u Makefile -i Makefile.patch >${OUT1:?} 2>${OUT2:?} ||
    { fail_sclass "SCRIPT FILE PATCH (Makefile.patch)"; return 1; }
  
  # ---------------------------------------------------------------------------
  cdfolder "${PACKDIR:?}/python" || return 1;
  
  cp "${CHANGES:?}/python/setup.patch" . 2>${OUT2:?} ||
    { fail_sclass "CP FILE PATCH (setup.patch)"; return 1; }

  patch -u setup.py -i setup.patch >${OUT1:?} 2>${OUT2:?} || 
    { fail_sclass "SCRIPT FILE PATCH (setup.patch)"; return 1; }
  
  # ---------------------------------------------------------------------------
  cdfolder "${ROOTDIR}" || return 1;

  pbottom 'INSTALLING CLASS'

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  unset_env_vars_sclass || return 1

  pbottom2 'SETUP_CLASS'
fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
