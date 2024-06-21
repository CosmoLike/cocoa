#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_CAMB_COMPILATION}" ]; then

  pfail() {
    echo -e \
    "\033[0;31m\t\t ERROR ENV VARIABLE ${1:-"empty arg"} NOT DEFINED \033[0m"
    unset pfail
  }
  
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1;
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
  
  unset_env_vars_scamb () {
    unset OUT1
    unset OUT2
    unset URL
    unset pfail
    unset CHANGES
    unset ECODEF
    unset CAMBF
    unset PACKDIR
    unset unset_env_vars_scamb
    cdroot || return 1;
  }

  fail_scb () {
    local MSG="\033[0;31m\t\t (setup_camb.sh) we cannot run \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG}${1:-"empty arg"}${MSG2}"
    unset fail_scb
    unset_env_vars_scamb
  }

  if [ -z "${COCOA_OUTPUT_VERBOSE}" ]; then
    export OUT1="/dev/null"; export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"; export OUT2="/dev/tty"
  fi
  
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { fail_scb "CD FOLDER: ${1}"; return 1; }
  }

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  
  ptop2 'SETUP_CAMB'

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  
  ptop 'INSTALLING CAMB'

  export URL="https://github.com/cmbant/CAMB"

  export CHANGES="${ROOTDIR:?}/../cocoa_installation_libraries/camb_changes"

  export ECODEF="${ROOTDIR:?}/external_modules/code"

  export CAMBF=${CAMB_NAME:-"CAMB"}

  export PACKDIR="${ECODEF:?}/${CAMBF:?}"

  # ---------------------------------------------------------------------------
  # in case this script is called twice
  # ---------------------------------------------------------------------------
  rm -rf "${PACKDIR:?}"

  # ---------------------------------------------------------------------------
  # clone from original repo
  # ---------------------------------------------------------------------------
  cdfolder "${ECODEF:?}" || return 1;

  "${GIT:?}" clone $URL --recursive "${CAMBF:?}" >${OUT1:?} 2>${OUT2:?} || 
    { fail_scb "GIT CLONE FROM CAMB REPO"; return 1; }
  
  cdfolder "${PACKDIR}" || return 1;

  if [ -n "${CAMB_GIT_COMMIT}" ]; then
    "${GIT:?}" checkout "${CAMB_GIT_COMMIT:?}" >${OUT1:?} 2>${OUT2:?} ||
      { fail_scb "GIT CHECKOUT CAMB"; return 1; }
  fi
  
  # ---------------------------------------------------------------------------
  # patch CAMB to be compatible w/ COCOA
  # ---------------------------------------------------------------------------
  
  # ---------------------------------------------------------------------------
  cdfolder "${PACKDIR:?}/camb/" || return 1;

  cp "${CHANGES:?}"/camb/_compilers.patch . 2>${OUT2:?} ||
    { fail_scb "CP FILE PATCH (_compilers)"; return 1; }
  
  patch -u _compilers.py -i _compilers.patch >${OUT1:?} 2>${OUT2:?} ||
    { fail_scb "SCRIPT FILE PATCH (_compilers)"; return 1; }

  # ---------------------------------------------------------------------------
  cdfolder "${PACKDIR:?}/fortran/" || return 1;
  
  cp "${CHANGES:?}"/fortran/Makefile.patch . 2>${OUT2:?} ||
    { fail_scb "CP FILE PATCH (Makefile)"; return 1; }
  
  patch -u Makefile -i Makefile.patch >${OUT1:?} 2>${OUT2:?} ||
    { fail_scb "SCRIPT FILE PATCH (Makefile)"; return 1; }

  # ---------------------------------------------------------------------------
  cdfolder "${PACKDIR:?}/forutils/" || return 1;
  
  cp "${CHANGES:?}/forutils/Makefile_compiler.patch" . 2>${OUT2:?} ||
    { fail_scb "CP FILE PATCH (Makefile_compiler)"; return 1; }
  
  patch -u Makefile_compiler -i Makefile_compiler.patch >${OUT1:?} \
    2>${OUT2:?} || { fail_scb "SCRIPT FILE PATCH (Makefile)"; return 1; }
  
  # ---------------------------------------------------------------------------
  cdfolder "${ROOTDIR}" || return 1;

  pbottom 'INSTALLING CAMB'
  
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  unset_env_vars_scamb || return 1
  
  pbottom2 'SETUP_CAMB'

fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------