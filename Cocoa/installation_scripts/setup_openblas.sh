#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_OPENBLAS_INSTALLATION}" ]; then
  
  pfail() {
    echo -e "\033[0;31m\t\t ERROR ENV VARIABLE ${1} IS NOT DEFINED \033[0m"
    unset pfail
  }
  
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1
  fi
  
  cdroot() {
    cd "${ROOTDIR}" 2>"/dev/null" || { echo -e \
      "\033[0;31m\t\t CD ROOTDIR (${ROOTDIR}) FAILED \033[0m"; return 1; }
    unset cdroot
  }
  
  if [ -z "${CXX_COMPILER}" ]; then
    pfail 'CXX_COMPILER'; cdroot; return 1;
  fi
  
  if [ -z "${C_COMPILER}" ]; then
    pfail 'C_COMPILER'; cdroot; return 1;
  fi
  
  if [ -z "${FORTRAN_COMPILER}" ]; then
    pfail 'FORTRAN_COMPILER'; cdroot; return 1;
  fi
  
  unset_env_vars_sopb () {
    unset MAKE_NB_JOBS
    unset OUT1
    unset OUT2
    unset OPBMNT
    unset unset_env_vars_sopb
    cdroot || return 1;
  }
  
  fail_sopb () {
    local MSG="\033[0;31m\t\t (setup_openblas.sh) WE CANNOT RUN \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"
    unset fail_sopb
    unset_env_vars_sftrp
  }
  
  if [ -z "${DEBUG_OPENBLAS_PACKAGE}" ]; then
    if [ -z "${MAKE_NUM_THREADS}" ]; then
      pfail 'MAKE_NUM_THREADS'; cdroot; return 1;
    fi
    export OUT1="/dev/null"; export OUT2="/dev/null"
    export OPBMNT="${MAKE_NUM_THREADS}"
  else
    export OUT1="/dev/tty"; export OUT2="/dev/tty"
    export OPBMNT=1
  fi

  cdfolder() {
    cd "${1}" 2>"/dev/null" || { fail_sopb "CD FOLDER: ${1}"; return 1; }
  }
  
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  ptop2 'SETUP_OPENBLAS'
  
  export CCIL="${ROOTDIR}/../cocoa_installation_libraries"

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  ptop  'INSTALLING OPENBLAS LIBRARY'
  
  if [ -z "${COCOA_OPENBLAS_DIR}" ]; then
    pfail 'COCOA_OPENBLAS_DIR'; cdroot; return 1;
  fi

  cdfolder "${CCIL}/${COCOA_OPENBLAS_DIR}" || return 1;

  export MAKE_NB_JOBS=$OPBMNT
  
  make clean >${OUT1} 2>${OUT2} || { fail_sopb "MAKE CLEAN"; return 1; }

  make CC=$C_COMPILER FC=$FORTRAN_COMPILER USE_OPENMP=1 >${OUT1} 2>${OUT2} || 
    { fail_sopb "MAKE"; return 1; }
  
  make install PREFIX="${ROOTDIR}/.local" >${OUT1} 2>${OUT2} || 
    { fail_sopb "MAKE INSTALL"; return 1; }

  cdfolder "${ROOTDIR}" || return 1;

  pbottom  'INSTALLING OPENBLAS LIBRARY'
  
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  unset_env_vars_sopb || return 1;
  
  pbottom2 'SETUP_OPENBLAS'

fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------