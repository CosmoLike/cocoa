#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_FORTRAN_INSTALLATION}" ]; then
  
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
  
  if [ -z "${CMAKE}" ]; then
    pfail 'CMAKE'; cdroot; return 1;
  fi
  
  unset_env_vars_sftrp () {
    unset OUT1
    unset OUT2
    unset FMNT
    unset CCIL
    unset pfail
    unset unset_env_vars_sftrp
    cdroot || return 1;
  }
  
  fail_sftrp () {
    local MSG="\033[0;31m\t\t (setup_fortran_packages.sh) WE CANNOT RUN \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"
    unset fail_sftrp
    unset_env_vars_sftrp
  }

  if [ -z "${DEBUG_FORTRAN_PACKAGES}" ]; then
    if [ -z "${MAKE_NUM_THREADS}" ]; then
      pfail 'MAKE_NUM_THREADS'; cdroot; return 1;
    fi
    export OUT1="/dev/null"; export OUT2="/dev/null"
    export FMNT="${MAKE_NUM_THREADS}"
  else
    export OUT1="/dev/tty"; export OUT2="/dev/tty"
    export FMNT=1
  fi

  cdfolder() {
    cd "${1}" 2>"/dev/null" || { fail_sftrp "CD FOLDER: ${1}"; return 1; }
  }

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  ptop2 'SETUP_FORTRAN_PACKAGES'

  export CCIL="${ROOTDIR}/../cocoa_installation_libraries"

  # ----------------------------------------------------------------------------
  # ------------------------------- FORTRAN LAPACK -----------------------------
  # ----------------------------------------------------------------------------
  if [ -z "${IGNORE_FORTRAN_LAPACK_INSTALLATION}" ]; then
    
    ptop 'INSTALLING LAPACK FORTRAN LIBRARY'

    if [ -z "${COCOA_LAPACK_DIR}" ]; then
      pfail 'COCOA_LAPACK_DIR'; cdroot; return 1;
    fi

    cdfolder "${CCIL}" || return 1;

    rm -rf ./lapack-build
    
    mkdir lapack-build >${OUT1} 2>${OUT2} || 
      { fail_sftrp "(LAPACK) MKDIR LAPACK-BUILD"; return 1; }

    cdfolder "${CCIL}/lapack-build" || return 1;

    $CMAKE -DBUILD_SHARED_LIBS=TRUE \
      -DCMAKE_INSTALL_PREFIX="${ROOTDIR}/.local" \
      -DCMAKE_C_COMPILER=$C_COMPILER \
      --log-level=ERROR ../"${COCOA_LAPACK_DIR}" \
      >${OUT1} 2>${OUT2} || { fail_sftrp "(LAPACK) CMAKE"; return 1; }

    make -j $FMNT all >${OUT1} 2>${OUT2} || { fail_sftrp "MAKE"; return 1; }

    make install >${OUT1} 2>${OUT2} || { fail_sftrp "MAKE INSTALL"; return 1; }
    
    cdfolder "${ROOTDIR}" || return 1;

    pbottom 'INSTALLING LAPACK FORTRAN LIBRARY'

  fi
  # ----------------------------------------------------------------------------
  # ---------------------------------------------------------------------------- 
  
  unset_env_vars_sftrp || return 1; 
  
  pbottom2 'SETUP_FORTRAN_PACKAGES'
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------