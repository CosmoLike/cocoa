#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_DISTUTILS_INSTALLATION}" ]; then
  
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
  
  if [ -z "${CXX_COMPILER}" ]; then
    pfail 'CXX_COMPILER'; cdroot; return 1;
  fi
  
  if [ -z "${C_COMPILER}" ]; then
    pfail 'C_COMPILER'; cdroot; return 1;
  fi
  
  if [ -z "${FORTRAN_COMPILER}" ]; then
    pfail 'FORTRAN_COMPILER'; cdroot; return 1;
  fi 
    
  unset_env_vars_sbinutils () {
    unset OUT1
    unset OUT2
    unset SBMNT
    unset pfail
    unset CCIL
    unset PACKDIR
    unset unset_env_vars_sbinutils
    cdroot || return 1;
  }
  
  fail_sbu () {
    local MSG="\033[0;31m\t\t (setup_binutils.sh) WE CANNOT RUN \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG} ${1:-"empty arg"} ${MSG2}"
    unset fail_sbu
    unset_env_vars_sbinutils
  }

  if [ -z "${COCOA_OUTPUT_VERBOSE}" ]; then
    export OUT1="/dev/null"; export OUT2="/dev/null"
    export SBMNT="${MAKE_NUM_THREADS:-1}"
    [[ ${SBMNT} == +([0-9]) ]] || export SBMNT=1
  else
    export OUT1="/dev/tty"; export OUT2="/dev/tty"
    export SBMNT=1
  fi

  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { fail_sbu "CD FOLDER: ${1}"; return 1; }
  }

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  ptop2 'SETUP_BINUTILS'
  
  export CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

  # ----------------------------------------------------------------------------
  # ----------------------------- TEXINFO LIBRARY  -----------------------------
  # ----------------------------------------------------------------------------
  
  ptop 'INSTALLING TEXINFO LIBRARY'

  export PACKDIR="${CCIL:?}/${COCOA_TEXINFO_DIR:-"texinfo-7.0.3/"}"

  cdfolder "${PACKDIR}" || return 1;

  FC="${FORTRAN_COMPILER:?}" CC="${C_COMPILER:?}" ./configure \
    --prefix="${ROOTDIR:?}/.local" --disable-perl-xs \
    >${OUT1:?} 2>${OUT2:?} || { fail_sbu "CONFIGURE"; return 1; }

  make -j $SBMNT all >${OUT1:?} 2>${OUT2:?} || { fail_sbu "MAKE ALL"; return 1; }
    
  make install >${OUT1:?} 2>${OUT2:?} || { fail_sbu "MAKE INSTALL"; return 1; }

  cdfolder "${ROOTDIR}" || return 1;
  
  pbottom 'INSTALLING TEXINFO LIBRARY'
  
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  # ----------------------------------------------------------------------------
  # ----------------------------- DISTUTILS LIBRARY  ---------------------------
  # ----------------------------------------------------------------------------
  
  ptop 'INSTALLING BINUTILS LIBRARY'
  
  export PACKDIR="${CCIL:?}/${COCOA_BINUTILS_DIR:-"binutils-2.37/"}"

  cdfolder "${PACKDIR}" || return 1;

  FC="${FORTRAN_COMPILER:?}" CC="${C_COMPILER:?}" ./configure \
    --prefix="${ROOTDIR:?}/.local" >${OUT1:?} 2>${OUT2:?} || 
    { fail_sbu "CONFIGURE"; return 1; }

  make -j $SBMNT >${OUT1:?} 2>${OUT2:?} || { fail_sbu "MAKE"; return 1; }

  make install >${OUT1:?} 2>${OUT2:?} || { fail_sbu "MAKE INSTALL"; return 1; }

  cdfolder "${ROOTDIR}" || return 1;

  pbottom 'INSTALLING BINUTILS LIBRARY'
  
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  unset_env_vars_sbinutils || return 1;
  
  pbottom2 'SETUP_BINUTILS'
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------