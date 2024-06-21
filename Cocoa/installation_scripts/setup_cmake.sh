#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_CMAKE_INSTALLATION}" ]; then
  
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
    
  if [ -z "${CPP_COMPILER}" ]; then
    pfail 'CPP_COMPILER'; cdroot; return 1
  fi

  if [ -z "${C_COMPILER}" ]; then
    pfail 'C_COMPILER'; cdroot; return 1
  fi
    
  unset_env_vars_scmk () {
    unset OUT1
    unset OUT2
    unset CMKMNT
    unset pfail
    unset unset_env_vars_scmk
    cdroot || return 1;
  }
  
  fail_scmk () {
    local MSG="\033[0;31m (setup_cmake.sh) we cannot run \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG} ${1:-"empty arg"} ${MSG2}"
    unset fail_scmk
    unset_env_vars_scmk
  }
  
  if [ -z "${DEBUG_CMAKE_PACKAGE}" ]; then
    export OUT1="/dev/null"; export OUT2="/dev/null"
    export CMKMNT="${MAKE_NUM_THREADS:-1}"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
    export CMKMNT=1
  fi

  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { fail_scmk "CD FOLDER: ${1}"; return 1; }
  }

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  ptop2 'SETUP_CMAKE'

  export CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  ptop 'INSTALLING CMAKE LIBRARY'

  if [ -z "${COCOA_CMAKE_DIR}" ]; then
    pfail 'COCOA_CMAKE_DIR'; cdroot; return 1;
  fi

  cdfolder "${CCIL:?}/${COCOA_CMAKE_DIR:?}" || return 1;

  env CC="${C_COMPILER:?}" CXX="${CPP_COMPILER:?}" ./bootstrap \
    --prefix="${ROOTDIR:?}/.local" >${OUT1:?} 2>${OUT2:?} || 
    { fail_scmk "BOOTSTRAP"; return 1; }

  make -j $CMKMNT >${OUT1:?} 2>${OUT2:?} || { fail_scmk "MAKE"; return 1; }

  make install >${OUT1:?} 2>${OUT2:?} || { fail_scmk "INSTALL"; return 1 }

  cdfolder "${ROOTDIR}" || return 1;

  pbottom 'INSTALLING CMAKE LIBRARY'
  
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  unset_env_vars_scmk
  pbottom2 'SETUP_CMAKE'
fi
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------