#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_CMAKE_INSTALLATION}" ]; then
  
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1;
  fi

  source "${ROOTDIR:?}/installation_scripts/.check_flags.sh" || return 1;
        
  unset_env_vars () {
    unset -v PACKDIR
    cdroot || return 1;
  }
  
  error () {
    fail_script_msg "setup_cmake.sh" "${1}"
    unset -f error
    unset_env_vars || return 1
  }
  
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
  }

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  ptop2 'SETUP_CMAKE'|| return 1

  unset_env_vars || return 1; 

  CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  ptop 'INSTALLING CMAKE LIBRARY' || return 1

  PACKDIR=${COCOA_CMAKE_DIR:-"cmake-3.26.4/"}

  cdfolder "${CCIL:?}/${PACKDIR:?}" || return 1;

  env CC="${C_COMPILER:?}" CXX="${CXX_COMPILER:?}" \
    ./bootstrap --prefix="${ROOTDIR:?}/.local" \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC19:?}"; return 1; }

  make -j $MNT >${OUT1:?} 2>${OUT2:?} || { error "${EC8:?}"; return 1; }

  make install >${OUT1:?} 2>${OUT2:?} || { error "${EC10:?}"; return 1; }

  cdfolder "${ROOTDIR}" || return 1;

  pbottom 'INSTALLING CMAKE LIBRARY' || return 1
  
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  unset_env_vars || return 1

  pbottom2 'SETUP_CMAKE' || return 1

fi

# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------