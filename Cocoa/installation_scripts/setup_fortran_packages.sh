#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_FORTRAN_INSTALLATION}" ]; then

  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1
  fi

  source "${ROOTDIR:?}/installation_scripts/.check_flags.sh" || return 1;
      
  unset_env_vars () {
    unset -v CCIL PACKDIR PACKDIR
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
    fail_script_msg "setup_fortran_packages.sh" "${1}"
    unset_all || return 1
  }

  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
  }

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  ptop2 'SETUP_FORTRAN_PACKAGES' || return 1;

  unset_env_vars || return 1; 

  CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

  # ----------------------------------------------------------------------------
  # ------------------------------- FORTRAN LAPACK -----------------------------
  # ----------------------------------------------------------------------------
  if [ -z "${IGNORE_FORTRAN_LAPACK_INSTALLATION}" ]; then
    
    ptop 'INSTALLING LAPACK FORTRAN LIBRARY' || return 1;

    PACKDIR="${COCOA_LAPACK_DIR:-"lapack-3.11.0/"}"
    BUILDIR="lapack-build"

    rm -rf "${CCIL:?}/${BUILDIR:?}"
    
    mkdir "${CCIL:?}/${BUILDIR:?}" 2>${OUT2:?} || 
      { error "(LAPACK) MKDIR LAPACK-BUILD"; return 1; }

    cdfolder "${CCIL:?}/${BUILDIR:?}" || return 1;

    ${CMAKE:?} -DBUILD_SHARED_LIBS=TRUE \
      -DCMAKE_INSTALL_PREFIX="${ROOTDIR:?}/.local" \
      -DCMAKE_C_COMPILER="${C_COMPILER:?}" \
      --log-level=ERROR "../${PACKDIR:?}" >${OUT1:?} 2>${OUT2:?} ||
      { error "(LAPACK) CMAKE"; return 1; }

    make -j $MNT all >${OUT1:?} 2>${OUT2:?} || { error "MAKE"; return 1; }

    make install >${OUT1:?} 2>${OUT2:?} || { error "MAKE INSTALL"; return 1; }
    
    cdfolder "${ROOTDIR}" || return 1;

    pbottom 'INSTALLING LAPACK FORTRAN LIBRARY' || return 1;

  fi

  # ---------------------------------------------------------------------------- 
  
  unset_all || return 1; 
  
  pbottom2 'SETUP_FORTRAN_PACKAGES' || return 1;

fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------