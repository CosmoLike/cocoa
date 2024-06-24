#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_OPENBLAS_INSTALLATION}" ]; then

  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1
  fi

  source "${ROOTDIR:?}/installation_scripts/.check_flags.sh" || return 1;

  
  unset_env_vars () {
    unset -v MAKE_NB_JOBS PACKDIR
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
    fail_script_msg "setup_openblas.sh" "${1}"
    unset_all || return 1
  }
  
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
  }
  
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  ptop2 'SETUP_OPENBLAS' || return 1;

  unset_env_vars || return 1; 
  
  CCIL="${ROOTDIR:?}/../cocoa_installation_libraries" || return 1;

  # ----------------------------------------------------------------------------  
  
  ptop  'INSTALLING OPENBLAS LIBRARY' || return 1;
  
  export PACKDIR="${COCOA_OPENBLAS_DIR:-"OpenBLAS-0.3.23/"}"

  cdfolder "${CCIL:?}/${PACKDIR:?}" || return 1;

  export MAKE_NB_JOBS=$MNT
  
  make clean >${OUT1:?} 2>${OUT2:?} || { error "${EC2:?}"; return 1; }

  make CC="${C_COMPILER:?}" FC="${FORTRAN_COMPILER:?}" USE_OPENMP=1 \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC8:?}"; return 1; }
  
  make install PREFIX="${ROOTDIR:?}/.local" \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC10:?}"; return 1; }

  cdfolder "${ROOTDIR:?}" || return 1;

  pbottom  'INSTALLING OPENBLAS LIBRARY' || return 1;
  
  # ----------------------------------------------------------------------------
  
  unset_all || return 1;
  
  pbottom2 'SETUP_OPENBLAS' || return 1;

fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------