#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_HDF5_INSTALLATION}" ]; then
  
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1
  fi

  source "${ROOTDIR:?}/installation_scripts/.check_flags.sh" || return 1;

  unset_env_vars () {
    unset -v CCIL BDF
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
    fail_script_msg "setup_hdf5.sh" "${1}"
    unset_all || return 1
  }
    
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { fail_shdf5 "CD FOLDER: ${1}"; return 1; }
  }
  
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  ptop2 'SETUP_HDF5' || return 1;
  
  unset_env_vars || return 1; 

  CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

  # ----------------------------------------------------------------------------

  ptop 'INSTALLING HFD5 LIBRARY' || return 1;
  
  PACKDIR="${CCIL:?}/${COCOA_HDF5_DIR:-"hdf5-1.12.3/"}"
  BDF="cocoa_HDF5_build"

  rm -f  "${PACKDIR:?}/CMakeCache.txt"
  rm -rf "${PACKDIR:?}/CMakeFiles/"
  rm -rf "${PACKDIR:?}/${BDF:?}/"
  
  mkdir "${PACKDIR:?}/${BDF:?}" 2>${OUT2:?} || 
    { fail_shdf5 "MKDIR COCOA_HDF5_BUILD"; return 1; }

  cdfolder "${PACKDIR:?}/${BDF:?}" || return 1;

  "${CMAKE:?}" -DBUILD_SHARED_LIBS=TRUE \
    -DCMAKE_INSTALL_PREFIX="${ROOTDIR:?}/.local" \
    -DCMAKE_C_COMPILER="${C_COMPILER:?}" \
    -DCMAKE_CXX_COMPILER="${CXX_COMPILER:?}" \
    -DCMAKE_FC_COMPILER="${FORTRAN_COMPILER:?}" \
    --log-level=ERROR .. >${OUT1:?} 2>${OUT2:?} || 
    { fail_shdf5 "CMAKE"; return 1; }

  make -j $MNT >${OUT1:?} 2>${OUT2:?} || { fail_shdf5 "MAKE"; return 1; }

  make install >${OUT1:?} 2>${OUT2:?} || 
    { fail_shdf5 "MAKE INSTALL"; return 1; }

  cdfolder "${ROOTDIR}" || return 1;

  pbottom 'INSTALLING HDF5 LIBRARY DONE' || return 1;

  # ---------------------------------------------------------------------------- 
  
  unset_all || return 1; 

  pbottom2 'SETUP_HDF5'  || return 1;

fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------