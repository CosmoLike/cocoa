#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_HDF5_INSTALLATION}" ]; then
  
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
  
  if [ -z "${CMAKE}" ]; then
    pfail 'CMAKE'; cdroot; return 1;
  fi
  
  unset_env_vars_shdf5 () {
    unset OUT1
    unset OUT2
    unset HDF5MNT
    unset CCIL
    unset BDF
    unset pfail
    unset unset_env_vars_shdf5
    cdroot || return 1;
  }

  fail_shdf5 () {
    local MSG="\033[0;31m\t\t (setup_hdf5.sh) WE CANNOT RUN \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG}${1:-"empty arg"}${MSG2}"
    unset fail_shdf5
    unset_env_vars_shdf5
  }
  
  if [ -z "${COCOA_OUTPUT_VERBOSE}" ]; then
    export OUT1="/dev/null"; export OUT2="/dev/null"
    export HDF5MNT="${MAKE_NUM_THREADS:-1}"
    [[ ${HDF5MNT} == +([0-9]) ]] || export HDF5MNT=1
  else
    export OUT1="/dev/tty"; export OUT2="/dev/tty"
    export HDF5MNT=1
  fi
  
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { fail_shdf5 "CD FOLDER: ${1}"; return 1; }
  }
  
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  ptop2 'SETUP_HDF5'
  
  export CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  ptop 'INSTALLING HFD5 LIBRARY'
  
  export PACKDIR="${CCIL:?}/${COCOA_HDF5_DIR:-"hdf5-1.12.3/"}"
  export BDF="cocoa_HDF5_build"

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

  make -j $HDF5MNT >${OUT1:?} 2>${OUT2:?} || { fail_shdf5 "MAKE"; return 1; }

  make install >${OUT1:?} 2>${OUT2:?} || 
    { fail_shdf5 "MAKE INSTALL"; return 1; }

  cdfolder "${ROOTDIR}" || return 1;

  pbottom 'INSTALLING HDF5 LIBRARY DONE'

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  unset_env_vars_shdf5 || return 1;

  pbottom2 'SETUP_HDF5'

fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------