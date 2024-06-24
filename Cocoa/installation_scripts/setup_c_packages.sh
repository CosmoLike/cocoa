#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_C_INSTALLATION}" ]; then

  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1;
  fi

  # Parenthesis = run in a subshell
  ( source "${ROOTDIR:?}/installation_scripts/.check_flags.sh" ) || return 1;
    
  unset_env_vars () {
    unset -v CCIL BFD PACKDIR
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
    fail_script_msg "setup_c_packages.sh" "${1}"
    unset_all || return 1
  }
    
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
  }
  
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  ptop2 'SETUP_C_PACKAGES' || return 1

  unset_env_vars || return 1; 

  CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

  # ----------------------------------------------------------------------------
  # -------------------------------- FFTW --------------------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_C_FFTW_INSTALLATION}" ]; then
    
    ptop 'INSTALLING FFTW C LIBRARY' || return 1
    
    PACKDIR="${CCIL:?}/${COCOA_FFTW_DIR:-"fftw-3.3.10/"}"

    cdfolder "${PACKDIR}" || return 1;

    FC="${FORTRAN_COMPILER:?}" CC="${C_COMPILER:?}" ./configure \
      --enable-openmp \
      --prefix="${ROOTDIR:?}/.local" \
      --enable-shared=yes \
      --enable-static=yes \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC11:?}"; return 1; }

    make -j $MNT all >${OUT1:?} 2>${OUT2:?} || { error "${EC8:?}"; return 1; }

    make install >${OUT1:?} 2>${OUT2:?} || { error "${EC10:?}"; return 1; }

    cdfolder "${ROOTDIR}" || return 1;

    pbottom 'INSTALLING FFTW C LIBRARY' || return 1

  fi

  # ----------------------------------------------------------------------------
  # ------------------------------- CFITSIO ------------------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_C_CFITSIO_INSTALLATION}" ]; then
    
    ptop 'INSTALLING CFITSIO C LIBRARY' || return 1

    PACKDIR="${CCIL:?}/${COCOA_CFITSIO_DIR:-"cfitsio-4.0.0/"}"
    
    BFD="CFITSIOBUILD"

    rm -f  "${PACKDIR:?}/CMakeCache.txt"
    rm -rf "${PACKDIR:?}/${BFD:?}"
    
    mkdir "${PACKDIR}/${BFD:?}/" 2>${OUT2:?} || { error "${EC14:?}"; return 1; }
    
    cdfolder "${PACKDIR:?}/${BFD:?}" || return 1;

    ${CMAKE:?} -DBUILD_SHARED_LIBS=TRUE \
      -DCMAKE_INSTALL_PREFIX="${ROOTDIR:?}/.local" \
      -DCMAKE_C_COMPILER="${C_COMPILER:?}" \
      -DCMAKE_CXX_COMPILER="${CXX_COMPILER:?}" \
      -DCMAKE_FC_COMPILER="${FORTRAN_COMPILER:?}" \
      --log-level=ERROR .. \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC12:?}"; return 1; }

    make -j $MNT all >${OUT1:?} 2>${OUT2:?} || { error "${EC8:?}"; return 1; }

    make install >${OUT1:?} 2>${OUT2:?} || { error "${EC10:?}"; return 1; }
    
    rm -rf "${PACKDIR:?}/${BFD:?}"

    cdfolder "${ROOTDIR}" || return 1;

    pbottom 'INSTALLING CFITSIO C LIBRARY' || return 1

  fi

  # ----------------------------------------------------------------------------
  # ---------------------------------- GSL -------------------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_C_GSL_INSTALLATION}" ]; then
    
    ptop 'INSTALLING GSL C LIBRARY' || return 1

    PACKDIR="${CCIL:?}/${COCOA_GSL_DIR:-"gsl-2.7/"}" 

    cdfolder "${PACKDIR}" || return 1;

    CC="${C_COMPILER:?}" ./configure \
      --prefix="${ROOTDIR:?}/.local" \
      --enable-shared=yes \
      --enable-static=yes \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC11:?}"; return 1; }
 
    make -j $MNT all >${OUT1:?} 2>${OUT2:?} || { error "${EC8:?}"; return 1; }

    make install >${OUT1:?} 2>${OUT2:?} || { error "${EC10:?}"; return 1; }

    cdfolder "${ROOTDIR}" || return 1;

    pbottom 'INSTALLING GSL C LIBRARY' || return 1

  fi

  # ----------------------------------------------------------------------------
  # ----------------------------------- EUCLID EMU -----------------------------
  # ----------------------------------------------------------------------------
  # WE MIGRATED euclidemu2 TO setup_c_packages: IT DEPENDS ON GSL-GNU LIB
  if [ -z "${IGNORE_ALL_PIP_INSTALLATION}" ]; then
            
    ptop 'INSTALLING EUCLIDEMU2' || return 1

    env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" ${PIP3:?} install \
      --global-option=build_ext "${CCIL:?}/euclidemu2-1.2.0" \
      --no-dependencies \
      --prefix="${ROOTDIR:?}/.local" \
      --no-index 
      >${OUT1:?} 2>${OUT2:?} || { error "${EC13:?}"; return 1; }

    cdfolder "${ROOTDIR}" || return 1;  
    
    pbottom 'INSTALLING EUCLIDEMU2' || return 1
  
  fi

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  unset_all || return 1; 

  pbottom2 'SETUP_C_PACKAGES DONE' || return 1

fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------