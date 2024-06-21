#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_C_INSTALLATION}" ]; then

  pfail() {
    echo -e \
    "\033[0;31m\t\t ERROR ENV VARIABLE ${1:-"empty arg"} NOT DEFINED \033[0m"
    unset pfail
  }
  
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1;
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
  
  unset_env_vars_scp () {
    unset OUT1
    unset OUT2
    unset CMNT
    unset CCIL
    unset pfail
    unset BFD
    unset unset_env_vars_scp
    cdroot || return 1;
  }
  
  fail_scp () {
    local MSG="\033[0;31m\t\t (setup_c_packages.sh) WE CANNOT RUN \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG} ${1:-"empty arg"} ${MSG2}"
    unset fail_scp
    unset_env_vars_scp
  }
  
  if [ -z "${DEBUG_C_PACKAGES}" ]; then
    export OUT1="/dev/null"; export OUT2="/dev/null"
    export CMNT="${MAKE_NUM_THREADS:-1}"
  else
    export OUT1="/dev/tty"; export OUT2="/dev/tty"
    export CMNT=1
  fi
  
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { fail_scp "CD FOLDER: ${1}"; return 1; }
  }
  
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  ptop2 'SETUP_C_PACKAGES'

  export CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

  # ----------------------------------------------------------------------------
  # -------------------------------- FFTW --------------------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_C_FFTW_INSTALLATION}" ]; then
    
    ptop 'INSTALLING FFTW C LIBRARY'
    
    if [ -z "${COCOA_FFTW_DIR}" ]; then
      pfail 'COCOA_FFTW_DIR'; cdroot; return 1;
    fi
    
    cdfolder "${CCIL:?}/${COCOA_FFTW_DIR:?}" || return 1;

    FC="${FORTRAN_COMPILER:?}" CC="${C_COMPILER:?}" ./configure \
      --enable-openmp \
      --prefix="${ROOTDIR:?}/.local" \
      --enable-shared=yes \
      --enable-static=yes \
      >${OUT1:?} 2>${OUT2:?} || { fail_scp "(FFTW) CONFIGURE"; return 1; }

    make -j $CMNT all \
      >${OUT1:?} 2>${OUT2:?} || { fail_scp "(FFTW) MAKE"; return 1; }

    make install \
      >${OUT1:?} 2>${OUT2:?} || { fail_scp "(FFTW) MAKE INSTALL"; return 1; }

    cdfolder "${ROOTDIR}" || return 1;

    pbottom 'INSTALLING FFTW C LIBRARY'

  fi

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  # ----------------------------------------------------------------------------
  # ------------------------------- CFITSIO ------------------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_C_CFITSIO_INSTALLATION}" ]; then
    
    ptop 'INSTALLING CFITSIO C LIBRARY'

    if [ -z "${COCOA_CFITSIO_DIR}" ]; then
      pfail 'COCOA_CFITSIO_DIR'; cdroot; return 1;
    fi

    cdfolder "${CCIL:?}/${COCOA_CFITSIO_DIR:?}" || return 1;
    
    rm -f CMakeCache.txt
    
    export BFD="CFITSIOBUILD"
    
    rm -rf "./${BFD:?}"

    mkdir "${BFD:?}/" >${OUT1:?} 2>${OUT2:?} || 
      { fail_scp "(CFITSIO) MKDIR BUILD FOLDER"; return 1; }
    
    cdfolder "${CCIL:?}/${COCOA_CFITSIO_DIR:?}/${BFD:?}" || return 1;

    "${CMAKE:?}" -DBUILD_SHARED_LIBS=TRUE \
      -DCMAKE_INSTALL_PREFIX="${ROOTDIR:?}/.local" \
      -DCMAKE_C_COMPILER="${C_COMPILER:?}" \
      -DCMAKE_CXX_COMPILER="${CXX_COMPILER:?}" \
      -DCMAKE_FC_COMPILER="${FORTRAN_COMPILER:?}" \
      --log-level=ERROR .. \
      >${OUT1:?} 2>${OUT2:?} || { fail_scp "(CFITSIO) CMAKE"; return 1; }

    make -j $CMNT all >${OUT1:?} 2>${OUT2:?} || 
      { fail_scp "(CFITSIO) MAKE"; return 1; }

    make install >${OUT1:?} 2>${OUT2:?} || 
      { fail_scp "(CFITSIO) MAKE INSTALL"; return 1; }

    cdfolder "${CCIL:?}/${COCOA_CFITSIO_DIR:?}" || return 1;
    
    rm -rf "./${BFD:?}"

    cdfolder "${ROOTDIR}" || return 1;

    pbottom 'INSTALLING CFITSIO C LIBRARY'

  fi

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  # ----------------------------------------------------------------------------
  # ---------------------------------- GSL -------------------------------------
  # ----------------------------------------------------------------------------
  
  if [ -z "${IGNORE_C_GSL_INSTALLATION}" ]; then
    
    ptop 'INSTALLING GSL C LIBRARY'

    if [ -z "${COCOA_GSL_DIR}" ]; then
      pfail 'COCOA_GSL_DIR'; cdroot; return 1;
    fi

    cdfolder "${CCIL:?}/${COCOA_GSL_DIR:?}" || return 1;

    CC=$C_COMPILER ./configure \
      --prefix="${ROOTDIR:?}/.local" \
      --enable-shared=yes \
      --enable-static=yes \
      >${OUT1:?} 2>${OUT2:?} || { fail_scp "(GSL) CONFIGURE"; return 1; }
 
    make -j $CMNT all \
      >${OUT1:?} 2>${OUT2:?} || { fail_scp "(GSL) MAKE"; return 1; }

    make install \
      >${OUT1:?} 2>${OUT2:?} || { fail_scp "(GSL) MAKE INSTALL"; return 1; }

    cdfolder "${ROOTDIR}" || return 1;

    pbottom 'INSTALLING GSL C LIBRARY'

  fi

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  unset_env_vars_scp || return 1; 

  pbottom2 'SETUP_C_PACKAGES DONE'
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# ----------------------------------- EUCLID EMU -------------------------------
# ------------------------------------------------------------------------------
# WE MIGRATED euclidemu2 TO setup_c_packages: IT DEPENDS ON GSL-GNU LIB
if [ -z "${IGNORE_ALL_PIP_INSTALLATION}" ]; then
  
  if [ -z "${PIP3}" ]; then
    pfail 'PIP3'; cdroot; return 1;
  fi
  
  unset_env_vars_scp_eemul2 () {
    unset OUT1
    unset OUT2
    cdroot || return 1;
  }
  
  fail_scp () {
    local MSG="\033[0;31m\t\t (setup_c_packages.sh) we cannot run \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG} ${1:-"empty arg"} ${MSG2}"
    unset fail_scp
    unset_env_vars_scp_eemul2
  }
  
  if [ -z "${DEBUG_PIP_PACKAGES}" ]; then
    export OUT1="/dev/null"; export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"; export OUT2="/dev/tty"
  fi

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  ptop2 'INSTALLING EUCLIDEMU2 DONE'

  env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" "${PIP3:?}" install \
    --global-option=build_ext \
    "${ROOTDIR:?}/../cocoa_installation_libraries/euclidemu2-1.2.0" \
    --no-dependencies \
    --prefix="${ROOTDIR:?}/.local" \
    --no-index >${OUT1:?} 2>${OUT2:?} || 
    { fail_scp "(EUCLIDEMU2) PIP INSTALL EUCLUDEMUL2"; return 1; }

  cdfolder "${ROOTDIR}" || return 1;  
  
  unset_env_vars_scp_eemul2 || return 1;

  pbottom2 'INSTALLING EUCLIDEMU2 DONE'
  
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------