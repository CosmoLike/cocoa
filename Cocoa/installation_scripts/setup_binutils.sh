#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_DISTUTILS_INSTALLATION}" ]; then
    
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1
  fi
  
  # Parenthesis = run in a subshell
  ( source "${ROOTDIR:?}/installation_scripts/.check_flags.sh" ) || return 1;
     
  unset_env_vars () {
    unset -v CCIL PACKDIR
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
    fail_script_msg "setup_binutils.sh" "${1}"
    unset_all || return 1
  }

  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
  }

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  ptop2 'SETUP_BINUTILS' || return 1
  
  unset_env_vars || return 1;

  CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

  # ----------------------------------------------------------------------------
  # ----------------------------- TEXINFO LIBRARY  -----------------------------
  # ----------------------------------------------------------------------------
  
  ptop 'INSTALLING TEXINFO LIBRARY' || return 1;

  PACKDIR="${CCIL:?}/${COCOA_TEXINFO_DIR:-"texinfo-7.0.3/"}"

  cdfolder "${PACKDIR}" || return 1;

  FC="${FORTRAN_COMPILER:?}" CC="${C_COMPILER:?}" \
    ./configure \
    --prefix="${ROOTDIR:?}/.local" \
    --disable-perl-xs \
    >${OUT1:?} 2>${OUT2:?} || { error "CONFIGURE"; return 1; }

  make -j $MNT all \
    >${OUT1:?} 2>${OUT2:?} || { error "(TEXINFO) ${EC7:?}"; return 1; }
    
  make install \
    >${OUT1:?} 2>${OUT2:?} || { error "(TEXINFO) ${EC10:?}"; return 1; }

  cdfolder "${ROOTDIR}" || return 1;
  
  pbottom 'INSTALLING TEXINFO LIBRARY' || return 1;
  
  # ----------------------------------------------------------------------------
  # ----------------------------- DISTUTILS LIBRARY  ---------------------------
  # ----------------------------------------------------------------------------
  
  ptop 'INSTALLING BINUTILS LIBRARY' || return 1;
  
  PACKDIR="${CCIL:?}/${COCOA_BINUTILS_DIR:-"binutils-2.37/"}"

  cdfolder "${PACKDIR}" || return 1;

  FC="${FORTRAN_COMPILER:?}" CC="${C_COMPILER:?}" \
    ./configure \
    --prefix="${ROOTDIR:?}/.local" \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC11:?}"; return 1; }

  make -j $MNT \
    >${OUT1:?} 2>${OUT2:?} || { error "(BINUTILS) ${EC7:?}"; return 1; }

  make install \
    >${OUT1:?} 2>${OUT2:?} || { error "(BINUTILS) ${EC10:?}"; return 1; }

  cdfolder "${ROOTDIR}" || return 1;

  pbottom 'INSTALLING BINUTILS LIBRARY' || return 1;
  
  # ----------------------------------------------------------------------------
  
  unset_all || return 1;
  
  pbottom2 'SETUP_BINUTILS' || return 1

fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------