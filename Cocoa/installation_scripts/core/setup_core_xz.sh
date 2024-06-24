#!/bin/bash
# ----------------------------------------------------------------------------
# -------------------------- XZ COMPRESSION Library --------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_XZ_INSTALLATION}" ]; then
  
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1
  fi

  # parenthesis = run in a subshell
  ( source "${ROOTDIR:?}/installation_scripts/.check_flags.sh" ) || return 1;

  unset_env_vars () {
    unset -v CCIL
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
    fail_script_msg "setup_xz.sh" "${1}"
    unset_env_vars
  }

  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
  }

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  ptop2 "SETUP_CORE_XZ" || return 1;

  unset_env_vars || return 1;
  
  CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

  # ----------------------------------------------------------------------------
  
  ptop "INSTALLING AND COMPILING XZ LIBRARY" || return 1;

  cdfolder "${CCIL}" || return 1;

  #False xz file: just to trigger GIT LFS
  cp xz-5.2.5.tar.gz.xz xz-5.2.5.tar.gz \
    2>${OUT2} ||  { error "CP XZ TAR"; return 1; }

  tar -xf xz-5.2.5.tar.gz.xz \
    >${OUT1} 2>${OUT2} ||  { error "TAR XZ TAR"; return 1; }

  cdfolder "${CCIL}/xz-5.2.5/" || return 1;

  CC=${C_COMPILER:?} ./configure --prefix="${ROOTDIR:?}/.local" \
    >${OUT1} 2>${OUT2} || { error "${EC11:?}"; return 1; }

  make -j $MNT all >${OUT1} 2>${OUT2} || { error "${:EC8?}"; return 1; }

  make install >${OUT1} 2>${OUT2} || { error "${EC10:?}"; return 1; }

  cdfolder "${ROOTDIR}" || return 1;

  pbottom "INSTALLING XZ LIBRARY" || return 1;
  
  # ---------------------------------------------------------------------------
  
  unset_all || return 1;

  pbottom2 "INSTALLING AND COMPILING XZ LIBRARY" || return 1;

fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------