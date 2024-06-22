#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_CAMB_COMPILATION}" ]; then

  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1;
  fi

  source "${ROOTDIR:?}/installation_scripts/.check_flags.sh" || return 1;

  unset_env_vars () {
    unset -v URL CHANGES ECODEF CAMBF PACKDIR
    cdroot || return 1;
  }

  error () {
    fail_script_msg "setup_camb.sh" "${1}"
    unset error
    unset_env_vars
  }
  
  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
  }

  cpfolder() {
    cp -r "${1:?}" 2>"/dev/null" || { error "CP FOLDER ${1} on ${2}"; return 1; }
  }

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  
  ptop2 'SETUP_CAMB' || return 1;

  unset_env_vars || return 1;

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  
  ptop 'INSTALLING CAMB' || return 1;

  URL="https://github.com/cmbant/CAMB"

  CHANGES="${ROOTDIR:?}/../cocoa_installation_libraries/camb_changes"

  ECODEF="${ROOTDIR:?}/external_modules/code"

  CAMBF=${CAMB_NAME:-"CAMB"}

  PACKDIR="${ECODEF:?}/${CAMBF:?}"

  # ---------------------------------------------------------------------------
  # in case this script is called twice
  # ---------------------------------------------------------------------------
  rm -rf "${PACKDIR:?}"

  # ---------------------------------------------------------------------------
  # clone from original repo
  # ---------------------------------------------------------------------------
  cdfolder "${ECODEF:?}" || return 1;

  ${GIT:?} clone $URL --recursive "${CAMBF:?}" \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC15:?}"; return 1; }
  
  cdfolder "${PACKDIR}" || return 1;

  if [ -n "${CAMB_GIT_COMMIT}" ]; then
    ${GIT:?} checkout "${CAMB_GIT_COMMIT:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC16:?}"; return 1; }
  fi
  
  # ---------------------------------------------------------------------------
  # patch CAMB to be compatible w/ COCOA
  # ---------------------------------------------------------------------------
  cdfolder "${PACKDIR:?}/camb/" || return 1;

  cpfolder "${CHANGES:?}"/camb/_compilers.patch . 2>${OUT2:?} || return 1;
  
  patch -u _compilers.py -i _compilers.patch \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC17:?} (_compilers)"; return 1; }

  # ---------------------------------------------------------------------------
  cdfolder "${PACKDIR:?}/fortran/" || return 1;
  
  cpfolder "${CHANGES:?}"/fortran/Makefile.patch . 2>${OUT2:?} || return 1;
  
  patch -u Makefile -i Makefile.patch \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC17:?} (Makefile)"; return 1; }

  # ---------------------------------------------------------------------------
  cdfolder "${PACKDIR:?}/forutils/" || return 1;
  
  cpfolder "${CHANGES:?}/forutils/Makefile_compiler.patch" . 2>${OUT2:?} || return 1;
  
  patch -u Makefile_compiler -i Makefile_compiler.patch \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC17:?} (Makefile_compiler)"; return 1; }
  
  # ---------------------------------------------------------------------------
  cdfolder "${ROOTDIR}" || return 1;

  pbottom 'INSTALLING CAMB' || return 1;

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  unset_env_vars || return 1;
  
  pbottom2 'SETUP_CAMB' || return 1;
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------