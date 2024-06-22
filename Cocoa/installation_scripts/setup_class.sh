#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_CLASS_COMPILATION}" ]; then
  
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1;
  fi
  
  source "${ROOTDIR:?}/installation_scripts/.check_flags.sh" || return 1;
    
  unset_env_vars () {
    unset -v URL CHANGES ECODEF PACKDIR CLNAME
    cdroot || return 1;
  }
  
  error () {
    fail_script_msg "setup_class.sh" "${1}"
    unset -f error
    unset_env_vars || return 1
  }

  cdfolder() {
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
  }

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  
  ptop2 'SETUP_CLASS' || return 1
  
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------

  ptop 'INSTALLING CLASS' || return 1

  URL="https://github.com/lesgourg/class_public.git"
  
  CHANGES="${ROOTDIR:?}/../cocoa_installation_libraries/class_changes"

  ECODEF="${ROOTDIR:?}/external_modules/code"

  CLNAME=${CLASS_NAME:-"class_public"}
  
  PACKDIR="${ECODEF:?}/${CLNAME:?}/"

  # ---------------------------------------------------------------------------
  # in case this script is called twice
  # ---------------------------------------------------------------------------
  rm -rf "${PACKDIR:?}"

  # ---------------------------------------------------------------------------
  # clone from original repo
  # ---------------------------------------------------------------------------
  cdfolder "${ECODEF}" || return 1;

  ${GIT:?} clone "${URL:?}" --recursive "${CLNAME:?}" \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC15:?}"; return 1; }
  
  cdfolder "${PACKDIR}" || return 1;

  if [ -n "${CLASS_GIT_COMMIT}" ]; then
    ${GIT:?} checkout "${CLASS_GIT_COMMIT:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC16\?}"; return 1; }
  fi
  
  # ---------------------------------------------------------------------------
  # historical reasons (we used to save class_python on Cocoa Branch)
  # historical: Workaround Cocoa .gitignore entry on /include
  # ---------------------------------------------------------------------------
  mv ./include ./include2/ \
    >${OUT1:?} 2>${OUT2:?} || { error "MKDIR FROM INCLUDE TO INCLUDE2"; return 1; }
  
  # ---------------------------------------------------------------------------
  # patch CLASS to be compatible w/ COCOA
  # ---------------------------------------------------------------------------
  
  # ---------------------------------------------------------------------------
  cdfolder "${PACKDIR}" || return 1;

  cp "${CHANGES:?}/Makefile.patch" .  2>${OUT2:?} || 
    { error "CP FILE PATCH (Makefile.patch)"; return 1; }
  
  patch -u Makefile -i Makefile.patch >${OUT1:?} 2>${OUT2:?} ||
    { error "SCRIPT FILE PATCH (Makefile.patch)"; return 1; }
  
  # ---------------------------------------------------------------------------
  cdfolder "${PACKDIR:?}/python" || return 1;
  
  cp "${CHANGES:?}/python/setup.patch" . 2>${OUT2:?} ||
    { error "CP FILE PATCH (setup.patch)"; return 1; }

  patch -u setup.py -i setup.patch >${OUT1:?} 2>${OUT2:?} || 
    { error "SCRIPT FILE PATCH (setup.patch)"; return 1; }
  
  # ---------------------------------------------------------------------------
  cdfolder "${ROOTDIR}" || return 1;

  pbottom 'INSTALLING CLASS' || return 1

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  unset_env_vars || return 1

  pbottom2 'SETUP_CLASS' || return 1
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
