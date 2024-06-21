#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_POLYCHORD_COMPILATION}" ]; then
  
  pfail() {
    echo -e "\033[0;31m ERROR ENV VARIABLE ${1} IS NOT DEFINED \033[0m"
    unset pfail
  }
  
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1
  fi

  cdroot() {
    cd "${ROOTDIR}" 2>"/dev/null" || { echo -e \
      "\033[0;31m\t\t CD ROOTDIR (${ROOTDIR}) FAILED \033[0m"; return 1; }
    unset cdroot
  }

  if [ -z "${GIT}" ]; then
    pfail 'GIT'; cdroot; return 1;
  fi

  if [ -z "${POLYCHORD_GIT_COMMIT}" ]; then
    pfail 'POLYCHORD_GIT_COMMIT'; cdroot; return 1;
  fi

  unset_env_vars_spoly () {
    cd $ROOTDIR
    unset OUT1
    unset OUT2
    unset POLY_URL
    unset CIL
    unset ECODEF
    unset POLY_CHANGES
    unset pfail
    unset unset_env_vars_spoly
    cdroot || return 1;
  }

  fail_spoly () {
    local MSG="\033[0;31m (setup_polychord.sh) WE CANNOT RUN \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"  
    unset fail_spoly
    unset_env_vars_spoly
  }
  
  if [ -z "${DEBUG_POLY_OUTPUT}" ]; then
    export OUT1="/dev/null"; export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"; export OUT2="/dev/tty"
  fi

  cdfolder() {
    cd "${1}" 2>"/dev/null" || { fail_spoly "CD FOLDER: ${1}"; return 1; }
  }

  ptop2 'SETUP_POLYCHORD'

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  ptop  'INSTALLING POLYCHORD'

  if [ -z "${POLY_NAME}" ]; then
    pfail 'POLY_NAME'; cdroot; return 1;
  fi

  export POLY_URL="https://github.com/PolyChord/PolyChordLite.git"
  
  export CIL="${ROOTDIR}/../cocoa_installation_libraries"
  
  export POLY_CHANGES="${CIL}/polychord_changes"

  export ECODEF="${ROOTDIR}/external_modules/code"

  # ---------------------------------------------------------------------------
  # in case this script is called twice
  # ---------------------------------------------------------------------------
  rm -rf "${ECODEF}/${POLY_NAME}"

  cdfolder "${ECODEF}" || return 1;

  $GIT clone "${POLY_URL}" "${POLY_NAME}" >${OUT1} 2>${OUT2} || 
    { fail_spoly "GIT CLONE"; return 1; }
  
  cdfolder "${ECODEF}/${POLY_NAME}" || return 1;

  $GIT checkout $POLYCHORD_GIT_COMMIT >${OUT1} 2>${OUT2} ||
    { fail_spoly "GIT CHECKOUT"; return 1; }
  
  cp "${POLY_CHANGES}"/Makefile.patch . 2>${OUT2} ||
    { fail_spoly "CP FILE PATCH (Makefile.patch)"; return 1; }
  
  patch -u Makefile -i Makefile.patch >${OUT1} 2>${OUT2} ||
    { fail_spoly "PATCH FILE (Makefile.patch)"; return 1; }
  
  cp "${POLY_CHANGES}/setup.patch" . 2>${OUT2} ||
    { fail_spoly "CP FILE PATCH (SETUP)"; return 1; }
  
  patch -u setup.py -i setup.patch >${OUT1} 2>${OUT2} ||
    { fail_spoly "PATCH FILE (SETUP)"; return 1; }

  cdfolder "${ROOTDIR}" || return 1;
  
  pbottom 'INSTALLING POLYCHORD'
  
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  unset_env_vars_spoly
  
  pbottom2 'SETUP_POLYCHORD'

fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------