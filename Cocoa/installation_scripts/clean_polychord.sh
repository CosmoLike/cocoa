#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_POLYCHORD_COMPILATION}" ]; then
  pfail() {
    echo -e "\033[0;31m ERROR ENV VARIABLE ${1} NOT DEFINED \033[0m"
    unset pfail
  }
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'
    return 1
  fi
  cdroot() {
    cd "${ROOTDIR}" 2>"/dev/null" || { echo -e \
      "\033[0;31m\t\t CD ROOTDIR (${ROOTDIR}) FAILED \033[0m"; return 1; }
    unset cdroot
  }
  if [ -z "${PYTHON3}" ]; then
    pfail "PYTHON3"; cdroot; return 1;
  fi
  if [ -z "${PYTHON_VERSION}" ]; then
    pfail "PYTHON_VERSION"; cdroot; return 1;
  fi
  if [ -z "${POLY_NAME}" ]; then
    pfail 'POLY_NAME'; cdroot; return 1;
  fi
  unset_env_vars_clean_poly () {
    unset OUT1
    unset OUT2
    unset PLIB
    unset pfail
    unset unset_env_vars_clean_poly
    cdroot || return 1;
  }
  fail_clean_poly () {
    local MSG="\033[0;31m\t\t (clean_polychord.sh) WE CANNOT RUN \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"
    unset fail_clean_poly
    unset_env_vars_clean_poly
    return 1
  }
  cdfolder() {
    cd "${1}" 2>"/dev/null" || { fail_clean_poly "CD FOLDER: ${1}"; return 1; }
  }
  if [ -z "${DEBUG_POLY_OUTPUT}" ]; then
    export OUT1="/dev/null"
    export OUT2="/dev/null"
  else
    export OUT1="/dev/tty"
    export OUT2="/dev/tty"
  fi

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  ptop 'CLEANING POLYCHORD'

  export PLIB="${ROOTDIR}/.local/lib/python${PYTHON_VERSION}/site-packages"
  rm -rf "${PLIB}"/pypolychord-*

  cdfolder "${ROOTDIR}/external_modules/code/${POLY_NAME}" || return 1

  make clean >${OUT1} 2>${OUT2} || { fail_clean_poly "MAKE CLEAN"; return 1; }

  rm -rf ./lib/*.a
  rm -rf ./lib/*.so

  unset_env_vars_clean_poly || return 1
  pbottom 'CLEANING POLYCHORD'
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------