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
    pfail 'ROOTDIR'
    return 1
  fi
  if [ -z "${PYTHON3}" ]; then
    pfail "PYTHON3"
    cd "${ROOTDIR}" || return 1
    return 1
  fi
  if [ -z "${PYTHON_VERSION}" ]; then
    pfail "PYTHON3"
    cd "${ROOTDIR}" || return 1
    return 1
  fi
  if [ -z "${POLY_NAME}" ]; then
    pfail 'POLY_NAME'
    cd "${ROOTDIR}" || return 1
    return 1
  fi
  unset_env_vars_clean_poly () {
    unset OUT1
    unset OUT2
    unset PLIB
    unset pfail
    unset unset_env_vars_clean_poly
    cd "${ROOTDIR}" || return 1
  }
  fail_clean_poly () {
    local MSG="\033[0;31m (clean_polychord.sh) WE CANNOT RUN \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"
    unset fail_clean_poly
    unset_env_vars_clean_poly
    return 1
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

  cd "${ROOTDIR}/external_modules/code/${POLY_NAME}" 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_clean_poly "CD POLY FOLDER"
  fi

  make clean > ${OUT1} 2> ${OUT2}
  if [ $? -ne 0 ]; then
    fail_clean_poly "MAKE CLEAN"
  fi

  rm -rf ./lib/*.a
  rm -rf ./lib/*.so

  unset_env_vars_clean_poly
  pbottom 'CLEANING POLYCHORD'
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
fi
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------