#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

if [ -z "${ROOTDIR}" ]; then
  ROOTDIR=$(pwd -P) || { echo -e \
  "\033[0;31m       ERROR ENV VARIABLE ROOTDIR NOT DEFINED \033[0m"; return 1; }
fi

fail_script_msg2 () {
  local MSG="\033[0;31m        (${1:-"empty arg"}) we cannot run \e[3m"
  local MSG2="\033[0m"
  echo -e "${MSG} ${2:-"empty arg"} ${MSG2}"
}

unset_env_vars () {
  unset -v TMP TMP2
  cd "${ROOTDIR:?}" 2>"/dev/null" || { echo -e \
    "\033[0;31m\t\t CD ROOTDIR (${ROOTDIR}) FAILED \033[0m"; return 1; }
}

unset_env_funcs () {
  unset -f cdfolder cpfolder error warning fail_script_msg2
  unset -f unset_env_funcs 
  cd "${ROOTDIR:?}" 2>"/dev/null" || { echo -e \
    "\033[0;31m\t\t CD ROOTDIR (${ROOTDIR}) FAILED \033[0m"; return 1; }
}

unset_all () {
  unset_env_vars
  unset_env_funcs
  unset -f unset_all
  cd "${ROOTDIR:?}" 2>"/dev/null" || { echo -e \
    "\033[0;31m\t\t CD ROOTDIR (${ROOTDIR}) FAILED \033[0m"; return 1; }
}

error () {
  fail_script_msg2 "$(basename "${BASH_SOURCE[0]}")" "${1}"
  unset_all || return 1
}

cdfolder() {
  cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER ${1}"; return 1; }
}

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

unset_env_vars || return 1

# ---------------------------------------------------------------------------

for TMP in $(find "${ROOTDIR:?}/projects" -mindepth 1 -maxdepth 1 -type d ! -name 'example'); do
  
  if [ -n ${TMP2} ]; then
    
    TMP2=$(echo "${TMP}" | sed -E "s@${ROOTDIR:?}/projects/@@")

    if [ -n ${TMP2} ]; then

      rm -f "${ROOTDIR:?}/external_modules/code/${TMP2}"

      rm -f "${ROOTDIR:?}/cobaya/cobaya/likelihoods/${TMP2}"

      rm -f "${ROOTDIR:?}/external_modules/data/${TMP2}"
    
    fi

  fi

done

# ------------------------------------------------------------------------------

cdfolder ${ROOTDIR:?}

unset_all

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------