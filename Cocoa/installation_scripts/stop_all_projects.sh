#!/bin/bash
# ------------------------------------------------------------------------------
# stop_all_projects.sh: run every project's stop script, undoing
# start_all_projects.sh.
#
# Sourced, never executed: use `return`, not `exit`, and unset everything
# defined here before returning.
# ------------------------------------------------------------------------------

if [ -z "${ROOTDIR}" ]; then
  ROOTDIR=$(pwd -P) || { echo -e \
  "\033[0;31m       ERROR ENV VARIABLE ROOTDIR NOT DEFINED \033[0m"; return 1; }
fi

fail_script_msg2 () {
  local MSG="\033[0;31m        (${1:-"empty arg"}) we cannot run \033[3m"
  local MSG2="\033[0m"
  echo -e "${MSG} ${2:-"empty arg"} ${MSG2}"
}

warning_script_msg2 () {
  local MSG="\033[0;31m        (${1:-"empty arg"}) warning: \033[3m"
  local MSG2="\033[0m"
  echo -e "${MSG} ${2:-"empty arg"} ${MSG2}"
}

# unset_env_vars: forget every variable the body defines. Sourced scripts
# share the caller's shell, so anything not unset here would leak into the
# user's environment. cdroot returns to the directory the user launched from.
unset_env_vars () {
  unset -v TMP TMP2 FILE FOLDER 
  cd "${ROOTDIR:?}" 2>"/dev/null" || { echo -e \
    "\033[0;31m\t\t CD ROOTDIR (${ROOTDIR}) FAILED \033[0m"; return 1; }
}

# unset_env_funcs: forget every helper function defined below (functions leak
# into the user's shell exactly like variables).
unset_env_funcs () {
  unset -f cdfolder cpfolder cpfile error warning fail_script_msg2
  unset -f unset_env_funcs warning_script_msg2
  cd "${ROOTDIR:?}" 2>"/dev/null" || { echo -e \
    "\033[0;31m\t\t CD ROOTDIR (${ROOTDIR}) FAILED \033[0m"; return 1; }
}

# unset_all: the single cleanup entry point: variables, then functions, then
# itself.
unset_all () {
  unset_env_vars
  unset_env_funcs
  unset -f unset_all
  cd "${ROOTDIR:?}" 2>"/dev/null" || { echo -e \
    "\033[0;31m\t\t CD ROOTDIR (${ROOTDIR}) FAILED \033[0m"; return 1; }
}

# error: print which step failed (fail_script_msg names this script and the
# message) and run the FULL cleanup. A call site whose helper already routes
# through error() therefore uses plain `|| return 1`; adding another unset_all
# there would run the cleanup twice.
error () {
  fail_script_msg2 "$(basename "${BASH_SOURCE[0]}")" "${1}"
  unset_all || return 1
}

# cdfolder: cd that reports the target and cleans up on failure.
cdfolder() {
  cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER ${1}"; return 1; }
}

warning () {
  warning_script_msg2 "$(basename "${BASH_SOURCE[0]}")" "${1}"
}




# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

unset_env_vars || return 1

# ---------------------------------------------------------------------------

for TMP in $(find "${ROOTDIR:?}/projects" -mindepth 1 -maxdepth 1 -type d ! -name 'example'); do
  TMP2=$(echo "${TMP:?}" | sed -E "s@${ROOTDIR:?}/projects/@@")
  
  FOLDER="${ROOTDIR:?}/projects/${TMP2:?}/scripts"

  if [ ! -d "${FOLDER:?}" ]; then
    warning "${EC31:?} (${FOLDER:?})" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC34:?}"; return 1; }
  else
    FILE="${FOLDER:?}/stop_${TMP2:?}.sh"

    if [ -f "${FILE:?}" ]; then
      source "${FILE:?}"  || { error "${EC31:?} (${FILE:?})"; return 1; }
    fi
  fi
done

unset_env_vars || return 1

( source "${ROOTDIR:?}/installation_scripts/delete_link_all_projects.sh" )

rm -f "${ROOTDIR:?}/projects/.gitignore"
rm -f "${ROOTDIR:?}/external_modules/data/.gitignore"
rm -f "${ROOTDIR:?}/external_modules/code/.gitignore"
rm -f "${ROOTDIR:?}/cobaya/cobaya/likelihoods/.gitignore"

# ------------------------------------------------------------------------------

cd ${ROOTDIR:?}

unset_all

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------