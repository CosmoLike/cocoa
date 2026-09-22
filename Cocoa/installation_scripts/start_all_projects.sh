#!/bin/bash
# ------------------------------------------------------------------------------
# start_all_projects.sh: run every project's start script (environment
# variables and symlinks).
#
# Sourced, never executed: use `return`, not `exit`, and unset everything
# defined here before returning.
# ------------------------------------------------------------------------------

# parenthesis = run in a subshell  
( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

# unset_env_vars: forget every variable the body defines. Sourced scripts
# share the caller's shell, so anything not unset here would leak into the
# user's environment. cdroot returns to the directory the user launched from.
unset_env_vars () {
  unset -v TMP TMP2 FILE FOLDER
  unset -v 
  cdroot || return 1;
}

# unset_env_funcs: forget every helper function defined below (functions leak
# into the user's shell exactly like variables).
unset_env_funcs () {
  unset -f cdfolder cpfolder cpfile error warning
  unset -f unset_env_funcs 
  cdroot || return 1;
}

# unset_all: the single cleanup entry point: variables, then functions, then
# itself.
unset_all () {
  unset_env_vars
  unset_env_funcs
  unset -f unset_all
  cdroot || return 1;
}

# error: print which step failed (fail_script_msg names this script and the
# message) and run the FULL cleanup. A call site whose helper already routes
# through error() therefore uses plain `|| return 1`; adding another unset_all
# there would run the cleanup twice.
error () {
  fail_script_msg "$(basename "${BASH_SOURCE[0]}")" "${1}"
  unset_all || return 1
}

# cdfolder: cd that reports the target and cleans up on failure.
cdfolder() {
  cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER ${1}"; return 1; }
}

warning () {
  warning_script_msg "$(basename "${BASH_SOURCE[0]}")" "${1}"
}

# cpfile: file copy that reports source/target and cleans up on failure.
cpfile() {
  cp -r "${1:?}" "${2:?}" \
    2>"/dev/null" || { error "CP FILE ${1} on ${2}"; return 1; }
}

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

unset_env_vars || return 1

# ---------------------------------------------------------------------------

rm -f "${ROOTDIR:?}/projects/.gitignore"
rm -f "${ROOTDIR:?}/external_modules/data/.gitignore"
rm -f "${ROOTDIR:?}/external_modules/code/.gitignore"
rm -f "${ROOTDIR:?}/cobaya/cobaya/likelihoods/.gitignore"

# create tmp .gitignore and include all projects in it
touch "${ROOTDIR:?}/projects/.gitignore"

for TMP in $(find "${ROOTDIR:?}/projects" -mindepth 1 -maxdepth 1 -type d ! -name 'example'); do
  TMP2=$(echo "${TMP:?}" | sed -E "s@${ROOTDIR:?}/projects/@@")

  FOLDER="${ROOTDIR:?}/projects/${TMP2:?}/scripts"

  if [ ! -d "${FOLDER:?}" ]; then
    warning "${EC31:?} (${FOLDER:?})" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC34:?}"; return 1; }
  else
    FILE="${FOLDER:?}/start_${TMP2:?}.sh"

    if [ -f "${FILE:?}" ]; then
      source "${FILE:?}" || { error "${EC31:?} (${FILE:?})"; return 1; }
    fi
  fi

  echo "${TMP2:?}" >> "${ROOTDIR:?}/projects/.gitignore" || 
    { error "${EC35:?}"; return 1; }
done

unset_env_vars || return 1

( source "${ROOTDIR:?}/installation_scripts/create_link_all_projects.sh" )

FILE="${ROOTDIR:?}"/projects/.gitignore

cpfile "${FILE:?}" "${ROOTDIR:?}"/external_modules/data/ || return 1;

cpfile "${FILE:?}" "${ROOTDIR:?}"/external_modules/code/ || return 1;

# ------------------------------------------------------------------------------

cdfolder ${ROOTDIR:?} || return 1;

unset_all

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------