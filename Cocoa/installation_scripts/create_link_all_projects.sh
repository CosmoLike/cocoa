#!/bin/bash
# ------------------------------------------------------------------------------
# create_link_all_projects.sh: create the symlinks that expose every project's
# likelihood and interface inside cobaya.
#
# Sourced, never executed: use `return`, not `exit`, and unset everything
# defined here before returning.
# ------------------------------------------------------------------------------

if [ -z "${ROOTDIR}" ]; then
  source start_cocoa.sh
  return 1;
fi

# parenthesis = run in a subshell  
( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

# unset_env_vars: forget every variable the body defines. Sourced scripts
# share the caller's shell, so anything not unset here would leak into the
# user's environment. cdroot returns to the directory the user launched from.
unset_env_vars () {
  unset -v TMP TMP2 TMP3 FOLDER TARGET ECODEF EDATAF COB COBLIKE CCLIKE
  unset -v 
  cdroot || return 1;
}

# unset_env_funcs: forget every helper function defined below (functions leak
# into the user's shell exactly like variables).
unset_env_funcs () {
  unset -f cdfolder cpfolder error clink warning
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

warning () {
  warning_script_msg "$(basename "${BASH_SOURCE[0]}")" "${1}"
}

# cdfolder: cd that reports the target and cleans up on failure.
cdfolder() {
  cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER ${1}"; return 1; }
}

clink() {
  if [ ! -d "${1:?}" ]; then
    warning "${EC31:?} (${1:?})" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC34:?}"; return 1; }
  else
    if [ ! -d "${4:?}/${2:?}" ]; then
      ln -sf "${1:?}" "${4:?}" \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC34:?}"; return 1; }
      
      cdfolder "${4:?}" || return 1

      mv -T "${3:?}" "${2:?}" \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC31:?}"; return 1; }
    else
      rm -f "${4:?}/${2:?}"

      ln -sf "${1:?}" "${4:?}" \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC34:?}"; return 1; }
      
      cdfolder "${4:?}" || return 1

      mv -T "${3:?}" "${2:?}" \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC31:?}"; return 1; }
    fi
  fi
}

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

unset_env_vars || return 1

# ------------------------------------------------------------------------------

# E = EXTERNAL, CODE, F=FODLER
ECODEF="${ROOTDIR:?}/external_modules/code"

# E = EXTERNAL, DATA, F=FODLER
EDATAF="${ROOTDIR:?}/external_modules/data"

COB="${ROOTDIR:?}/cobaya" 

COBLIKE="cobaya/likelihoods"

COBTH="cobaya/theories"

for TMP in $(find "${ROOTDIR:?}/projects" -mindepth 1 -maxdepth 1 -type d ! -name 'example'); do
  
  TMP2=$(echo "${TMP:?}" | sed -E "s@${ROOTDIR:?}/projects/@@")

  TMP3=interface
  
  FOLDER="${ROOTDIR:?}/projects/${TMP2:?}/${TMP3:?}"

  TARGET="${ECODEF:?}"

  clink "${FOLDER:?}" "${TMP2:?}" "${TMP3:?}" "${TARGET:?}" || return 1;

  # ----------------------------------------------------------------------------

  TMP3=likelihood
  
  FOLDER="${ROOTDIR:?}/projects/${TMP2:?}/${TMP3:?}"

  TARGET="${COB:?}/${COBLIKE:?}"

  clink "${FOLDER:?}" "${TMP2:?}" "${TMP3:?}" "${TARGET:?}" || return 1;

  # ----------------------------------------------------------------------------

  TMP3=data
  
  FOLDER="${ROOTDIR:?}/projects/${TMP2:?}/${TMP3:?}"

  TARGET=${EDATAF:?}

  clink "${FOLDER:?}" "${TMP2:?}" "${TMP3:?}" "${TARGET:?}" || return 1;

  # ----------------------------------------------------------------------------

  TMP3=theory

  FOLDER="${ROOTDIR:?}/projects/${TMP2:?}/${TMP3:?}"

  TARGET="${COB:?}/${COBTH:?}"

  clink "${FOLDER:?}" "${TMP2:?}" "${TMP3:?}" "${TARGET:?}" || return 1;

done

# ------------------------------------------------------------------------------

cdfolder ${ROOTDIR:?}

unset_all

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------