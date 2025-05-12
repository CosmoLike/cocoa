#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

if [ -z "${ROOTDIR}" ]; then
  source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
fi

# parenthesis = run in a subshell  
( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

unset_env_vars () {
  unset -v TMP TMP2 FILE FOLDER
  cdroot || return 1;
}

unset_env_funcs () {
  unset -f cdfolder cpfolder error warning
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
  fail_script_msg "$(basename "${BASH_SOURCE[0]}")" "${1}"
  unset_all || return 1
}

cdfolder() {
  cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER ${1}"; return 1; }
}

warning () {
  warning_script_msg "$(basename "${BASH_SOURCE[0]}")" "${1}"
}
 
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

unset_env_vars || return 1

# ------------------------------------------------------------------------------

for TMP in $(find "${ROOTDIR:?}/projects/" -mindepth 1 -maxdepth 1 -type d ! -name 'example'); do

  TMP2=$(echo "${TMP:?}" | sed -E "s@${ROOTDIR:?}/projects/@@")

  FOLDER="${ROOTDIR:?}/projects/${TMP2:?}/scripts"

  if [ ! -d "${FOLDER:?}" ]; then
    
    warning "${EC31:?} (${FOLDER:?})";
  
  else

    FILE="${FOLDER:?}/compile_${TMP2:?}.sh"

    if [ -f "${FILE:?}" ]; then
      ( export LD_LIBRARY_PATH=${CONDA_PREFIX:?}/lib:$LD_LIBRARY_PATH && \
        source "${FILE:?}" ) || { error "${EC31:?} (${FILE:?})"; return 1; }
    fi

  fi

done

# ------------------------------------------------------------------------------

unset_all

cd ${ROOTDIR:?}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------