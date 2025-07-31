#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_DERIVKIT_CODE}" ]; then

  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi

  # parenthesis = run in a subshell
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

  unset_env_vars () {
    unset -v URL URL CCIL ECODEF FOLDER PACKDIR  
    cdroot || return 1;
  }

  unset_env_funcs () {
    unset -f cdfolder cpfolder error cpfile
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
    cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
  }

  cpfolder() {
    cp -r "${1:?}" "${2:?}"  \
      2>"/dev/null" || { error "CP FOLDER ${1} on ${2}"; return 1; }
  }

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------

  unset_env_vars || return 1

  CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

  # E = EXTERNAL, CODE, F=FODLER
  ECODEF="${ROOTDIR:?}/external_modules/code"

  URL="${DERIVKIT_URL:-"https://github.com/nikosarcevic/derivkit.git"}"
  
  FOLDER="${DERIVKIT_NAME:-"derivkit"}"

  PACKDIR="${ECODEF:?}/${FOLDER:?}"
  
  ptop "INSTALLING DERIVKIT DERIVATIVE CALCULATOR" || return 1;

  if [ -n "${OVERWRITE_EXISTING_DERIVKIT_CODE}" ]; then
    rm -rf "${PACKDIR:?}"
  fi

  if [ ! -d "${PACKDIR:?}" ]; then
    cdfolder "${ECODEF:?}" || return 1;

    "${GIT:?}" clone "${URL:?}" --depth ${GIT_CLONE_MAXIMUM_DEPTH:?} --recursive \
      "${FOLDER:?}" >${OUT1:?} 2>${OUT2:?} || { error "${EC15:?}"; return 1; }
  
    cdfolder "${PACKDIR:?}" || { cdroot; return 1; }

    if [ -n "${DERIVKIT_GIT_COMMIT}" ]; then
      "${GIT:?}" checkout "${DERIVKIT_GIT_COMMIT:?}" \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi  
  fi
  
  cdfolder "${ROOTDIR}" || return 1
  
  pbottom "INSTALLING DERIVKIT DERIVATIVE CALCULATOR" || return 1

  # ---------------------------------------------------------------------------

  unset_all || return 1
  
fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------