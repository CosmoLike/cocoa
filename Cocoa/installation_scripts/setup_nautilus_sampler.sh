#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_NAUTILUS_SAMPLER_CODE:-}" ]; then

  if [ -z "${ROOTDIR:-}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi

  # parenthesis = run in a subshell
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

  unset_env_vars () {
    unset -v URL CCIL ECODEF FOLDER PACKDIR  
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

  URL="${NAUTILUS_SAMPLER_URL:-"https://github.com/johannesulf/nautilus.git"}"
  
  FOLDER="${NAUTILUS_SAMPLER_NAME:-"nautilus_sampler"}"

  PACKDIR="${ECODEF:?}/${FOLDER:?}"
  
  ptop "INSTALLING NAUTILUS SAMPLER" || { unset_all; return 1; }

  # ---------------------------------------------------------------------------
  # in case this script is called twice
  # ---------------------------------------------------------------------------
  if [ -n "${OVERWRITE_EXISTING_NAUTILUS_CODE:-}" ]; then
  
    rm -rf "${PACKDIR:?}"
  
  fi

  if [ ! -d "${PACKDIR:?}" ]; then
    
    cdfolder "${ECODEF:?}" || { unset_all; return 1; }

    "${GIT:?}" clone "${URL:?}" --depth ${GIT_CLONE_MAXIMUM_DEPTH:-1000} \
      --recursive "${FOLDER:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }
  
    cdfolder "${PACKDIR:?}" || { unset_all; return 1; }

    if [ -n "${NAUTILUS_SAMPLER_GIT_COMMIT:-}" ]; then
  
      if [ "$("${GIT:?}" rev-parse --is-shallow-repository)" = "true" ]; then
    
        "${GIT:?}" fetch --unshallow --all --tags --prune \
          >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    
      else
    
        "${GIT:?}" fetch --all --tags --prune \
          >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    
      fi

      "${GIT:?}" checkout "${NAUTILUS_SAMPLER_GIT_COMMIT:?}" \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    
    fi  
  fi
  
  cdfolder "${ROOTDIR}" || { unset_all; return 1; }
  
  pbottom "INSTALLING NAUTILUS SAMPLER" || { unset_all; return 1; }

  # ---------------------------------------------------------------------------

  unset_all || return 1
  
fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------