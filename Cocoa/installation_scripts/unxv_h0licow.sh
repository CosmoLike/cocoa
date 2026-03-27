#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_HOLICOW_STRONG_LENSING_DATA:-}" ]; then

  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi

  # parenthesis = run in a subshell 
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;
  
  unset_env_vars () {
    unset -v EDATAF FOLDER URL TMP PRINTNAME
    cdroot || return 1;
  }

  unset_env_funcs () {
    unset -f cdfolder cpfolder error
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

  # --------------------------------------------------------------------------- 
  # --------------------------------------------------------------------------- 
  # ---------------------------------------------------------------------------

  unset_env_vars || return 1

  # E = EXTERNAL, DATA, F=FODLER
  EDATAF="${ROOTDIR:?}/external_modules/data/"
  
  FOLDER="h0licow_distance_chains"

  TMP="holicow_tmp" # = TMP
  
  URL="${HOLICOW_DATA_URL:-"https://github.com/shsuyu/H0LiCOW-public.git"}"

  ptop "SETUP/UNXV H0LICOW DATA" || { unset_all; return 1; }

  # ---------------------------------------------------------------------------
  # in case this script is called twice
  # ---------------------------------------------------------------------------
  if [ -n "${OVERWRITE_EXISTING_HOLICOW_DATA:-}" ]; then
  
    rm -rf "${EDATAF:?}/${TMP:?}"
    rm -rf "${EDATAF:?}/${FOLDER:?}"
  
  fi
  
  if [ ! -d "${EDATAF:?}/${FOLDER:?}" ]; then

    cdfolder "${EDATAF:?}" || { unset_all; return 1; }
      
    "${GIT:?}" clone "${URL:?}" --depth ${GIT_CLONE_MAXIMUM_DEPTH:-1000} \
      --recursive  "${TMP:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }

    cdfolder "${EDATAF:?}/${TMP:?}" || { unset_all; return 1; }

    if [ -n "${HOLICOW_DATA_GIT_COMMIT:-}" ]; then

      if [ "$("${GIT:?}" rev-parse --is-shallow-repository)" = "true" ]; then
      
        "${GIT:?}" fetch --unshallow --all --tags --prune \
          >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
      
      else
      
        "${GIT:?}" fetch --all --tags --prune \
          >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
      
      fi
      
      "${GIT:?}" checkout "${HOLICOW_DATA_GIT_COMMIT:?}" \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi

    mv "${FOLDER:?}" "${EDATAF:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "MV H0LICOW DATA"; return 1; }
    
    rm -rf "${EDATAF:?}/${TMP:?}"
  
  fi

  pbottom "SETUP/UNXV H0LICOW DATA" || { unset_all; return 1; }

  cdfolder "${ROOTDIR}" || { unset_all; return 1; }
  
  unset_all || return 1

fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------