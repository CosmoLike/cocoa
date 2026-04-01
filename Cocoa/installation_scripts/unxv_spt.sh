#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -n "${IGNORE_SPT_CMB_DATA:-}" ]; then
  return 99
fi

if [ -z "${ROOTDIR}" ]; then
  source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
fi

# parenthesis = run in a subshell 
( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;
  
unset_env_vars () { 
  unset -v EDATAF FOLDER URL PACKDIR PRINTNAME
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
EDATAF="${ROOTDIR:?}/external_modules/data"

FOLDER="${SPT_3G_NAME:-"spt_3g"}"

# PACK = PACKAGE, DIR = DIRECTORY
PACKDIR="${EDATAF:?}/${FOLDER:?}"

URL="${SPT3G_DATA_URL:-"https://github.com/SouthPoleTelescope/spt3g_y1_dist.git"}"

# Name to be printed on this shell script messages
PRINTNAME="SPT-3G Y1"

ptop "SETUP/UNXV ${PRINTNAME:?} DATA" || { unset_all; return 1; }

# ---------------------------------------------------------------------------
# in case this script is called twice
# ---------------------------------------------------------------------------
if [ -n "${OVERWRITE_EXISTING_SPT3G_CMB_DATA:-}" ]; then

  rm -rf "${PACKDIR:?}"

fi

if [ ! -d "${PACKDIR:?}" ]; then
  
  cdfolder "${EDATAF:?}" || { unset_all; return 1; }
    
  "${GIT:?}" clone "${URL:?}" --depth ${GIT_CLONE_MAXIMUM_DEPTH:-1000} \
    --recursive "${FOLDER:?}" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }

  if [ -n "${SPT3G_DATA_GIT_COMMIT:-}" ]; then
    
    cdfolder "${PACKDIR:?}" || { unset_all; return 1; }
    
    if [ "$("${GIT:?}" rev-parse --is-shallow-repository)" = "true" ]; then
    
      "${GIT:?}" fetch --unshallow --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    
    else
    
      "${GIT:?}" fetch --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    
    fi

    "${GIT:?}" checkout "${SPT3G_DATA_GIT_COMMIT:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  
  fi

fi

cdfolder "${ROOTDIR}" || { unset_all; return 1; }

pbottom "SETUP/UNXV ${PRINTNAME:?} DATA" || { unset_all; return 1; }

#-------------------------------------------------------------------------------

unset_all || return 1; 

#-------------------------------------------------------------------------------

return 55; # why this odd number? Setup_cocoa will cache this installation only
           #   if this script runs entirely. What if the user close the terminal 
           #   or the system shuts down in the middle of a git clone?  
           #   In this case, PACKDIR would exists, but it is corrupted

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------