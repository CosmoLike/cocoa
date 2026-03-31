#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -n "${IGNORE_ACTDR6_CODE}" ]; then
  return 99
fi

if [ -z "${ROOTDIR}" ]; then
  source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
fi

# parenthesis = run in a subshell
( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

unset_env_vars () {
  unset -v URL CCIL ECODEF FOLDER PACKDIR CHANGES TFOLDER TFILE TFILEP AL
  unset -v EDATAF COB COBLIKE
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

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

unset_env_vars || return 1

CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

# E = EXTERNAL, CODE, F=FODLER
ECODEF="${ROOTDIR:?}/external_modules/code"

COB="${ROOTDIR:?}/cobaya"        # COB = Cobaya

COBLIKE="cobaya/likelihoods"      # COB = Cobaya, LIKE = likelihoods

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

ptop 'SETUP ACTDR6 (CMBONLY)' || { unset_all; return 1; }

URL="${ACTDR6_CMBONLY_URL:-"https://github.com/ACTCollaboration/DR6-ACT-lite.git"}"

CHANGES="${CCIL:?}/act_dr6_cmbonly_changes"

FOLDER="${ACTDR6_CMBONLY_NAME:-"act_dr6_cmbonly"}"

PACKDIR="${ECODEF:?}/${FOLDER:?}"




if [[ -n "${OVERWRITE_EXISTING_ACTDR6_CMB_CODE:-}" ]]; then

  rm -rf "${PACKDIR:?}"

fi

if [[ ! -d "${PACKDIR:?}" ]]; then

  cdfolder "${ECODEF:?}" || { unset_all; return 1; }

  "${GIT:?}" clone "${URL:?}" --depth ${GIT_CLONE_MAXIMUM_DEPTH:-1000} \
    --recursive "${FOLDER:?}" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }
  
  cdfolder "${PACKDIR:?}" || { unset_all; return 1; }

  if [ -n "${ACTDR6_CMBONLY_GIT_COMMIT:-}" ]; then

    if [ "$("${GIT:?}" rev-parse --is-shallow-repository)" = "true" ]; then
    
      "${GIT:?}" fetch --unshallow --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    
    else
    
      "${GIT:?}" fetch --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    
    fi
    
    "${GIT:?}" checkout "${ACTDR6_CMBONLY_GIT_COMMIT:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  
  fi

  # --------------------------------------------------------------------------
  # Patch code to be compatible w/ COCOA environment -------------------------
  # --------------------------------------------------------------------------
  # T = TMP
  declare -a TFOLDER=("act_dr6_cmbonly/") # If nonblank, path must include /
  declare -a TFILE=("act_dr6_cmbonly.py")
  declare -a TFILEP=("act_dr6_cmbonly.patch")
  # AL = Array Length
  AL=${#TFOLDER[@]}

  for (( i=0; i<${AL}; i++ ));
  do
    cdfolder "${PACKDIR:?}/${TFOLDER[$i]}" || { unset_all; return 1; }

    cpfolder "${CHANGES:?}/${TFOLDER[$i]}${TFILEP[$i]:?}" . \
      2>>${OUT2:?} || { unset_all; return 1; }

    patch -u "${TFILE[$i]:?}" -i "${TFILEP[$i]:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC17:?} (${TFILE[$i]:?})"; return 1; }
  done
fi

pbottom 'SETUP ACTDR6 (CMBONLY)' || { unset_all; return 1; }

cdfolder "${ROOTDIR}" || { unset_all; return 1; }

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

ptop 'SETUP ACTDR6 (MFLIKE)' || { unset_all; return 1; }

URL="${ACTDR6_MFLIKE_URL:-"https://github.com/ACTCollaboration/act_dr6_mflike.git"}"

CHANGES="${CCIL:?}/act_dr6_mflike_changes"

FOLDER="${ACTDR6_MFLIKE_NAME:-"act_dr6_mflike"}"

PACKDIR="${ECODEF:?}/${FOLDER:?}"

if [[ -n "${OVERWRITE_EXISTING_ACTDR6_CMB_CODE:-}" ]]; then
  
  rm -rf "${PACKDIR:?}"

fi

if [[ ! -d "${PACKDIR:?}" ]]; then

  cdfolder "${ECODEF:?}" || { unset_all; return 1; }

  "${GIT:?}" clone "${URL:?}" --depth ${GIT_CLONE_MAXIMUM_DEPTH:-1000} \
    --recursive "${FOLDER:?}" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }
  
  cdfolder "${PACKDIR:?}" || { unset_all; return 1; }

  if [ -n "${ACTDR6_MFLIKE_GIT_COMMIT:-}" ]; then
  
    if [ "$("${GIT:?}" rev-parse --is-shallow-repository)" = "true" ]; then
    
      "${GIT:?}" fetch --unshallow --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    
    else
    
      "${GIT:?}" fetch --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    
    fi

    "${GIT:?}" checkout "${ACTDR6_MFLIKE_GIT_COMMIT:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
  
  fi
fi

pbottom 'SETUP ACTDR6 (MFLIKE)' || { unset_all; return 1; }

cdfolder "${ROOTDIR}" || { unset_all; return 1; }

unset_all || return 1

#-----------------------------------------------------------------------------

return 55; # why this odd number? Setup_cocoa will cache this installation only
           #   if this script runs entirely. What if the user close the terminal 
           #   or the system shuts down in the middle of a git clone?  
           #   In this case, PACKDIR would exists, but it is corrupted

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------