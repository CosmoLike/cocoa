#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_ACTDR6_CODE}" ]; then

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
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  ptop 'SETUP ACTDR6 (CMB ONLY LIKE)' || return 1;

  URL="${ACTDR6_CMBONLY_URL:-"https://github.com/ACTCollaboration/DR6-ACT-lite.git"}"

  CHANGES="${CCIL:?}/act_dr6_cmbonly_changes"

  FOLDER="act_dr6_cmbonly"

  PACKDIR="${ECODEF:?}/${FOLDER:?}"

  # ---------------------------------------------------------------------------
  # In case this script is called twice ---------------------------------------
  # ---------------------------------------------------------------------------
  if [[ -n "${OVERWRITE_EXISTING_ACTDR6_CMB_CODE}" ]]; then

    rm -rf "${PACKDIR:?}"

  fi

  if [[ ! -d "${PACKDIR:?}" ]]; then

    # --------------------------------------------------------------------------
    # Clone from original repo -------------------------------------------------
    # --------------------------------------------------------------------------
    cdfolder "${ECODEF:?}" || { cdroot; return 1; }

    "${CURL:?}" -fsS "${URL:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC27:?} (URL=${URL:?})"; return 1; }

    "${GIT:?}" clone "${URL:?}" "${FOLDER:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC15:?}"; return 1; }
    
    cdfolder "${PACKDIR:?}" || { cdroot; return 1; }

    if [ -n "${ACTDR6_CMBONLY_GIT_COMMIT}" ]; then
      "${GIT:?}" checkout "${ACTDR6_CMBONLY_GIT_COMMIT:?}" \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi

    # --------------------------------------------------------------------------
    # Patch code to be compatible w/ COCOA environment -------------------------
    # --------------------------------------------------------------------------
    # T = TMP
    declare -a TFOLDER=("act_dr6_cmbonly/") # If nonblank, path must include /
    
    # T = TMP
    declare -a TFILE=("act_dr6_cmbonly.py")

    #T = TMP, P = PATCH
    declare -a TFILEP=("act_dr6_cmbonly.patch")
    
    # AL = Array Length
    AL=${#TFOLDER[@]}

    for (( i=0; i<${AL}; i++ ));
    do
    
      cdfolder "${PACKDIR:?}/${TFOLDER[$i]}" || return 1

      cpfolder "${CHANGES:?}/${TFOLDER[$i]}${TFILEP[$i]:?}" . \
        2>${OUT2:?} || return 1;

      patch -u "${TFILE[$i]:?}" -i "${TFILEP[$i]:?}" >${OUT1:?} \
        2>${OUT2:?} || { error "${EC17:?} (${TFILE[$i]:?})"; return 1; }
    
    done

    # need to create symlinks for this likelihood to work w/ COBAYA.
    # we copied the code below to start_cocoa.sh shell script as well
    # we added corresponding code to stop_cocoa.sh that delete these symlinks
    if [[ -d "${PACKDIR:?}/act_dr6_cmbonly" ]]; then
      
      if [[ -L "${COB:?}/${COBLIKE:?}/act_dr6_cmbonly" ]]; then
      
        rm -f "${COB:?}/${COBLIKE:?}/act_dr6_cmbonly"
      
      fi
      
      ln -s "${PACKDIR:?}/act_dr6_cmbonly" "${COB:?}/${COBLIKE:?}" \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC34:?}"; return 1; }

      # it seems this is synthetic data (we will download real data from lambda)
      if [[ -d "${PACKDIR:?}/act_dr6_cmbonly/data" ]]; then
      
        rm -rf "${PACKDIR:?}/act_dr6_cmbonly/data"
      
      fi
    
    fi

  fi

  cdfolder "${ROOTDIR}" || return 1;

  pbottom 'SETUP ACTDR6 (CMB ONLY LIKE)' || return 1

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  ptop 'SETUP ACTDR6 (MFLIKE)' || return 1;

  URL="${ACTDR6_MFLIKE_URL:-"https://github.com/ACTCollaboration/act_dr6_mflike.git"}"

  CHANGES="${CCIL:?}/act_dr6_mflike_changes"

  FOLDER="act_dr6_mflike"

  PACKDIR="${ECODEF:?}/${FOLDER:?}"

  # ---------------------------------------------------------------------------
  # In case this script is called twice ---------------------------------------
  # ---------------------------------------------------------------------------
  if [[ -n "${OVERWRITE_EXISTING_ACTDR6_CMB_CODE}" ]]; then

    rm -rf "${PACKDIR:?}"

  fi

  if [[ ! -d "${PACKDIR:?}" ]]; then
    
    # --------------------------------------------------------------------------
    # Clone from original repo -------------------------------------------------
    # --------------------------------------------------------------------------
    cdfolder "${ECODEF:?}" || { cdroot; return 1; }

    "${CURL:?}" -fsS "${URL:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC27:?} (URL=${URL:?})"; return 1; }

    "${GIT:?}" clone "${URL:?}" "${FOLDER:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC15:?}"; return 1; }
    
    cdfolder "${PACKDIR:?}" || { cdroot; return 1; }

    if [ -n "${ACTDR6_MFLIKE_GIT_COMMIT}" ]; then
      "${GIT:?}" checkout "${ACTDR6_MFLIKE_GIT_COMMIT:?}" \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi

    # --------------------------------------------------------------------------
    # need to create symlinks for this likelihood to work w/ COBAYA.
    # we copied the code below to start_cocoa.sh shell script as well
    # we added corresponding code to stop_cocoa.sh that delete these symlinks
    # --------------------------------------------------------------------------
    if [[ -d "${PACKDIR:?}/act_dr6_mflike" ]]; then
      
      if [[ -L "${COB:?}/${COBLIKE:?}/act_dr6_mflike" ]]; then
      
        rm -f "${COB:?}/${COBLIKE:?}/act_dr6_mflike"
      
      fi
      
      ln -s "${PACKDIR:?}/act_dr6_mflike" "${COB:?}/${COBLIKE:?}" \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC34:?}"; return 1; }
    
    fi

  fi

  cdfolder "${ROOTDIR}" || return 1;
  
  pbottom 'SETUP ACTDR6 (MFLIKE)' || return 1

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  cdfolder "${ROOTDIR}" || return 1;
  
  unset_all || return 1
  
fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------