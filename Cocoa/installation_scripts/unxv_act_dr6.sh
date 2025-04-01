#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_ACTDR6_DATA}" ]; then
  
  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi

  # parenthesis = run in a subshell 
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

  unset_env_vars () {
    unset -v EDATAF FOLDER PACKDIR FILE URL_BASE URL PRINTNAME TMP
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

  unset_env_vars || return 1

  # E = EXTERNAL, DATA, F=FODLER
  EDATAF="${ROOTDIR:?}/external_modules/data"

  # --------------------------------------------------------------------------- 
  # ---------------------------------------------------------------------------

  ptop "SETUP/UNXV ACT-DR6 DATA (LENSING)" || return 1

  FOLDER="act"

  # PACK = PACKAGE, DIR = DIRECTORY
  PACKDIR="${EDATAF:?}/${FOLDER:?}"

  FILE="${ACTDR6_LENSING_DATA_FILE:-"ACT_dr6_likelihood_v1.2.tgz"}"

  URL_BASE="https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr6/likelihood/data"

  URL="${ACTDR6_LENSING_DATA_URL:-"${URL_BASE:?}"}/${FILE:?}"

  if [ -n "${OVERWRITE_EXISTING_ACTDR6_CMB_DATA}" ]; then
    rm -rf "${PACKDIR:?}"
    if [ -n "${REDOWNLOAD_EXISTING_ACTDR6_CMB_DATA}" ]; then
      rm -f "${EDATAF:?}/${FILE:?}"
    fi
  fi 

  if [[ ! -d "${PACKDIR:?}" || ! -d "${PACKDIR:?}/lensing" ]]; then
    mkdir -p "${PACKDIR:?}" >${OUT1:?} 2>${OUT2:?} || { error "${EC20:?}"; return 1; }
    mkdir -p "${PACKDIR:?}/lensing" >${OUT1:?} 2>${OUT2:?} || { error "${EC20:?}"; return 1; }
    
    cdfolder "${EDATAF:?}" || return 1

    if [ ! -e "${FILE:?}" ]; then
      "${WGET:?}" "${URL:?}" -q --show-progress --progress=bar:force || 
        { error "${EC24:?}"; return 1; }
    fi

    TMP=$(tar -tf "${FILE:?}" | head -1 | cut -f1 -d"/")
    tar -zxvf "${FILE:?}" >${OUT1:?} 2>${OUT2:?} || { error "${EC25:?}"; return 1; }
    mv "${TMP:?}" "${PACKDIR:?}/lensing"
  fi

  cdfolder "${ROOTDIR}" || return 1;
  
  pbottom "SETUP/UNXV ACT-DR6 DATA (LENSING)" || return 1

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  ptop "SETUP/UNXV ACT-DR6 DATA (CMBONLY)" || return 1

  FOLDER="act_dr6_cmbonly"

  # PACK = PACKAGE, DIR = DIRECTORY
  PACKDIR="${EDATAF:?}/${FOLDER:?}"

  FILE="${ACTDR6_CMBONLY_DATA_FILE:-"dr6_data_cmbonly.tar.gz"}"

  URL_BASE="https://lambda.gsfc.nasa.gov/data/act/pspipe/sacc_files/"

  URL="${ACTDR6_CMBONLY_DATA_URL:-"${URL_BASE:?}"}/${FILE:?}"

  if [ -n "${OVERWRITE_EXISTING_ACTDR6_CMB_DATA}" ]; then    
    rm -rf "${PACKDIR:?}"
    if [ -n "${REDOWNLOAD_EXISTING_ACTDR6_CMB_DATA}" ]; then
      rm -f "${EDATAF:?}/${FILE:?}"
    fi
  fi

  if [ ! -d "${PACKDIR:?}" ]; then
    mkdir -p "${PACKDIR:?}" >${OUT1:?} 2>${OUT2:?} || { error "${EC20:?}"; return 1; }
          
    cdfolder "${EDATAF:?}" || return 1

    if [ ! -e "${FILE:?}" ]; then
      "${WGET:?}" "${URL:?}" -q --show-progress --progress=bar:force || 
        { error "${EC24:?}"; return 1; }
    fi

    TMP=$(tar -tf "${FILE:?}" | head -1 | cut -f1 -d"/")
    tar -zxvf "${FILE:?}" >${OUT1:?} 2>${OUT2:?} || { error "${EC25:?}"; return 1; }
    mv "${TMP:?}" "${PACKDIR:?}"
  fi

  cdfolder "${ROOTDIR}" || return 1;

  pbottom "SETUP/UNXV ACT-DR6 DATA (CMBONLY)" || return 1

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------

  ptop "SETUP/UNXV ACT-DR6 DATA (MFLIKE)" || return 1

  FOLDER="act_dr6_mflike"

  PACKDIR="${EDATAF:?}/${FOLDER:?}"

  FILE="${ACTDR6_MFLIKE_DATA_FILE:-"dr6_data.tar.gz"}"

  URL_BASE="https://lambda.gsfc.nasa.gov/data/act/pspipe/sacc_files/"

  URL="${ACTDR6_MFLIKE_DATA_URL:-"${URL_BASE:?}"}/${FILE:?}"

  if [ -n "${OVERWRITE_EXISTING_ACTDR6_CMB_DATA}" ]; then 
    rm -rf "${PACKDIR:?}"
    if [ -n "${REDOWNLOAD_EXISTING_ACTDR6_CMB_DATA}" ]; then
      rm -f "${EDATAF:?}/${FILE:?}"
    fi
  fi

  if [ ! -d "${PACKDIR:?}" ]; then
    mkdir -p "${PACKDIR:?}" >${OUT1:?} 2>${OUT2:?} || { error "${EC20:?}"; return 1; }
          
    cdfolder "${EDATAF:?}" || return 1

    if [ ! -e "${FILE:?}" ]; then
      "${WGET:?}" "${URL:?}" -q --show-progress \
        --progress=bar:force:noscroll || { error "${EC24:?}"; return 1; }
    fi

    TMP=$(tar -tf "${FILE:?}" | head -1 | cut -f1 -d"/")
    tar -zxvf "${FILE:?}" >${OUT1:?} 2>${OUT2:?} || { error "${EC25:?}"; return 1; }
    mv "${TMP:?}" "${PACKDIR:?}"
  fi

  pbottom "SETUP/UNXV ACT-DR6 DATA (MFLIKE)" || return 1

  cdfolder "${ROOTDIR}" || return 1;

  unset_all || return 1
  
fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------