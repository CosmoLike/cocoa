#!/bin/bash
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
if [ -z "${IGNORE_CAMSPEC_LIKELIHOOD_COMPILATION}" ]; then

  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi

  # parenthesis = run in a subshell
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

  unset_env_vars () {
    unset -v COB CCCOB COBLIKE URL TFILE TFOLDER
    cdroot || return 1;
  }

  unset_env_funcs () {
    unset -f cdfolder cpfolder error cpfile cppatch cppatchfolder
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

  cpfile() {
    cp "${1:?}" "${2:?}" \
      2>"/dev/null" || { error "CP FILE ${1} on ${2}"; return 1; }
  }

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  unset_env_vars || return 1

  CCIL="${ROOTDIR:?}/../cocoa_installation_libraries" # IL = installation lib

  COB="${ROOTDIR:?}/cobaya"        # COB = Cobaya

  CCCOB="${CCIL:?}/cobaya_changes"  # CC = CoCoA, COB = Cobaya (Cocoa Cobaya)
  
  COBLIKE="cobaya/likelihoods"      # COB = Cobaya, LIKE = likelihoods
  
  cppatch() {
    cp "${CCCOB:?}/${1:?}/${2:?}" "${COB:?}/${1:?}" \
      2>"/dev/null" || 
      { error "CP FILE ${CCCOB:?}/${1:?}/${2:?} on ${COB:?}/${1:?}"; return 1; }
  }

  cppatchfolder() {
    cp -r "${CCCOB:?}/${1:?}/${2:?}" "${COB:?}/${1:?}" \
      2>"/dev/null" || 
      { error "CP FOLDER ${CCCOB:?}/${1:?}/${2:?} on ${COB:?}/${1:?}"; \
      return 1; }
  }

  #-----------------------------------------------------------------------------
  # INSTALL CAMSPEC LIKELIHOOD -------------------------------------------------
  #-----------------------------------------------------------------------------

  ptop "SETUP CAMSPEC-2021 LIKELIHOOD" || return 1;

  # note: we always delete CAMSPEC2018 likelihood
  rm -rf "${COB:?}/${COBLIKE:?}"/planck_2018_highl_CamSpec

  TFOLDER="${COBLIKE}/base_classes"
  
  TFILE="InstallableLikelihood"
  
  cppatch "${TFOLDER:?}" "${TFILE:?}.patch"
  
  cdfolder "${COB:?}/${TFOLDER}/" || return 1;

  patch -u "${TFILE:?}.py" -i "${TFILE:?}.patch" \
    >${OUT1:?} 2>${OUT2:?} || { error "${EC17:?}"; return 1; }
  
  unset -v TFOLDER TFILE
  
  cdfolder "${ROOTDIR:?}" || return 1;

  pbottom "SETUP CAMSPEC-2021 LIKELIHOOD" || return 1;

  #-----------------------------------------------------------------------------
  
  unset_all || return 1;
  
fi

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------