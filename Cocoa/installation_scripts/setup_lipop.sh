#!/bin/bash
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
if [ -z "${IGNORE_LIPOP_LIKELIHOOD_COMPILATION}" ]; then

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
    unset -f cdfolder cpfolder error cpfile flipop cppatch cppatchfolder
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

  flipop() {
    local TFOLDER="${COBLIKE:?}/${1:?}"
    local URL="${2:?}"

    # in case you run this script more than once
    rm -rf "${COB:?}/${TFOLDER:?}"
    rm -rf "${COB:?}/${COBLIKE:?}/tmp"

    cdfolder "${COB:?}/${COBLIKE:?}" || return 1;

    "${CURL:?}" -fsS "${URL:?}" >${OUT1:?} \
      2>${OUT2:?} || { error "${EC27:?} (URL=${URL:?})"; return 1; }

    "${GIT:?}" clone "${URL:?}" "tmp" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC15:?}"; return 1; }
  
    cdfolder "${COB:?}/${COBLIKE}/tmp" || return 1;

    "${GIT:?}" reset --hard "${3:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC23:?}"; return 1; }

    cpfolder ${1:?}/ "${COB:?}/${COBLIKE:?}" || return 1;
    
    cdfolder "${COB:?}/${TFOLDER:?}" || return 1;
    
    rm -rf "${COB:?}/${COBLIKE:?}/tmp"

    cppatch "${TFOLDER:?}" "init.patch" || return 1

    patch -u '__init__.py' -i 'init.patch' >${OUT1:?} 2>${OUT2:?} || 
      { error "${EC17:?}"; return 1; }

    cdfolder "${ROOTDIR:?}" || return 1; 
  }
  
  #---------------------------------------------------------------------------
  # HILLIPOP_URL LIKELIHOOD --------------------------------------------------
  #---------------------------------------------------------------------------

  ptop "SETUP HILLIPOP LIKELIHOOD" || return 1

  TFOLDER="planck_2020_hillipop"
  
  URL="${HILLIPOP_URL:-"https://github.com/planck-npipe/hillipop.git"}"
  
  flipop "${TFOLDER:?}" "${URL:?}" "${HILLIPOP_GIT_COMMIT:-"HEAD~"}" || return 1;
  
  unset -v TFOLDER URL

  pbottom "SETUP HILLIPOP LIKELIHOOD" || return 1

  #---------------------------------------------------------------------------
  # LOLLIPOP LIKELIHOOD ------------------------------------------------------
  #---------------------------------------------------------------------------

  ptop "SETUP LOLLIPOP LIKELIHOOD" || return 1

  TFOLDER="planck_2020_lollipop"
  
  URL="${LOLLIPOP_URL:-"https://github.com/planck-npipe/lollipop.git"}"
  
  flipop "${TFOLDER:?}" "${URL:?}" "${LOLLIPOP_GIT_COMMIT:-"HEAD~"}" || return 1;
  
  unset -v TFOLDER URL

  pbottom "SETUP LOLLIPOP LIKELIHOOD" || return 1

  #-----------------------------------------------------------------------------
  
  unset_all || return 1;
  
fi

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------