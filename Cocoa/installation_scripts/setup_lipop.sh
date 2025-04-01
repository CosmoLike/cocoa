#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_LIPOP_LIKELIHOOD_CODE}" ]; then

  if [ -z "${ROOTDIR}" ]; then
    source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
  fi

  # parenthesis = run in a subshell
  ( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

  unset_env_vars () {
    unset -v COB CCCOB COBLIKE URL TFILE TFOLDER ECODEF CCIL
    cdroot || return 1;
  }

  unset_env_funcs () {
    unset -f cdfolder error flipop
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
  
  # E = EXTERNAL, CODE, F=FODLER
  ECODEF="${ROOTDIR:?}/external_modules/code"

  CCIL="${ROOTDIR:?}/../cocoa_installation_libraries" # IL = installation lib

  COBLIKE="${ROOTDIR:?}/cobaya/cobaya/likelihoods" # COB = Cobaya, LIKE = likelihoods
  
  flipop() {
    local TF="${ECODEF:?}/${1:?}"
    
    if [ -n "${OVERWRITE_EXISTING_LIPOP_CMB_DATA}" ]; then
      rm -rf "${TF:?}"
    fi

    if [ ! -d "${TF:?}" ]; then

      cdfolder "${ECODEF:?}" || return 1;

      "${GIT:?}" clone --depth ${GIT_CLONE_MAXIMUM_DEPTH:?} "${3:?}" --recursive \
        ${1:?} >${OUT1:?} 2>${OUT2:?} || { error "${EC15:?}"; return 1; }
    
      cdfolder "${TF:?}" || return 1;

      "${GIT:?}" reset --hard "${4:?}" \
        >${OUT1:?} 2>${OUT2:?} || { error "${EC23:?}"; return 1; }

      cp "${CCIL:?}/${2:?}_changes/init.patch" "${TF:?}/${2:?}" 2>"/dev/null" \
        || { error "CP FILE init.patch on ${TF:?}/${2:?}"; return 1; }

      cdfolder "${TF:?}/${2:?}" || return 1;

      patch -u '__init__.py' -i 'init.patch' >${OUT1:?} 2>${OUT2:?} || 
        { error "${EC17:?}"; return 1; }

      # ------------------------------------------------------------------------
      # Symlinks for likelihood to work w/ COBAYA.(also @start/stop_cocoa.sh)
      # ------------------------------------------------------------------------
      if [[ -d "${TF:?}/${2:?}" ]]; then
        if [[ -L "${COBLIKE:?}/${2:?}" ]]; then
          rm -f "${COBLIKE:?}/${2:?}"
        fi
        ln -s "${TF:?}/${2:?}" "${COBLIKE:?}" \
          >${OUT1:?} 2>${OUT2:?} || { error "${EC34:?}"; return 1; }
      fi
    fi

    cdfolder "${ROOTDIR:?}" || return 1; 
  }
  
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------

  ptop "SETUP HILLIPOP LIKELIHOOD" || return 1

  TFOLDER="${PL2020_HILLIPOP_NAME:-"planck_2020_hillipop"}"
  
  URL="${HILLIPOP_URL:-"https://github.com/planck-npipe/hillipop.git"}"
  
  flipop "${TFOLDER:?}" "planck_2020_hillipop" "${URL:?}" \
    "${HILLIPOP_GIT_COMMIT:-"HEAD~"}" || return 1;
  
  pbottom "SETUP HILLIPOP LIKELIHOOD" || return 1

  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------

  ptop "SETUP LOLLIPOP LIKELIHOOD" || return 1

  TFOLDER="${PL2020_LOLLIPOP_NAME:-"planck_2020_lollipop"}"
  
  URL="${LOLLIPOP_URL:-"https://github.com/planck-npipe/lollipop.git"}"
  
  flipop "${TFOLDER:?}" "planck_2020_lollipop" "${URL:?}" \
    "${LOLLIPOP_GIT_COMMIT:-"HEAD~"}" || return 1;
  
  pbottom "SETUP LOLLIPOP LIKELIHOOD" || return 1

  #-----------------------------------------------------------------------------
  
  cdfolder "${ROOTDIR:?}" || return 1;

  unset_all || return 1;
  
fi

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------