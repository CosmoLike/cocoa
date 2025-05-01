#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${ROOTDIR}" ]; then
  ROOTDIR=$(pwd -P) || { echo -e \
  "\033[0;31m       ERROR ENV VARIABLE ROOTDIR NOT DEFINED \033[0m"; return 1; }
fi
  
if command -v deactivate &> "/dev/null"; then
  deactivate 2>"/dev/null"
fi

error_stop_cocoa () {
  local FILE="$(basename ${BASH_SOURCE[0]})"
  local MSG="\033[0;31m (${FILE}) we cannot run "
  local MSG2="\033[0m"
  echo -e "${MSG}${1:?}${MSG2}" 2>"/dev/null"
  unset -f error_stop_cocoa
  cd $(pwd -P) 2>"/dev/null"
  source "${ROOTDIR:?}/installation_scripts/flags_impl_unset_keys.sh" 2>"/dev/null"
  unset -v ROOTDIR SETUP_COBAYA START_COCOA_DONE fail
  unset -v SETUP_PREREQUISITE_DONE SET_INSTALLATION_OPTIONS
  return 1
}

source "$(pwd -P)/installation_scripts/flags_recover_old.sh"
if [ $? -ne 0 ]; then
  error_stop_cocoa 'script installation_scripts/flags_recover_old.sh'; 
  return 1;
fi

# ----------------------------------------------------------------------------
# ----------------------------- PLANCK LIKELIHOOD ----------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_PLANCK_COMPILATION}" ]; then
  if [ -n "${OLD_CLIK_PATH}" ]; then
    if [ "${OLD_CLIK_PATH}" != "x" ]; then
      export CLIK_PATH=$OLD_CLIK_PATH
    else
      unset CLIK_PATH
    fi
    unset OLD_CLIK_PATH
  fi

  if [ -n "${OLD_CLIK_DATA}" ]; then
    if [ "${OLD_CLIK_DATA}" != "x" ]; then
      export CLIK_DATA=$OLD_CLIK_DATA
    else
      unset CLIK_DATA
    fi
    unset OLD_CLIK_DATA
  fi

  if [ -n "${OLD_CLIK_PLUGIN}" ]; then
    if [ "${OLD_CLIK_PLUGIN}" != "x" ]; then
      export CLIK_PLUGIN=$OLD_CLIK_PLUGIN
    else
      unset CLIK_PLUGIN
    fi
    unset OLD_CLIK_PLUGIN
  fi
fi

# ----------------------------------------------------------------------------
# ------------------------------- ACT LIKELIHOOD -----------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_ACTDR6_CODE}" ]; then
  COBLIKE="${ROOTDIR:?}/cobaya/cobaya/likelihoods"

  TMP="act_dr6_cmbonly"
  if [[ -L "${COBLIKE:?}/${TMP:?}" ]]; then
    rm -f "${COBLIKE:?}/${TMP:?}"
  fi

  TMP="act_dr6_mflike"
  if [[ -L "${COBLIKE:?}/${TMP:?}" ]]; then
    rm -f "${COBLIKE:?}/${TMP:?}"
  fi

  unset -v COBLIKE TMP
fi

# ----------------------------------------------------------------------------
# ------------------------------- SO LIKELIHOOD ------------------------------
# ----------------------------------------------------------------------------
if [[ -z "${IGNORE_SIMONS_OBSERVATORY_LIKELIHOOD_CODE}" ]]; then
  COBLIKE="${ROOTDIR:?}/cobaya/cobaya/likelihoods"
  COBTH="${ROOTDIR:?}/cobaya/cobaya/theories"
  TMP="mflike"
  
  if [[ -L "${COBLIKE:?}/${TMP:?}" ]]; then
    rm -f "${COBLIKE:?}/${TMP:?}"
  fi
  if [[ -L "${COBTH:?}/${TMP:?}" ]]; then
    rm -f "${COBTH:?}/${TMP:?}"
  fi

  unset -v COBLIKE COBTH TMP
fi

# ----------------------------------------------------------------------------
# ----------------------------- LIPOP LIKELIHOOD -----------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_ACTDR6_CODE}" ]; then
  COBLIKE="${ROOTDIR:?}/cobaya/cobaya/likelihoods"

  TMP="planck_2020_hillipop"
  if [[ -L "${COBLIKE:?}/${TMP:?}" ]]; then
    rm -f "${COBLIKE:?}/${TMP:?}"
  fi

  TMP="planck_2020_lollipop"
  if [[ -L "${COBLIKE:?}/${TMP:?}" ]]; then
    rm -f "${COBLIKE:?}/${TMP:?}"
  fi

  unset -v COBLIKE TMP
fi

# ----------------------------------------------------------------------------
# ----------------------------- COSMOPOWER THEORY ----------------------------
# ----------------------------------------------------------------------------
if [[ -z "${IGNORE_COSMOPOWER_CODE}" ]]; then
  COBTH="${ROOTDIR:?}/cobaya/cobaya/theories"

  TMP="cosmopower"
  if [[ -L "${COBTH:?}/${TMP:?}" ]]; then
    rm -f "${COBTH:?}/${TMP:?}"
  fi

  unset -v COBTH TMP
fi

# ----------------------------------------------------------------------------
# ----------------------------- CMB TRF THERY --------------------------------
# ----------------------------------------------------------------------------
if [[ -z "${IGNORE_EMULTRF_CODE}" ]]; then
  COBTH="${ROOTDIR:?}/cobaya/cobaya/theories"

  TMP="emulcmbtrf"
  if [[ -L "${COBTH:?}/${TMP:?}" ]]; then
    rm -f "${COBTH:?}/${TMP:?}"
  fi

  unset -v COBTH TMP
fi

# ----------------------------------------------------------------------------
# ------------------------ STOP EXTERNAL PROJECTS ---------------------------
# ----------------------------------------------------------------------------


if [ -n "${ROOTDIR}" ]; then
  source "${ROOTDIR:?}/installation_scripts/stop_all_projects.sh"
  if [ $? -ne 0 ]; then
    error_stop_cocoa 'script stop_all_projects.sh'; return 1
  fi
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

source "${ROOTDIR:?}/installation_scripts/flags_impl_unset_keys.sh" 2>"/dev/null"
unset -v ROOTDIR SETUP_COBAYA START_COCOA_DONE fail
unset -v SETUP_PREREQUISITE_DONE SET_INSTALLATION_OPTIONS

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------