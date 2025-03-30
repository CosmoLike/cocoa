#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if [ -n "${ROOTDIR}" ]; then
  source stop_cocoa.sh 
fi

error_start_cocoa () {
  local FILE="$(basename ${BASH_SOURCE[0]})"
  local MSG="\033[0;31m (${FILE}) we cannot run "
  local MSG2="\033[0m"
  echo -e "${MSG}${1:?}${MSG2}" 2>"/dev/null"
  unset -f error_cip
  cd $(pwd -P) 2>"/dev/null"
  source stop_cocoa 2>"/dev/null"
  return 1
}

# ----------------------------------------------------------------------------
# ------------------------------ Basic Settings ------------------------------
# ----------------------------------------------------------------------------

source $(pwd -P)/installation_scripts/flags_save_old.sh
if [ $? -ne 0 ]; then
  error_start_cocoa 'script flags_save_old.sh'; return 1;
fi

# note: here is where we define env flag ROOTDIR as $(pwd -P)
source ${SET_INSTALLATION_OPTIONS:-"set_installation_options.sh"}
if [ $? -ne 0 ]; then
  error_start_cocoa 'script set_installation_options.sh'; return 1;
fi

if [ -n "${MINICONDA_INSTALLATION}" ]; then
  if [ -z ${CONDA_PREFIX} ]; then
    error_start_cocoa "conda environment activation"; return 1;
  fi
fi

# ----------------------------------------------------------------------------
# ---------------------- Activate Virtual Environment ------------------------
# ----------------------------------------------------------------------------

source "${ROOTDIR:?}/.local/bin/activate"
if [ $? -ne 0 ]; then
  error_start_cocoa "cocoa private python environment activation"; return 1;
fi

source "${ROOTDIR:?}/installation_scripts/flags_set_new.sh"
if [ $? -ne 0 ]; then
  error_start_cocoa 'script flags_set_new.sh'; return 1;
fi

# ----------------------------------------------------------------------------
# ----------------------------- PLANCK LIKELIHOOD ----------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_PLANCK_COMPILATION}" ]; then
  
  if [ -n "${CLIK_PATH}" ]; then
    export OLD_CLIK_PATH=$CLIK_PATH
  else
    export OLD_CLIK_PATH="x"
  fi

  if [ -n "${CLIK_DATA}" ]; then
    export OLD_CLIK_DATA=$CLIK_DATA
  else
    export OLD_CLIK_DATA="x"
  fi

  if [ -n "${CLIK_PLUGIN}" ]; then
    export OLD_CLIK_PLUGIN=$CLIK_PLUGIN
  else
    export OLD_CLIK_PLUGIN="x"
  fi
  
  export CLIK_PATH="${ROOTDIR:?}/.local"

  export CLIK_DATA="${ROOTDIR:?}/.local/share/clik"

  export CLIK_PLUGIN=rel2015

fi

# ----------------------------------------------------------------------------
# ------------------------------- ACT LIKELIHOOD -----------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_ACTDR6_CODE}" ]; then

  if [[ ! -L "${ROOTDIR:?}/cobaya/cobaya/likelihoods/act_dr6_cmbonly" ]]; then
    ln -s "${ROOTDIR:?}/external_modules/code/act_dr6_cmbonly/act_dr6_cmbonly" \
          "${ROOTDIR:?}/cobaya/cobaya/likelihoods" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC34:?}"; return 1; }
  fi

  if [[ ! -L "${ROOTDIR:?}/cobaya/cobaya/likelihoods/act_dr6_mflike" ]]; then
    ln -s "${ROOTDIR:?}/external_modules/code/act_dr6_mflike/act_dr6_mflike" \
          "${ROOTDIR:?}/cobaya/cobaya/likelihoods" \
      >${OUT1:?} 2>${OUT2:?} || { error "${EC34:?}"; return 1; }
  fi

fi

# ----------------------------------------------------------------------------
# ------------------------ START EXTERNAL PROJECTS ---------------------------
# ----------------------------------------------------------------------------

source "${ROOTDIR:?}/installation_scripts/start_all_projects.sh"

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

unset -f error_start_cocoa 

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------