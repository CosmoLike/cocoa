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
  echo -e "${MSG}${1:?}${MSG2}"
  unset -f error_start_cocoa
  cd $(pwd -P) 
  source stop_cocoa.sh
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
if [[ -z "${IGNORE_ACTDR6_CODE}" ]]; then
  ECODEF="${ROOTDIR:?}/external_modules/code"
  COBLIKE="${ROOTDIR:?}/cobaya/cobaya/likelihoods"

  TMP="${ACTDR6_CMBONLY_NAME:-"act_dr6_cmbonly"}"
  TMP2="act_dr6_cmbonly"
  if [[ ! -L "${COBLIKE:?}/${TMP2}" ]]; then
    ln -s "${ECODEF:?}/${TMP}/${TMP2}" "${COBLIKE:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error_start_cocoa "${EC34:?}"; return 1; }
  fi

  TMP="${ACTDR6_MFLIKE_NAME:-"act_dr6_mflike"}"
  TMP2="act_dr6_mflike"
  if [[ ! -L "${COBLIKE:?}/${TMP2}" ]]; then
    ln -s "${ECODEF:?}/${TMP}/${TMP2}" "${COBLIKE:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error_start_cocoa "${EC34:?}"; return 1; }
  fi

  unset -v ECODEF COBLIKE TMP TMP2
fi

# ----------------------------------------------------------------------------
# ------------------------------- SO LIKELIHOOD ------------------------------
# ----------------------------------------------------------------------------
if [[ -z "${IGNORE_SIMONS_OBSERVATORY_LIKELIHOOD_CODE}" ]]; then
  ECODEF="${ROOTDIR:?}/external_modules/code"
  COBLIKE="${ROOTDIR:?}/cobaya/cobaya/likelihoods"
  COBTH="${ROOTDIR:?}/cobaya/cobaya/theories"

  TMP="${SO_MFLIKE_NAME:-"mflike"}"
  TMP2="mflike"
  
  if [[ ! -L "${COBLIKE:?}/${TMP2}" ]]; then
    ln -s "${ECODEF:?}/${TMP}/${TMP2}" "${COBLIKE:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error_start_cocoa "${EC34:?}"; return 1; }
  fi
  if [[ ! -L "${COBTH:?}/${TMP2}" ]]; then
    ln -s "${ECODEF:?}/${TMP}/${TMP2}" "${COBTH:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error_start_cocoa "${EC34:?}"; return 1; }
  fi

  unset -v ECODEF COBLIKE TMP TMP2 COBTH
fi

# ----------------------------------------------------------------------------
# ----------------------------- LIPOP LIKELIHOOD -----------------------------
# ----------------------------------------------------------------------------
if [[ -z "${IGNORE_LIPOP_LIKELIHOOD_CODE}" ]]; then
  ECODEF="${ROOTDIR:?}/external_modules/code"
  COBLIKE="${ROOTDIR:?}/cobaya/cobaya/likelihoods"

  TMP="${PL2020_HILLIPOP_NAME:-"planck_2020_hillipop"}"
  TMP2="planck_2020_hillipop"
  if [[ ! -L "${COBLIKE:?}/${TMP2}" ]]; then
    ln -s "${ECODEF:?}/${TMP}/${TMP2}" "${COBLIKE:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error_start_cocoa "${EC34:?}"; return 1; }
  fi

  TMP="${PL2020_LOLLIPOP_NAME:-"planck_2020_lollipop"}"
  TMP2="planck_2020_lollipop"
  if [[ ! -L "${COBLIKE:?}/${TMP2}" ]]; then
    ln -s "${ECODEF:?}/${TMP}/${TMP2}" "${COBLIKE:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error_start_cocoa "${EC34:?}"; return 1; }
  fi

  unset -v ECODEF COBLIKE TMP TMP2
fi

# ----------------------------------------------------------------------------
# ----------------------------- COSMOPOWER THEORY ----------------------------
# ----------------------------------------------------------------------------
if [[ -z "${IGNORE_COSMOPOWER_CODE}" ]]; then
  ECODEF="${ROOTDIR:?}/external_modules/code"
  COBTH="${ROOTDIR:?}/cobaya/cobaya/theories"

  TMP="${COSMOPOWER_SOLIKET_NAME:-"soliket"}"
  TMP2="soliket/cosmopower"
  TMP3="cosmopower"

  if [[ ! -L "${COBTH:?}/${TMP3:?}" ]]; then
    ln -s "${ECODEF:?}/emulators/${TMP}/${TMP2}" "${COBTH:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error_start_cocoa "${EC34:?}"; return 1; }
  fi

  unset -v ECODEF COBTH TMP TMP2 TMP3
fi

# ----------------------------------------------------------------------------
# ----------------------------- EMUL TRF THEORY -------------------------------
# ----------------------------------------------------------------------------
if [[ -z "${IGNORE_EMULTRF_CODE}" ]]; then
  ECODEF="${ROOTDIR:?}/external_modules/code"
  COBTH="${ROOTDIR:?}/cobaya/cobaya/theories"
  TMP="${EMULTRF_NAME:-"emultrf"}"

  TMP2="emulcmb"
  if [[ ! -L "${COBTH:?}/${TMP2}" ]]; then
    ln -s "${ECODEF:?}/emulators/${TMP}/${TMP2}" "${COBTH:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error_start_cocoa "${EC34:?}"; return 1; }
  fi

  TMP2="emulbaosn"
  if [[ ! -L "${COBTH:?}/${TMP2}" ]]; then
    ln -s "${ECODEF:?}/emulators/${TMP}/${TMP2}" "${COBTH:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error_start_cocoa "${EC34:?}"; return 1; }
  fi

  TMP2="emultheta"
  if [[ ! -L "${COBTH:?}/${TMP2}" ]]; then
    ln -s "${ECODEF:?}/emulators/${TMP}/${TMP2}" "${COBTH:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error_start_cocoa "${EC34:?}"; return 1; }
  fi

  TMP2="emulrdrag"
  if [[ ! -L "${COBTH:?}/${TMP2}" ]]; then
    ln -s "${ECODEF:?}/emulators/${TMP}/${TMP2}" "${COBTH:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error_start_cocoa "${EC34:?}"; return 1; }
  fi

  TMP2="emul_cosmic_shear"
  if [[ ! -L "${COBTH:?}/${TMP2}" ]]; then
    ln -s "${ECODEF:?}/emulators/${TMP}/${TMP2}" "${COBTH:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error_start_cocoa "${EC34:?}"; return 1; }
  fi

  TMP2="emul_ggl"
  if [[ ! -L "${COBTH:?}/${TMP2}" ]]; then
    ln -s "${ECODEF:?}/emulators/${TMP}/${TMP2}" "${COBTH:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error_start_cocoa "${EC34:?}"; return 1; }
  fi

  TMP2="emul_wtheta"
  if [[ ! -L "${COBTH:?}/${TMP2}" ]]; then
    ln -s "${ECODEF:?}/emulators/${TMP}/${TMP2}" "${COBTH:?}" \
      >${OUT1:?} 2>${OUT2:?} || { error_start_cocoa "${EC34:?}"; return 1; }
  fi

  unset -v ECODEF COBTH TMP TMP2
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

unset -f error_start_cocoa 

# ----------------------------------------------------------------------------
# ------------------------ START EXTERNAL PROJECTS ---------------------------
# ----------------------------------------------------------------------------

source "${ROOTDIR:?}/installation_scripts/start_all_projects.sh"

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------