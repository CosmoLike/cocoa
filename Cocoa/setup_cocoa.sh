#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

if [ -n "${ROOTDIR}" ]; then
  source stop_cocoa.sh
fi

error_cip_msg () {
  local FILE="$(basename "${BASH_SOURCE[0]}")"
  local MSG="\033[0;31m\t\t (${FILE}) we cannot run "
  local MSG2="\033[0m"
  echo -e "${MSG}${1:?}${MSG2}" 2>"/dev/null"
}

error_cip () {
  error_cip_msg ${1:?}
  unset -v SCRIPTS
  unset -f error_cip error_cip_msg
  cd $(pwd -P) 2>"/dev/null"
  source stop_cocoa 2>"/dev/null"
  return 1
}

# ------------------------------------------------------------------------------
# ------------------------------ Basic Settings --------------------------------
# ------------------------------------------------------------------------------
source $(pwd -P)/installation_scripts/flags_save_old.sh
if [ $? -ne 0 ]; then
  error_cip 'script flags_save_old.sh' 
  return 1
fi

# note: here is where we define env flag ROOTDIR as $(pwd -P)
source ${SET_INSTALLATION_OPTIONS:-"set_installation_options.sh"}
if [ $? -ne 0 ]; then
  error_cip 'script set_installation_options.sh'; return 1;
fi

if [ -n "${MINICONDA_INSTALLATION}" ]; then
  if [ -z ${CONDA_PREFIX} ]; then
    error_cip "conda environment activation"; return 1;
  fi
fi

ptop2 'SETUP COCOA INSTALLATION PACKAGES'

# ------------------------------------------------------------------------------
# ---------------------- Activate Virtual Environment --------------------------
# ------------------------------------------------------------------------------
if [ -n "${OVERWRITE_EXISTING_COCOA_PRIVATE_PYTHON_ENV}" ]; then
  rm -rf ${ROOTDIR:?}/.local/
fi

if [ ! -d "${ROOTDIR:?}/.local/" ]; then
  ptop 'SETUP COCOA PRIVATE PYTHON ENV'

  cd ${ROOTDIR:?}/../

  if [ -n "${DONT_USE_SYSTEM_PIP_PACKAGES}" ]; then
    ${GLOBALPYTHON3:?} -m venv "${ROOTDIR:?}/.local/"
  else
    ${GLOBALPYTHON3:?} -m venv "${ROOTDIR:?}/.local/" --system-site-packages
  fi

  pbottom 'SETUP COCOA PRIVATE PYTHON ENV'
fi

source "${ROOTDIR:?}/.local/bin/activate"
if [ $? -ne 0 ]; then
  error_cip "cocoa private python environment activation"; return 1;
fi

source "${ROOTDIR:?}/installation_scripts/flags_set_new.sh"
if [ $? -ne 0 ]; then
  error_cip 'script flags_set_new.sh'
  return 1
fi

# ------------------------------------------------------------------------------
# ------------------------------ INSTALL PACKAGES ------------------------------
# ------------------------------------------------------------------------------
declare -i ERRORCODE=0

declare -a SCRIPTS=( "setup_core_packages.sh" 
                     "setup_pip_core_packages.sh"
                     "setup_cobaya.sh"
                     "setup_fgspectra.sh"
                     "setup_simons_observatory.sh"
                     "setup_lipop.sh"
                     "setup_act_dr4.sh"
                     "setup_act_dr6.sh"
                     "setup_polychord.sh"
                     "setup_hyrec2.sh"
                     "setup_cosmorec.sh"
                     "setup_camb.sh"
                     "setup_mgcamb.sh"
                     "setup_class.sh"
                     "setup_velocileptors.sh"
                     "setup_ee2.sh"
                     "setup_cosmopower.sh"
                     "setup_emultrf.sh"
                     "setup_cosmolike_projects.sh"
                     "setup_darkemulator.sh"
                     "setup_private_projects.sh"
                     "setup_nautilus_sampler.sh"
                     "setup_tensiometer.sh"
                     "setup_getdist.sh"
                     "setup_derivkit.sh"
                     "unxv_core_packages.sh"  
                     "unxv_sn.sh"
                     "unxv_bao.sh"
                     "unxv_h0licow.sh" 
                     "unxv_act_dr6.sh"
                     "unxv_simons_observatory.sh"
                     "unxv_bicep.sh"
                     "unxv_spt.sh"
                     "unxv_planck2018_basic.sh"
                     "unxv_camspec.sh"
                     "unxv_lipop.sh"
                     "unxv_emultrf.sh"
                     "unxv_cosmopower.sh")

for (( i=0; i<${#SCRIPTS[@]}; i++ ));
do
  cdroot; 
  ( source "${ROOTDIR:?}/installation_scripts/${SCRIPTS[$i]}" )
  if [ $? -ne 0 ]; then
    if [ -n "${COCOA_OUTPUT_VERBOSE}" ]; then
      error_cip_msg "script ${SCRIPTS[$i]}"
      ERRORCODE=1
    else
      error_cip_msg "script ${SCRIPTS[$i]}"
      ERRORCODE=1
    fi
  fi
done

if [ ${ERRORCODE:?} -ne 0 ]; then
  error_cip "setup_cocoa.py"; return 1
  return 1
fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

unset -v ERRORCODE SCRIPTS
unset -f error_cip error_cip_msg
pbottom2 'SETUP COCOA INSTALLATION PACKAGES'
source stop_cocoa.sh || return 1;

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------