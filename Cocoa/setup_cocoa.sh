#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------ Basic Settings --------------------------------
# ------------------------------------------------------------------------------
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  FILE="$(basename "${BASH_SOURCE[0]}")"
  MSG="\033[0;31m ${FILE} must be sourced (not executed as program)"
  MSG2=", e.g.: \n source ${FILE}\033[0m"
  echo -e "${MSG}${MSG2}"
  unset FILE MSG MSG2
  exit 1
fi

if [ -n "${ROOTDIR}" ]; then
  source stop_cocoa.sh
fi

unset_all_cip () {
  unset -v ERRORCODE SCRIPTS CACHE CACHE_FILE mode_arg_seen mode INIT END
  unset -v CORE THEORY ML LIKELIHOOD DATA COSMOLIKE 
  unset -f init_cache reset_cache save_cache load_cache
  unset -f error_cip error_cip_msg unset_all_cip
}

error_cip_msg () {
  local FILE="$(basename "${BASH_SOURCE[0]}")"
  local MSG="\033[0;31m\t\t (${FILE}) we cannot run "
  local MSG2="\033[0m"
  echo -e "${MSG}${1:?}${MSG2}" 2>"/dev/null"
}

error_cip () {
  error_cip_msg "${1:?}"
  unset_all_cip
  cd "$(pwd -P)" 2>"/dev/null"
  source stop_cocoa 2>"/dev/null"
  return 1
}

source "$(pwd -P)/installation_scripts/flags_save_old.sh"
if [ $? -ne 0 ]; then
  error_cip 'script flags_save_old.sh' 
  return 1
fi

# note: here is where we define env flag ROOTDIR as $(pwd -P)
source "${SET_INSTALLATION_OPTIONS:-"set_installation_options.sh"}"
if [ $? -ne 0 ]; then
  error_cip 'script set_installation_options.sh'; return 1;
fi

if [ -n "${MINICONDA_INSTALLATION}" ]; then
  if [ -z "${CONDA_PREFIX:-}" ]; then
    error_cip "conda environment activation"; return 1;
  fi
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# --------------------------- SET PACKAGES TO INSTALL --------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
declare -a CORE=("setup_core_packages.sh" 
                 "setup_pip_core_packages.sh"
                 "setup_cobaya.sh"
                 "setup_polychord.sh"
                 "setup_nautilus_sampler.sh"
                 "setup_tensiometer.sh"
                 "setup_getdist.sh"
                 "setup_derivkit.sh"
                 "setup_ee2.sh"
                 "setup_simde.sh"
                )

declare -a THEORY=("setup_hyrec2.sh"
                   "setup_cosmorec.sh"
                   "setup_camb.sh"
                   "setup_mgcamb.sh"
                   "setup_class.sh"
                   "setup_velocileptors.sh"
                   "setup_pyfastpt.sh"
                  )

declare -a ML=("setup_cosmopower.sh"
               "setup_emultrf.sh"
               "setup_darkemulator.sh"
              )

declare -a LIKELIHOOD=("setup_fgspectra.sh"
                       "setup_simons_observatory.sh"
                       "setup_lipop.sh"
                       "setup_act_dr4.sh"
                       "setup_act_dr6.sh"
                      )

declare -a DATA=("unxv_core_packages.sh"  
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
                 "unxv_cosmopower.sh"
                )

declare -a COSMOLIKE=("setup_cosmolike.sh"
                      "setup_private_projects.sh"
                      "setup_cosmolike_projects.sh"
                     )
declare -a SCRIPTS=()
SCRIPTS+=("${CORE[@]}")
SCRIPTS+=("${THEORY[@]}")
SCRIPTS+=("${ML[@]}")
SCRIPTS+=("${LIKELIHOOD[@]}")
SCRIPTS+=("${DATA[@]}")
SCRIPTS+=("${COSMOLIKE[@]}")

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------ 
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------                    
# Below here should not require user input
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------ 
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------ 

# ------------------------------------------------------------------------------
# ------------------------------ CHOOSE RUN MODES ------------------------------
# ------------------------------------------------------------------------------
CACHE_FILE="${ROOTDIR:?}/.setup_local_packages.txt"
mode="default"        # default
mode_arg_seen=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --soft|--hard|--aggressive|--extreme|--purge)
      if [[ $mode_arg_seen -eq 1 ]]; then
        error_cip "choose only one mode"
        return 1
      fi
      mode_arg_seen=1
      mode="${1#--}"
      ;;
    *)
      error_cip "unknown run mode option $1"
      return 1
      ;;
  esac
  shift
done

# ------------------------------------------------------------------------------
# -------------------------- SET CACHE  ----------------------------------------
# ------------------------------------------------------------------------------

declare -a CACHE=()

load_cache() {
  local file="$1"
  local line
  CACHE=()
  # read line and continue if line is not empty --------------------------------
  while IFS= read -r line || [[ -n "$line" ]]; do
    CACHE+=("$line")
  done < "$file"
}

save_cache() {
  printf '%s\n' "${CACHE[@]}" > "${CACHE_FILE:?}"
}

reset_cache() {
  local i
  CACHE=()
  for (( i=0; i<${#SCRIPTS[@]}; i++ )); do
    CACHE[i]=0
  done
  save_cache
}

init_cache() {
  local i
  local LOCAL_DIR="${ROOTDIR:?}"
  if [[ ! -d "${LOCAL_DIR:?}" ]]; then
    error_cip "(.local)  does not exist: ${LOCAL_DIR:?}"
    return 1 
  fi
  if [[ -f "${CACHE_FILE:?}" ]]; then
    load_cache "${CACHE_FILE:?}"
    # If the file size does not match the number of scripts, rebuild it --------
    if [[ ${#CACHE[@]} -ne ${#SCRIPTS[@]} ]]; then
      reset_cache
      return
    fi
    # Validate entries: only 0 or 1 are allowed --------------------------------
    for (( i=0; i<${#CACHE[@]}; i++ )); do
      if [[ "${CACHE[i]}" != "0" && "${CACHE[i]}" != "1" ]]; then
        reset_cache
        return
      fi
    done
    return
  else
    reset_cache
  fi
}

init_cache

# ------------------------------------------------------------------------------
# -------------------------- CONFIG MODE ---------------------------------------
# ------------------------------------------------------------------------------

case "$mode" in
  default)
    source "${ROOTDIR:?}/installation_scripts/flags_derived.sh"
    if [ $? -ne 0 ]; then
      return 1;
    fi
    unset -v OVERWRITE_EXISTING_COCOA_PRIVATE_PYTHON_ENV
    unset -v OVERWRITE_EXISTING_COSMOLIKE_CODE # ADDITIONAL PROTECTION
    unset -v OVERWRITE_EXISTING_PRIVATE_CODE   # ADDITIONAL PROTECTION
    ;;
  soft)
    # FORCE DOWNLOAD OF THEORY + ML
    INIT=${#CORE[@]}
    END=$(( ${#CORE[@]} + ${#THEORY[@]} + ${#ML[@]} ))
    for (( i=$INIT; i<$END; i++ ));
    do
      CACHE[i]=0
    done
    source "${ROOTDIR:?}/installation_scripts/flags_derived.sh"
    if [ $? -ne 0 ]; then
      return 1;
    fi
    unset -v OVERWRITE_EXISTING_COCOA_PRIVATE_PYTHON_ENV
    unset -v OVERWRITE_EXISTING_COSMOLIKE_CODE # ADDITIONAL PROTECTION
    unset -v OVERWRITE_EXISTING_PRIVATE_CODE   # ADDITIONAL PROTECTION
    ;;
  hard)
    # FORCE DOWNLOAD OF THEORY + ML + LIKELIHOOD
    INIT=${#CORE[@]}
    END=$(( ${#THEORY[@]} + ${#CORE[@]} + ${#ML[@]} + ${#LIKELIHOOD[@]} ))
    for (( i=$INIT; i<$END; i++ ));
    do
      CACHE[i]=0
    done
    source "${ROOTDIR:?}/installation_scripts/flags_derived.sh"
    if [ $? -ne 0 ]; then
      return 1;
    fi
    unset -v OVERWRITE_EXISTING_COCOA_PRIVATE_PYTHON_ENV
    unset -v OVERWRITE_EXISTING_COSMOLIKE_CODE # ADDITIONAL PROTECTION
    unset -v OVERWRITE_EXISTING_PRIVATE_CODE   # ADDITIONAL PROTECTION
    ;;
  aggressive)
    # FORCE DOWNLOAD OF CORE + THEORY + ML + LIKELIHOOD
    INIT=0
    END=$(( ${#THEORY[@]} + ${#CORE[@]} + ${#ML[@]} + ${#LIKELIHOOD[@]} ))
    for (( i=$INIT; i<$END; i++ ));
    do
      CACHE[i]=0
    done
    source "${ROOTDIR:?}/installation_scripts/flags_derived.sh"
    if [ $? -ne 0 ]; then
      return 1;
    fi
    unset -v OVERWRITE_EXISTING_COCOA_PRIVATE_PYTHON_ENV
    unset -v OVERWRITE_EXISTING_COSMOLIKE_CODE # ADDITIONAL PROTECTION
    unset -v OVERWRITE_EXISTING_PRIVATE_CODE   # ADDITIONAL PROTECTION
    ;;
  extreme)
    # FORCE DOWNLOAD OFF ALL CODE (EXCEPT COSMOLIKE PROJECTS)
    INIT=0
    END=$(( ${#THEORY[@]} + ${#CORE[@]} + ${#ML[@]} + ${#DATA[@]} ))
    for (( i=$INIT; i<$END; i++ ));
    do
      CACHE[i]=0
    done
    source "${ROOTDIR:?}/installation_scripts/flags_derived.sh"
    if [ $? -ne 0 ]; then
      return 1;
    fi
    export OVERWRITE_EXISTING_COCOA_PRIVATE_PYTHON_ENV=1
    unset -v OVERWRITE_EXISTING_COSMOLIKE_CODE # ADDITIONAL PROTECTION
    unset -v OVERWRITE_EXISTING_PRIVATE_CODE   # ADDITIONAL PROTECTION
    ;;
  purge)
    # FORCE DOWNLOAD OFF ALL CODE (INCLUDING COSMOLIKE PROJECTS)
    INIT=0
    END=${#SCRIPTS[@]}
    for (( i=$INIT; i<$END; i++ ));
    do
      CACHE[i]=0
    done
    source "${ROOTDIR:?}/installation_scripts/flags_derived.sh"
    if [ $? -ne 0 ]; then
      return 1;
    fi
    export OVERWRITE_EXISTING_COCOA_PRIVATE_PYTHON_ENV=1
    export OVERWRITE_EXISTING_COSMOLIKE_CODE=1
    export OVERWRITE_EXISTING_PRIVATE_CODE=1
    ;;
  *)
    echo "Error: invalid mode: $mode" >&2
    exit 1
    ;;
esac

unset mode

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

ptop2 'SETUP COCOA & EXTERNAL MODULES' || { unset_all; return 1; }

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
# -------------------------- RUN INSTALL PACKAGES ------------------------------
# ------------------------------------------------------------------------------

declare -i ERRORCODE=0

for (( i=0; i<${#SCRIPTS[@]}; i++ ));
do
  cdroot; 

  [[ "${CACHE[i]}" -eq 1 ]] && continue
  
  if [ -n "${COCOA_OUTPUT_DEBUG}" ]; then
    # bash strict mode explanation
    # http://redsymbol.net/articles/unofficial-bash-strict-mode/
    ( set -exo pipefail; source "${ROOTDIR:?}/installation_scripts/${SCRIPTS[$i]}" )
  else
    ( source "${ROOTDIR:?}/installation_scripts/${SCRIPTS[$i]}" )
  fi
  
  rc=$?

  if [[ ${rc:?} -eq 55 ]]; then
    CACHE[i]=1
    save_cache
  elif [[ ${rc:?} -eq 99 ]]; then
    continue
  else
    if [ -n "${COCOA_OUTPUT_VERBOSE}" ]; then
      error_cip "script ${SCRIPTS[$i]}"
      return 1
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


unset_all_cip
pbottom2 'SETUP EXTERNAL & EXTERNAL MODULES' || return 1;
source stop_cocoa.sh || return 1;

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------