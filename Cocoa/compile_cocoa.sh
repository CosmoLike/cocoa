#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------ Basic Settings --------------------------------
# ------------------------------------------------------------------------------

if [[ ! "${BASH_SOURCE[0]}" != "$0" ]]; then
  FILE="$(basename ${BASH_SOURCE[0]})"
  MSG="\033[0;31m ${FILE} must be sourced (not executed as program)"
  MSG2=", e.g.: \n source ${FILE}\033[0m"
  echo -e "${MSG}${MSG2}"
  unset FILE MSG MSG2
  exit 1
fi

if [ -n "${ROOTDIR}" ]; then
  source stop_cocoa.sh
fi

source start_cocoa.sh || { error_cem "start_cocoa.sh"; return 1; }

# ------------------------------------------------------------------------------
# ----------------- LIST OF SCRIPTS TO COMPILE ---------------------------------
# ------------------------------------------------------------------------------

declare -a CORE=("compile_core_packages.sh"
                 "compile_getdist.sh"
                 "compile_derivkit.sh"
                 "compile_polychord.sh"
                 "compile_nautilus_sampler.sh"
                 "compile_tensiometer.sh"
                 "compile_ee2.sh"
                )

declare -a THEORY=("compile_hyrec2.sh"
                   "compile_cosmorec.sh"
                   "compile_camb.sh"
                   "compile_class.sh"
                   "compile_mgcamb.sh"
                   "compile_velocileptors.sh" 
                  )

declare -a LIKELIHOOD=("compile_planck.sh"
                       "compile_fgspectra.sh"
                       "compile_simons_observatory.sh"
                       "compile_act_dr4.sh"
                       "compile_act_dr6.sh"
                       
                       "compile_lipop.sh"

                    )

declare -a ML=("compile_cosmopower.sh"
               "compile_darkemulator.sh"
              )

declare -a SCRIPTS=()
SCRIPTS+=("${CORE[@]}")
SCRIPTS+=("${THEORY[@]}")
SCRIPTS+=("${ML[@]}")
SCRIPTS+=("${LIKELIHOOD[@]}")
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------ 
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------                    
# Below here should not require users input
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------ 
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------ 

unset_all () {
  unset -v ERRORCODE SCRIPTS ML LIKELIHOOD THEORY CORE 
  unset -v CACHE CACHE_FILE mode_arg_seen mode
  unset -f init_cache reset_cache save_cache load_cache 
  unset -f error_cem error_cem_msg
}

error_cem_msg () {
  local FILE="$(basename ${BASH_SOURCE[0]})"
  local MSG="\033[0;31m\t\t (${FILE}) we cannot run "
  local MSG2="\033[0m"
  echo -e "${MSG}${1:?}${MSG2}" 2>"/dev/null"
}

error_cem () {
  error_cem_msg ${1:?}
  unset_all
  cd $(pwd -P) 2>"/dev/null"
  source stop_cocoa 2>"/dev/null"
  return 1
}

# ------------------------------------------------------------------------------
# ------------------------------ SET RUN MODES ---------------------------------
# ------------------------------------------------------------------------------
CACHE_FILE="${ROOTDIR:?}/.local/compile_local_packages.txt"
mode="default"        # default
mode_arg_seen=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --soft|--hard|--aggressive|--purge)
      if [[ $mode_arg_seen -eq 1 ]]; then
        error_cem "choose only one mode"
        return 1
      fi
      mode_arg_seen=1
      mode="${1#--}"
      ;;
    *)
      error_cem "Error: unknown option: $1"
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
  local LOCAL_DIR="${ROOTDIR:?}/.local"
  
  if [[ ! -d "${LOCAL_DIR:?}" ]]; then
    error_cem "(.local)  does not exist: ${LOCAL_DIR:?}"
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

case "$mode" in
  default) ;;
  soft)
    # FORCE RECOMPILE THEORY
    INIT=${#CORE[@]}
    END=$(( ${#CORE[@]} + ${#THEORY[@]} ))
    for (( i=$INIT; i<$END; i++ ));
    do
      CACHE[i]=0
    done
    ;;
  hard)
    # FORCE RECOMPILE THEORY + CORE
    INIT=0
    END=$(( ${#THEORY[@]} + ${#CORE[@]} ))
    for (( i=$INIT; i<$END; i++ ));
    do
      CACHE[i]=0
    done
    ;;
  aggressive)
    # FORCE RECOMPILE THEORY + CORE + ML
    INIT=0
    END=$(( ${#THEORY[@]} + ${#CORE[@]} + ${#ML[@]} ))
    for (( i=$INIT; i<$END; i++ ));
    do
      CACHE[i]=0
    done
    ;;
  purge)
    # FORCE RECOMPILE THEORY + CORE + ML + LIKELIHOOD
    INIT=0
    END=${#SCRIPTS[@]}
    for (( i=$INIT; i<$END; i++ ));
    do
      CACHE[i]=0
    done
    ;;
  *)
    echo "Error: invalid mode: $mode" >&2
    exit 1
    ;;
esac

# ----------------------------------------------------------------------------
# ------------------------ COMPILE EXTERNAL MODULES --------------------------
# ----------------------------------------------------------------------------

ptop2 'COMPILING EXTERNAL MODULES' || return 1

declare -i ERRORCODE=0

for (( i=0; i<${#SCRIPTS[@]}; i++ ));
do
  cdroot;

  [[ "${CACHE[i]}" -eq 1 ]] && continue

  ( export LD_LIBRARY_PATH=${CONDA_PREFIX:?}/lib:$LD_LIBRARY_PATH && \
    export LD_LIBRARY_PATH=${ROOTDIR:?}/.local/lib:$LD_LIBRARY_PATH && \
    if [ -n "${COCOA_OUTPUT_DEBUG}" ]; then  
      set -exo pipefail; 
    fi
    source "${ROOTDIR:?}/installation_scripts/${SCRIPTS[$i]}" 
  )
  
  rc=$?
  
  if [[ ${rc:?} -eq 55 ]]; then
    CACHE[i]=1
    save_cache
  elif [[ ${rc:?} -eq 99 ]]; then
    continue
  else
    if [ -n "${COCOA_OUTPUT_VERBOSE}" ]; then
      error_cem "script ${SCRIPTS[$i]}";
      return 1
    else
      error_cem_msg "script ${SCRIPTS[$i]}";
      ERRORCODE=1
    fi
  fi
done

# ----------------------------------------------------------------------------
# ------------------------ COMPILE EXTERNAL PROJECTS -------------------------
# ----------------------------------------------------------------------------

( source "${ROOTDIR:?}/installation_scripts/compile_all_projects.sh" )
if [ $? -ne 0 ]; then
  if [ -n "${COCOA_OUTPUT_VERBOSE}" ]; then
    error_cem "script compile_all_projects.sh"
    return 1
  else
    error_cem_msg "script compile_all_projects.sh"
    ERRORCODE=1
  fi
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if [ ${ERRORCODE:?} -ne 0 ]; then
  error_cem "compile_cocoa.sh";
  return 1
fi

unset_all
pbottom2 "COMPILING EXTERNAL MODULES" || return 1
source stop_cocoa.sh

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------