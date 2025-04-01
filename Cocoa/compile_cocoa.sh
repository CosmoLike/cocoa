#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
error_cem_msg () {
  local FILE="$(basename ${BASH_SOURCE[0]})"
  local MSG="\033[0;31m\t\t (${FILE}) we cannot run "
  local MSG2="\033[0m"
  echo -e "${MSG}${1:?}${MSG2}" 2>"/dev/null"
}

error_cem () {
  error_cem_msg ${1:?}
  unset -v SCRIPTS ERRORCODE
  unset -f error_cem error_cem_msg
  cd $(pwd -P) 2>"/dev/null"
  source stop_cocoa 2>"/dev/null"
  return 1
}

if [ -n "${ROOTDIR}" ]; then
  source stop_cocoa.sh
fi
source start_cocoa.sh || { error_cem "start_cocoa.sh"; return 1; }

# ----------------------------------------------------------------------------
# ------------------------ COMPILE EXTERNAL MODULES --------------------------
# ----------------------------------------------------------------------------

ptop2 'COMPILING EXTERNAL MODULES' || return 1

declare -i ERRORCODE=0

declare -a SCRIPTS=("compile_core_packages.sh"
                     "compile_hyrec2.sh"
                     "compile_cosmorec.sh"
                     "compile_camb.sh"
                     "compile_mgcamb.sh"
                     "compile_class.sh" 
                     "compile_planck.sh" 
                     "compile_polychord.sh"
                     "compile_velocileptors.sh"
                     "compile_ee2.sh"
                     "compile_fgspectra.sh"
                     "compile_lipop.sh"
                     "compile_simons_observatory.sh"
                     "compile_act_dr4.sh"
                     "compile_act_dr6.sh"
                     ) # T = TMP

for (( i=0; i<${#SCRIPTS[@]}; i++ ));
do

  cdroot; 

  ( source "${ROOTDIR:?}/installation_scripts/${SCRIPTS[$i]}" )
  if [ $? -ne 0 ]; then
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
  error_cem "setup_cocoa.py";
  return 1
fi

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

unset -v ERRORCODE SCRIPTS
unset -f error_cem error_cem_msg

pbottom2 "COMPILING EXTERNAL MODULES" || return 1

source stop_cocoa.sh

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------