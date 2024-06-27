#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
error_cem () {
  local FILE="$(basename ${BASH_SOURCE[0]})"
  local MSG="\033[0;31m\t\t (${FILE}) we cannot run "
  local MSG2="\033[0m"
  echo -e "${MSG}${1:?}${MSG2}" 2>"/dev/null"
  unset -v TSCRIPTS
  unset -f error_cem
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

ptop2 'COMPILING ALL EXTERNAL MODULES' || return 1

declare -a TSCRIPTS=("compile_camb.sh"
                     "compile_class.sh" 
                     "compile_planck.sh" 
                     "compile_act.sh"
                     "compile_polychord.sh"
                     ) # T = TMP

for (( i=0; i<${#TSCRIPTS[@]}; i++ ));
do

  cdroot; 

  ( source "${ROOTDIR:?}/installation_scripts/${TSCRIPTS[$i]}" )
  if [ $? -ne 0 ]; then
    error_cem "script ${TSCRIPTS[$i]}"; return 1
    return 1
  fi

done

unset -f error_cem

pbottom2 "COMPILING ALL EXTERNAL MODULES" || return 1

source stop_cocoa.sh

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------