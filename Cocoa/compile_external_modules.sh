#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -n "${ROOTDIR}" ]; then
  source stop_cocoa.sh
fi
source start_cocoa.sh

ptop2 'COMPILING EXTERNAL MODULES' || return 1

# ----------------------------------------------------------------------------
# ------------------------------ COMPILE COBAYA ------------------------------
# ----------------------------------------------------------------------------

ptop "INSTALLING COBAYA (VIA PIP)" || return 1

if [ -z "${COCOA_OUTPUT_VERBOSE}" ]; then
  export OUT1="/dev/null"; export OUT2="/dev/null"
else
  export OUT1="/dev/tty"; export OUT2="/dev/tty"
fi

${PIP3:?} install --editable cobaya --prefix="${ROOTDIR:?}/.local" >"${OUT1:?}" 2>"${OUT2:?}"
if [ $? -ne 0 ]; then
  echo -e '\033[0;31m'"PIP COULD NOT RUN \e[3mPIP INSTALL COBAYA"'\033[0m'
  cd $ROOTDIR
  return 1
fi

pbottom "INSTALLING COBAYA (VIA PIP)" || return 1

# ----------------------------------------------------------------------------
# ------------------------ COMPILE EXTERNAL MODULES --------------------------
# ----------------------------------------------------------------------------

source "${ROOTDIR:?}"/installation_scripts/compile_camb.sh
if [ $? -ne 0 ]; then
  cd $ROOTDIR
  source stop_cocoa
  return 1
fi

source "${ROOTDIR:?}"/installation_scripts/compile_class.sh
if [ $? -ne 0 ]; then
  cd $ROOTDIR
  source stop_cocoa
  return 1
fi

source "${ROOTDIR:?}"/installation_scripts/compile_planck.sh
if [ $? -ne 0 ]; then
  cd $ROOTDIR
  source stop_cocoa
  return 1
fi

source "${ROOTDIR:?}"/installation_scripts/compile_act.sh
if [ $? -ne 0 ]; then
  cd $ROOTDIR
  source stop_cocoa
  return 1
fi

source "${ROOTDIR:?}"/installation_scripts/compile_polychord.sh
if [ $? -ne 0 ]; then
  cd $ROOTDIR
  source stop_cocoa
  return 1
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
pbottom2 "COMPILING EXTERNAL MODULES" || return 1
source stop_cocoa.sh
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------