#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
echo -e '\033[1;45m''COMPILING EXTERNAL MODULES''\033[0m'

if [ -n "${ROOTDIR}" ]; then
  source stop_cocoa
fi

source start_cocoa

# ----------------------------------------------------------------------------
# ------------------------------ COMPILE COBAYA ------------------------------
# ----------------------------------------------------------------------------
echo -e '\033[1;34m''\tINSTALLING COBAYA (VIA PIP)''\033[0m'

if [ -z "${DEBUG_PIP_OUTPUT}" ]; then
  export OUTPUT_PIP_1="/dev/null"
  export OUTPUT_PIP_2="/dev/null"
else
  export OUTPUT_PIP_1="/dev/tty"
  export OUTPUT_PIP_2="/dev/tty"
fi

$PIP3 install --editable cobaya --prefix=$ROOTDIR/.local > ${OUTPUT_PIP_1} 2> ${OUTPUT_PIP_2}
if [ $? -ne 0 ]; then
  echo -e '\033[0;31m'"PIP COULD NOT RUN \e[3mPIP INSTALL COBAYA"'\033[0m'
  cd $ROOTDIR
  return 1
else
  echo -e '\033[0;32m'"PIP RUN \e[3mPIP INSTALL COBAYA\e[0m\e\033[0;32m DONE"'\033[0m'
fi

echo -e '\033[1;34m''\t\e[4mINSTALLING COBAYA (VIA PIP) DONE''\033[0m'

# ----------------------------------------------------------------------------
# ------------------------ COMPILE EXTERNAL MODULES --------------------------
# ----------------------------------------------------------------------------

source $ROOTDIR/installation_scripts/compile_camb.sh

source $ROOTDIR/installation_scripts/compile_class.sh

source $ROOTDIR/installation_scripts/compile_polychord.sh

source $ROOTDIR/installation_scripts/compile_planck.sh

source $ROOTDIR/installation_scripts/compile_act.sh

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
source stop_cocoa
echo -e '\033[1;45m''\e[4mCOMPILING EXTERNAL MODULES DONE''\033[0m'
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------