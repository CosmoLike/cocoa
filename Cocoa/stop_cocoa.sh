#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if ! command -v deactivate &> /dev/null
then
  echo "WARNING: COCOA PRIVATE PYTHON ENV DEACTIVATION FAILED"
else
  deactivate
fi

fail () {
    export FAILMSG="\033[0;31m WE CANNOT RUN"
    export FAILMSG2="\033[0m"
    echo -e "${FAILMSG} ${1} ${FAILMSG2}"
    unset FAILMSG
    unset FAILMSG2
    unset fail
}

source $(pwd -P)/.recover_old_flags.sh
if [ -z "${MAKE_NUM_THREADS}" ]; then
  fail 'SCRIPT .RECOVER_OLD_FLAGS.SH'
  source stop_cocoa.sh
  return 1
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
# ------------------------ STOP EXTERNAL PROJECTS ---------------------------
# ----------------------------------------------------------------------------

if [ -n "${ROOTDIR}" ]; then
  source $ROOTDIR/projects/stop_all.sh
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
source ./installation_scripts/.impl_unset_keys.sh
unset -v ROOTDIR SETUP_COBAYA START_COCOA_DONE fail
unset -v SETUP_PREREQUISITE_DONE SET_INSTALLATION_OPTIONS
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------