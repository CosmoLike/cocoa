#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if command -v deactivate &> "/dev/null"; then
  deactivate 2>"/dev/null"
fi

error_stop_cocoa () {
  local FILE="$(basename ${BASH_SOURCE[0]})"
  local MSG="\033[0;31m (${FILE}) we cannot run "
  local MSG2="\033[0m"
  echo -e "${MSG}${1:?}${MSG2}" 2>"/dev/null"
  unset -f error
  cd $(pwd -P) 2>"/dev/null"
  source stop_cocoa 2>"/dev/null"
  return 1
}

source $(pwd -P)/.recover_old_flags.sh
if [ $? -ne 0 ]; then
  error_stop_cocoa 'script .recover_old_flags.sh'; return 1
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
  if [ $? -ne 0 ]; then
    error_stop_cocoa 'script projects/stop_all.sh'; return 1
  fi
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

source ./installation_scripts/.impl_unset_keys.sh
if [ $? -ne 0 ]; then
  error_stop_cocoa 'script .impl_unset_keys.sh'; return 1
fi

unset -v ROOTDIR SETUP_COBAYA START_COCOA_DONE fail

unset -v SETUP_PREREQUISITE_DONE SET_INSTALLATION_OPTIONS

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------