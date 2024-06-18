#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
if [ -n "${START_COCOA_DONE}" ]; then
  source stop_cocoa
fi

fail () {
    export FAILMSG="\033[0;31m WE CANNOT RUN "
    export FAILMSG2="\033[0m"
    echo -e "${FAILMSG} ${1} ${FAILMSG2}"
    unset FAILMSG
    unset FAILMSG2
    unset fail
    cd $(pwd -P)
    source stop_cocoa
}

# ----------------------------------------------------------------------------
# ------------------------------ Basic Settings ------------------------------
# ----------------------------------------------------------------------------

source $(pwd -P)/.save_old_flags.sh
if [ $? -ne 0 ]; then
  fail 'SCRIPT .SAVE_OLD_FLAGS.SH'
  source stop_cocoa
  return 1
fi

if [ -n "${SET_INSTALLATION_OPTIONS}" ]; then
  source $SET_INSTALLATION_OPTIONS
else
  source set_installation_options.sh
fi

# BASIC TEST IF A CONDA ENV IS ACTIVATED
if [ -n "${MINICONDA_INSTALLATION}" ]; then
  if [ -z ${CONDA_PREFIX} ]; then
    fail "start_cocoa.sh: CONDA ENVIRONMENT NOT ACTIVATED"
    return 1
  fi
fi

# ----------------------------------------------------------------------------
# ---------------------- Activate Virtual Environment ------------------------
# ----------------------------------------------------------------------------

source $ROOTDIR/.local/bin/activate
if [ $? -ne 0 ]; then
  echo "COCOA PRIVATE PYTHON ENV ACTIVATION"
  return 1
fi

source $(pwd -P)/.set_new_flags.sh
if [ -z "${MAKE_NUM_THREADS}" ]; then
  fail 'SCRIPT .SET_NEW_FLAGS.SH'
  return 1
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
  
  export CLIK_PATH=$ROOTDIR/.local

  export CLIK_DATA=$ROOTDIR/.local/share/clik

  export CLIK_PLUGIN=rel2015
fi

# ----------------------------------------------------------------------------
# ------------------------ START EXTERNAL PROJECTS ---------------------------
# ----------------------------------------------------------------------------
source $ROOTDIR/projects/start_all.sh

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
unset fail 
export START_COCOA_DONE=1
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------