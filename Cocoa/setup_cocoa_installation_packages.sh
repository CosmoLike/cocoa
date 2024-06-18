#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
echo -e '\033[1;45m''SETUP COCOA INSTALLATION PACKAGES''\033[0m'

if [ -n "${ROOTDIR}" ]; then
  source stop_cocoa.sh
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ------------------------------ Basic Settings ------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
fail () {
  export FAILMSG="\033[0;31m WE CANNOT RUN "
  export FAILMSG2="\033[0m"
  echo -e "${FAILMSG} ${1} ${FAILMSG2}"
  unset FAILMSG
  unset FAILMSG2
  unset fail
  cd $(pwd -P)
  source stop_cocoa.sh
}

source .save_old_flags.sh
if [ $? -ne 0 ]; then
  fail 'SCRIPT .SAVE_OLD_FLAGS.SH'
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
    fail "setup_cocoa_installation_package.sh: CONDA ENVIRONMENT NOT ACTIVATED"
    return 1
  fi
fi

# ----------------------------------------------------------------------------
# ---------------------- Activate Virtual Environment ------------------------
# ----------------------------------------------------------------------------
cd $ROOTDIR/../

if [ -n "${DONT_USE_SYSTEM_PIP_PACKAGES}" ]; then
  $GLOBALPYTHON3 -m venv $ROOTDIR/.local/
else
  $GLOBALPYTHON3 -m venv $ROOTDIR/.local/ --system-site-packages
fi

cd $ROOTDIR

source $ROOTDIR/.local/bin/activate

source $(pwd -P)/.set_new_flags.sh
if [ $? -ne 0 ]; then
  fail 'SCRIPT .SET_NEW_FLAGS.SH'
  return 1
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${MINICONDA_INSTALLATION}" ]; then
  source $ROOTDIR/installation_scripts/setup_xz.sh
  if [ $? -ne 0 ]; then
    cd $ROOTDIR
    source stop_cocoa.sh
    return 1
  fi
fi

source ./installation_scripts/setup_installation_libraries.sh
if [ $? -ne 0 ]; then
  cd $ROOTDIR
  source stop_cocoa.sh
  return 1
fi
  
if [ -z "${DEBUG_SKIP_FILE_DECOMPRESSION_SETUP_COCOA}" ]; then
  source $ROOTDIR/installation_scripts/setup_decompress_files.sh
  if [ $? -ne 0 ]; then
    cd $ROOTDIR
    source stop_cocoa.sh
    return 1
  fi
fi

if [ -z "${MINICONDA_INSTALLATION}" ]; then
  source ./installation_scripts/setup_cmake.sh
  if [ $? -ne 0 ]; then
    cd $ROOTDIR
    source stop_cocoa.sh
    return 1
  fi
fi

if [ -z "${MINICONDA_INSTALLATION}" ]; then
  source $ROOTDIR/installation_scripts/setup_binutils.sh
  if [ $? -ne 0 ]; then
    cd $ROOTDIR
    source stop_cocoa.sh
    return 1
  fi
fi

if [ -z "${MINICONDA_INSTALLATION}" ]; then
  source $ROOTDIR/installation_scripts/setup_hdf5.sh
  if [ $? -ne 0 ]; then
    cd $ROOTDIR
    source stop_cocoa.sh
    return 1
  fi
fi

if [ -z "${MINICONDA_INSTALLATION}" ]; then
  source $ROOTDIR/installation_scripts/setup_openblas.sh
  if [ $? -ne 0 ]; then
    cd $ROOTDIR
    source stop_cocoa.sh
    return 1
  fi
fi

source $ROOTDIR/installation_scripts/setup_pip_packages.sh
if [ $? -ne 0 ]; then
  cd $ROOTDIR
  source stop_cocoa.sh
  return 1
fi

source $ROOTDIR/installation_scripts/setup_fortran_packages.sh
if [ $? -ne 0 ]; then
  cd $ROOTDIR
  source stop_cocoa.sh
  return 1
fi

source $ROOTDIR/installation_scripts/setup_cpp_packages.sh
if [ $? -ne 0 ]; then
  cd $ROOTDIR
  source stop_cocoa.sh
  return 1
fi

source $ROOTDIR/installation_scripts/setup_c_packages.sh
if [ $? -ne 0 ]; then
  cd $ROOTDIR
  source stop_cocoa.sh
  return 1
fi

source $ROOTDIR/installation_scripts/setup_update_cobaya.sh
if [ $? -ne 0 ]; then
  cd $ROOTDIR
  source stop_cocoa.sh
  return 1
fi

source $ROOTDIR/installation_scripts/setup_polychord.sh
if [ $? -ne 0 ]; then
  cd $ROOTDIR
  source stop_cocoa.sh
  return 1
fi

source $ROOTDIR/installation_scripts/setup_camb.sh
if [ $? -ne 0 ]; then
  cd $ROOTDIR
  source stop_cocoa.sh
  return 1
fi

source $ROOTDIR/installation_scripts/setup_class.sh
if [ $? -ne 0 ]; then
  cd $ROOTDIR
  source stop_cocoa.sh
  return 1
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
unset fail
source stop_cocoa
pbottom2 'SETUP COCOA INSTALLATION PACKAGES'
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------