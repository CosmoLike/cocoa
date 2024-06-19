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
fail_scipack () {
  export MSG="\033[0;31m (setup_cocoa_installation_package.sh) WE CANNOT RUN "
  export MSG2="\033[0m"
  echo -e "${MSG} ${1} ${MSG2}"
  unset MSG
  unset MSG2
  unset fail_scipack
  cd $(pwd -P)
  source stop_cocoa.sh
}

source .save_old_flags.sh
if [ $? -ne 0 ]; then
  fail_scipack 'SCRIPT .SAVE_OLD_FLAGS.SH'
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
    fail_scipack "CONDA ENVIRONMENT NOT ACTIVATED"
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
  fail_scipack 'SCRIPT .SET_NEW_FLAGS.SH'
  return 1
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

if [ -z "${MINICONDA_INSTALLATION}" ]; then
  source $ROOTDIR/installation_scripts/setup_xz.sh
  if [ $? -ne 0 ]; then
    fail_scipack 'SCRIPT SETUP_XZ.SH'
    return 1
  fi
fi

source ./installation_scripts/setup_installation_libraries.sh
if [ $? -ne 0 ]; then
  fail_scipack 'SCRIPT SETUP_INSTALLATION_LIBRARIES.SH'
  return 1
fi
  
if [ -z "${DEBUG_SKIP_FILE_DECOMPRESSION_SETUP_COCOA}" ]; then
  source $ROOTDIR/installation_scripts/setup_decompress_files.sh
  if [ $? -ne 0 ]; then
    fail_scipack 'SCRIPT SETUP_DECOMPRESS_FILES.SH'
    return 1
  fi
fi

if [ -z "${MINICONDA_INSTALLATION}" ]; then
  source ./installation_scripts/setup_cmake.sh
  if [ $? -ne 0 ]; then
    fail_scipack 'SCRIPT SETUP_CMAKE.SH'
    return 1
  fi
fi

if [ -z "${MINICONDA_INSTALLATION}" ]; then
  source $ROOTDIR/installation_scripts/setup_binutils.sh
  if [ $? -ne 0 ]; then
    fail_scipack 'SCRIPT SETUP_BUNUTILS.SH'
    return 1
  fi
fi

if [ -z "${MINICONDA_INSTALLATION}" ]; then
  source $ROOTDIR/installation_scripts/setup_hdf5.sh
  if [ $? -ne 0 ]; then
    fail_scipack 'SCRIPT SETUP_HDF5.SH'
    return 1
  fi
fi

if [ -z "${MINICONDA_INSTALLATION}" ]; then
  source $ROOTDIR/installation_scripts/setup_openblas.sh
  if [ $? -ne 0 ]; then
    fail_scipack 'SCRIPT SETUP_OPENBLAS.SH'
    return 1
  fi
fi

source $ROOTDIR/installation_scripts/setup_pip_packages.sh
if [ $? -ne 0 ]; then
  fail_scipack 'SCRIPT SETUP_PIP_PACKAGES.SH'
  return 1
fi

source $ROOTDIR/installation_scripts/setup_fortran_packages.sh
if [ $? -ne 0 ]; then
  fail_scipack 'SCRIPT SETUP_FORTRAN_PACKAGES.SH'
  return 1
fi

source $ROOTDIR/installation_scripts/setup_cpp_packages.sh
if [ $? -ne 0 ]; then
  fail_scipack 'SCRIPT SETUP_CPP_PACKAGES.SH'
  return 1
fi

source $ROOTDIR/installation_scripts/setup_c_packages.sh
if [ $? -ne 0 ]; then
  fail_scipack 'SCRIPT SETUP_C_PACKAGES.SH'
  return 1
fi

source $ROOTDIR/installation_scripts/setup_update_cobaya.sh
if [ $? -ne 0 ]; then
  fail_scipack 'SCRIPT SETUP_UPDATE_COBAYA.SH'
  return 1
fi

source $ROOTDIR/installation_scripts/setup_polychord.sh
if [ $? -ne 0 ]; then
  fail_scipack 'SCRIPT SETUP_POLYCHORD.SH'
  return 1
fi

source $ROOTDIR/installation_scripts/setup_camb.sh
if [ $? -ne 0 ]; then
  fail_scipack 'SCRIPT SETUP_CAMB.SH'
  return 1
fi

source $ROOTDIR/installation_scripts/setup_class.sh
if [ $? -ne 0 ]; then
  fail_scipack 'SCRIPT SETUP_CLASS.SH'
  return 1
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
unset fail_scipack
source stop_cocoa.sh
pbottom2 'SETUP COCOA INSTALLATION PACKAGES'
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------