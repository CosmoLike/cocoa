#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
export ROOTDIR=$(pwd -P)
if [ $? -ne 0 ]; then
  return 1;
fi

source "${ROOTDIR:?}/installation_scripts/flags_impl_unset_keys.sh" 
if [ $? -ne 0 ]; then
  return 1;
fi

# ------------------------------------------------------------------------------
# HOW COCOA SHOULD BE INSTALLED? -----------------------------------------------
# ------------------------------------------------------------------------------
export MINICONDA_INSTALLATION=1
#export MANUAL_INSTALLATION=1

# ------------------------------------------------------------------------------
# VERBOSE AS DEBUG TOOL --------------------------------------------------------
# ------------------------------------------------------------------------------
#export COCOA_OUTPUT_VERBOSE=1

# ------------------------------------------------------------------------------
# If set, COSMOLIKE will compile with DEBUG flags ------------------------------
# ------------------------------------------------------------------------------
#export COSMOLIKE_DEBUG_MODE=1

# ------------------------------------------------------------------------------
# The flags below allow users to skip downloading specific datasets ------------
# ------------------------------------------------------------------------------
# export IGNORE_ACTDR6_DATA=1
# export IGNORE_BAO_DATA=1
# export IGNORE_BICEP_CMB_DATA=1
# export IGNORE_HOLICOW_STRONG_LENSING_DATA=1
# export IGNORE_SN_DATA=1
# export IGNORE_SPT_CMB_DATA=1
export IGNORE_SIMONS_OBSERVATORY_CMB_DATA=1
# export IGNORE_PLANCK_CMB_DATA=1
export IGNORE_CAMSPEC_CMB_DATA=1
export IGNORE_LIPOP_CMB_DATA=1

# ------------------------------------------------------------------------------
# If set, compile_planck.sh uses click like code from github.com/benabed/clik
# ------------------------------------------------------------------------------
export USE_SPT_CLIK_PLANCK=1

# ------------------------------------------------------------------------------
# THREADING COMPILATION/INSTALLATION OF LIBRARIES ------------------------------
# ------------------------------------------------------------------------------
export MAKE_NUM_THREADS=4

# ------------------------------------------------------------------------------
# If not set, pip_core_packages.sh will install several ML packages ------------
# ------------------------------------------------------------------------------
export IGNORE_EMULATOR_CPU_PIP_PACKAGES=1
export IGNORE_EMULATOR_GPU_PIP_PACKAGES=1

# ------------------------------------------------------------------------------
# DERIVED & RARELY CHANGED FLAGS -----------------------------------------------
# ------------------------------------------------------------------------------

if [ -n "${MANUAL_INSTALLATION}" ]; then      
  source "${ROOTDIR:?}/installation_scripts/flags_manual_installation.sh" 
  if [ $? -ne 0 ]; then
    return 1;
  fi
elif [ -n "${MINICONDA_INSTALLATION}" ]; then
  source "${ROOTDIR:?}/installation_scripts/flags_miniconda_installation.sh"
  if [ $? -ne 0 ]; then
    return 1;
  fi
fi

# `flags_derived.sh` also contains many rarely used flags (useful to debug)
source "${ROOTDIR:?}/installation_scripts/flags_derived.sh"
if [ $? -ne 0 ]; then
  return 1;
fi

# ------------------------------------------------------------------------------
# PACKAGE URL AND VERSIONS. CHANGES IN THE COMMIT ID MAY BREAK COCOA -----------
# ------------------------------------------------------------------------------

export COBAYA_URL="https://github.com/CobayaSampler/cobaya.git"
export COBAYA_GIT_COMMIT="2636ea9ed399c35c5d276de1acb15aaafbcab10c"

export HILLIPOP_URL="https://github.com/planck-npipe/hillipop.git"
export HILLIPOP_GIT_COMMIT="cc9cbe31991d4662522241543a46d44d2cdec251"

export LOLLIPOP_URL="https://github.com/planck-npipe/lollipop.git"
export LOLLIPOP_GIT_COMMIT="280a9c93d33bc6a058d6bf769ec82d9f7fdbd2b6"

export LIPOP_DATA_URL="https://portal.nersc.gov/cfs/cmb/planck2020/likelihoods"
export LIPOP_DATA_VERSION=4.2

export SPT3G_DATA_URL='https://github.com/SouthPoleTelescope/spt3g_y1_dist.git'
export SPT3G_DATA_GIT_COMMIT="66da8e9e2f325024566fe13245788bf8ede897bc"
export SPT_3G_NAME="spt_3g"

export HOLICOW_DATA_URL='https://github.com/shsuyu/H0LiCOW-public.git'
export HOLICOW_DATA_GIT_COMMIT="f792647d1fd6c09d9e052fef526669cbd702ab82"

export POLY_URL="https://github.com/PolyChord/PolyChordLite.git"
export POLYCHORD_GIT_COMMIT="daba49d1385d065122db76a2b384050f9e95d278"
export POLY_NAME="PolyChordLite"

export CAMB_URL="https://github.com/cmbant/CAMB"
export CAMB_GIT_COMMIT="45d1c3d27e7480c0f9a82c98522c17ed422dd408"
export CAMB_NAME='CAMB'

export CLASS_URL="https://github.com/lesgourg/class_public.git"
export CLASS_GIT_COMMIT="8df566c1ff2d0b3e40e106567c435575aea337be"
export CLASS_NAME="class_public"

export ACTDR4_URL="https://github.com/ACTCollaboration/pyactlike"
export ACTDR4_GIT_COMMIT="1cac8c5d047bc2cad991890f2ebf1d8e3fb483b3"
export ACTDR4_NAME="pyactlike"

export ACT_DR6_DATA_URL="https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr6/likelihood/data"
export ACT_DR6_DATA_FILE="ACT_dr6_likelihood_v1.2.tgz"

export SO_DATA_URL="https://portal.nersc.gov/cfs/sobs/users/MFLike_data"
# Cocoa can download multiple versions of the data (to reproduce existing work)
# This is only possible because each version is saved in a separate folder
export SO_DATA_VERSION="v0.7.1 v0.8"

export VELOCILEPTORS_URL="https://github.com/sfschen/velocileptors.git"
export VELOCILEPTORS_GIT_COMMIT="889a0c98895831eb23b250a26162cfb8a93237bd"
export VELOCILEPTORS_NAME="velocileptors"

# ------------------------------------------------------------------------------
# The keys below control which packages will be installed and compiled 
# ------------------------------------------------------------------------------
#export IGNORE_CAMB_CODE=1
#export IGNORE_CLASS_CODE=1
#export IGNORE_COSMOLIKE_CODE=1
#export IGNORE_POLYCHORD_SAMPLER_CODE=1
#export IGNORE_PLANCK_LIKELIHOOD_CODE=1
#export IGNORE_ACTDR4_CODE=1
#export IGNORE_ACTDR6_CODE=1
#export IGNORE_CPP_CUBA_INSTALLATION=1
#export IGNORE_VELOCILEPTORS_CODE=1
#export IGNORE_SIMONS_OBSERVATORY_LIKELIHOOD_CODE=1
#export IGNORE_CAMSPEC_LIKELIHOOD_CODE=1
#export IGNORE_LIPOP_LIKELIHOOD_CODE=1
#export IGNORE_COBAYA_CODE=1
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------