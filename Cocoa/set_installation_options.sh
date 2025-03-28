#!/bin/bash

# ------------------------------------------------------------------------------
# CLEANING PREVIOUSLY SET ENVIRONMENT (DO NOT CHANGE) --------------------------
# ------------------------------------------------------------------------------
source "$(pwd -P)/installation_scripts/flags_impl_unset_keys.sh" 

# ------------------------------------------------------------------------------
# DEFINING ROOTDIR (DO NOT CHANGE) ---------------------------------------------
# ------------------------------------------------------------------------------
export ROOTDIR=$(pwd -P)
if [ $? -ne 0 ]; then
  return 1;
fi

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
#export IGNORE_ACTDR6_DATA=1
# export IGNORE_BAO_DATA=1
export IGNORE_BICEP_CMB_DATA=1
# export IGNORE_HOLICOW_STRONG_LENSING_DATA=1
# export IGNORE_SN_DATA=1
# export IGNORE_SPT_CMB_DATA=1
export IGNORE_SIMONS_OBSERVATORY_CMB_DATA=1
# export IGNORE_PLANCK_CMB_DATA=1
export IGNORE_CAMSPEC_CMB_DATA=1
export IGNORE_LIPOP_CMB_DATA=1

# ------------------------------------------------------------------------------
# We download packages from servers that slowdown transfer after many runs. So, 
# if rerun setup_cocoa to download some non-core lib, these flags may be handy
# ------------------------------------------------------------------------------
#export IGNORE_CORE_INSTALLATION=1
#export IGNORE_PIP_CORE_INSTALLATION=1

# ------------------------------------------------------------------------------
# The keys below control which packages will be installed and compiled 
# ------------------------------------------------------------------------------
#export IGNORE_COBAYA_CODE=1

#export IGNORE_CAMB_CODE=1
#export IGNORE_CLASS_CODE=1
#export IGNORE_COSMOLIKE_CODE=1
#export IGNORE_POLYCHORD_SAMPLER_CODE=1
#export IGNORE_PLANCK_LIKELIHOOD_CODE=1
#export IGNORE_ACTDR4_CODE=1
#export IGNORE_ACTDR6_CODE=1
export IGNORE_CPP_CUBA_INSTALLATION=1
export IGNORE_VELOCILEPTORS_CODE=1
export IGNORE_SIMONS_OBSERVATORY_LIKELIHOOD_CODE=1
export IGNORE_CAMSPEC_LIKELIHOOD_CODE=1
export IGNORE_LIPOP_LIKELIHOOD_CODE=1
export IGNORE_HYREC_CODE=1
export IGNORE_COSMOREC_CODE=1
export IGNORE_MGCAMB_CODE=1

#Many cosmolike projects (including LSST-Y1) require euclid emulator
#export IGNORE_EUCLID_EMULATOR_V2_CODE=1
#export IGNORE_COSMOLIKE_LSSTY1_CODE=1

# ------------------------------------------------------------------------------
# If OVERWRITE_EXISTING_XXX_CODE is set, setup_XXX overwrites existing PACKAGES
# overwrite = delete existing PACKAGE folder and install it again --------------
# ------------------------------------------------------------------------------
#export OVERWRITE_EXISTING_COSMOLIKE_CODE=1
export OVERWRITE_EXISTING_CAMB_CODE=1
export OVERWRITE_EXISTING_MGCAMB_CODE=1
export OVERWRITE_EXISTING_CLASS_CODE=1
export OVERWRITE_EXISTING_HYREC_CODE=1
export OVERWRITE_EXISTING_COSMOREC_CODE=1

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
# Adopted Python version -------------------------------------------------------
# ------------------------------------------------------------------------------
export PYTHON_VERSION=3.9

# ------------------------------------------------------------------------------
# HOW COCOA CORE PACKAGES SHOULD BE INSTALLED? ---------------------------------
# ------------------------------------------------------------------------------
export MINICONDA_INSTALLATION=1
#export MANUAL_INSTALLATION=1

# ------------------------------------------------------------------------------
# DERIVED & RARELY CHANGED FLAGS (DO NOT CHANGE) -------------------------------
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

# We decided to install C++ Armadillo library locally 
# to link it against lapack & OpenBlas & arpack
unset IGNORE_CPP_ARMA_INSTALLATION

# ------------------------------------------------------------------------------
# PACKAGE URL AND VERSIONS. CHANGES IN THE COMMIT ID MAY BREAK COCOA -----------
# ------------------------------------------------------------------------------

# This flag saves a lot of time when running setup_cocoa.py 
# Why? Some git repos can be hundreds of MegaBytes (Class is 500 MegaBytes) 
# But, this can create problem if GIT_COMMIT < LAST COMMIT - GIT_MAXIMUM_DEPTH
export GIT_CLONE_MAXIMUM_DEPTH=40

export COBAYA_URL="https://github.com/CobayaSampler/cobaya.git"
#export COBAYA_GIT_COMMIT="2636ea9ed399c35c5d276de1acb15aaafbcab10c"
export COBAYA_GIT_COMMIT="86943d81d48d2edb2961b17077461df9e799f4d1"

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
export SO_DATA_VERSION="v0.8"

export VELOCILEPTORS_URL="https://github.com/sfschen/velocileptors.git"
export VELOCILEPTORS_GIT_COMMIT="889a0c98895831eb23b250a26162cfb8a93237bd"
export VELOCILEPTORS_NAME="velocileptors"

export FGSPECTRA_URL="https://github.com/simonsobs/fgspectra.git"
export FGSPECTRA_GIT_COMMIT="87206ac7f919bcf82b08dee180f0821d038d24d2"
export FGSPECTRA_NAME="fgspectra"

export SO_MFLIKE_URL="https://github.com/simonsobs/LAT_MFLike.git"
export SO_MFLIKE_GIT_COMMIT="660e0fa06e1c68335b7def416c256dfd316564d4"

export SO_SYSLIB_URL="https://github.com/simonsobs/syslibrary.git"
export SO_SYSLIB_GIT_COMMIT="c4a112dc2324647395b6c91d4747e3ab19754413"
export SO_SYSLIB_NAME="SOSYSLIB"

export EE2_URL="https://github.com/miknab/EuclidEmulator2.git"
export EE2_GIT_COMMIT="ff59f6683069417f6b4d2fb5d59197044d424445"
export EE2_NAME="euclidemu2"

export HYREC_URL="https://github.com/nanoomlee/HYREC-2.git"
export HYREC_GIT_COMMIT="09e8243d0e08edd3603a94dfbc445ae06cafe139"
export HYREC_NAME="hyrec2"

export COSMOREC_URL="https://www.cita.utoronto.ca/~jchluba/Recombination/_Downloads_/"
export COSMOREC_CODE_FILE="CosmoRec.v2.0.3b"
export COSMOREC_CODE_FILE_EXT="tar.gz"
export COSMOREC_NAME="cosmorec"

export MGCAMB_URL="https://github.com/sfu-cosmo/MGCobaya.git"
export MGCAMB_GIT_COMMIT="443c4a733db687ac18e918b8ed09b45003a8c4ca"
export MGCAMB_NAME='MGCAMB'

# --------------------------------------------------------------------
# --------------------------------------------------------------------
# Cosmolike projects below
# --------------------------------------------------------------------
# --------------------------------------------------------------------
export LSST_Y1_URL="https://github.com/CosmoLike/cocoa_lsst_y1.git"
export LSST_Y1_NAME="lsst_y1"
#export LSST_Y1_BRANCH="main"
#export LSST_Y1_COMMIT="8f7f0c40784fcae718006b6c6715c58fb0b085b0"
export LSST_Y1_TAG="v4.0-beta17"

export DES_Y3_URL="https://github.com/CosmoLike/cocoa_des_y3.git"
export DES_Y3_NAME="des_y3"
#export DES_Y3_BRANCH="main"
#export DES_Y3_COMMIT="2d7e05944bf9c4724dcb5407e7cb76c992993807"
export DES_Y3_TAG="v4.0-beta17"

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------