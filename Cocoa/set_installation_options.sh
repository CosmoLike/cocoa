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
export COCOA_OUTPUT_VERBOSE=1

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
#export IGNORE_SN_DATA=1
export IGNORE_SPT_CMB_DATA=1
#export IGNORE_SIMONS_OBSERVATORY_CMB_DATA=1
#export IGNORE_PLANCK_CMB_DATA=1
export IGNORE_CAMSPEC_CMB_DATA=1
#export IGNORE_LIPOP_CMB_DATA=1
export IGNORE_COSMOPOWER_DATA=1
#export IGNORE_EMULTRF_DATA=1 #SaraivanovZhongZhu (SZZ) transformer-based emul

# ------------------------------------------------------------------------------
# The keys below control which packages will be installed and compiled 
# ------------------------------------------------------------------------------
#export IGNORE_COBAYA_CODE=1
#export IGNORE_CAMB_CODE=1
export IGNORE_CLASS_CODE=1 # Default: we just use CAMB (reduces compilation time)
#export IGNORE_COSMOLIKE_CODE=1
#export IGNORE_POLYCHORD_SAMPLER_CODE=1
#export IGNORE_PLANCK_LIKELIHOOD_CODE=1
export IGNORE_ACTDR4_CODE=1
#export IGNORE_ACTDR6_CODE=1
export IGNORE_CPP_CUBA_INSTALLATION=1
export IGNORE_FGSPECTRA_CODE=1
export IGNORE_VELOCILEPTORS_CODE=1
#export IGNORE_SIMONS_OBSERVATORY_LIKELIHOOD_CODE=1
export IGNORE_CAMSPEC_LIKELIHOOD_CODE=1
#export IGNORE_LIPOP_LIKELIHOOD_CODE=1
export IGNORE_HYREC_CODE=1
export IGNORE_COSMOREC_CODE=1
export IGNORE_MGCAMB_CODE=1
#export IGNORE_EMULTRF_CODE=1 #SaraivanovZhongZhu (SZZ) transformer-based emul
export IGNORE_COSMOPOWER_CODE=1
#export IGNORE_EUCLID_EMULATOR_V2_CODE=1
export IGNORE_DARK_EMULATOR_CODE=1

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# The keys below control which cosmolike projects will be installed and compiled 
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#export IGNORE_COSMOLIKE_LSST_Y1_CODE=1
export IGNORE_COSMOLIKE_DES_Y3_CODE=1
export IGNORE_COSMOLIKE_ROMAN_FOURIER_CODE=1
export IGNORE_COSMOLIKE_ROMAN_REAL_CODE=1

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# The keys below control which private projects (not public repo)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
export INSTALL_PRIVATE_AXIONS_PROJECT=1

# ------------------------------------------------------------------------------
# OVERWRITE_EXISTING_XXX_CODE=1 -> setup_cocoa overwrites existing PACKAGES ----
# overwrite: delete the existing PACKAGE folder and install it again -----------
# redownload: delete the compressed file and download data again ---------------
# These keys are only relevant if you run setup_cocoa multiple times -----------
# ------------------------------------------------------------------------------
export OVERWRITE_EXISTING_ALL_PACKAGES=1    # except cosmolike projects
#export OVERWRITE_EXISTING_COSMOLIKE_CODE=1 # dangerous (possible lost of uncommit work)
                                            # if unset, users must manually delete
                                            # project if wants setup_cocoa to reclone it
#export REDOWNLOAD_EXISTING_ALL_DATA=1      # warning: some data are many GB

#export OVERWRITE_EXISTING_PRIVATE_CODE=1   # dangerous (possible lost of uncommit work)
                                            # if unset, users must manually delete
                                            # project if wants setup_cocoa to reclone it
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
#export IGNORE_EMULATOR_GPU_PIP_PACKAGES=1

# ------------------------------------------------------------------------------
# Adopted Python version -------------------------------------------------------
# ------------------------------------------------------------------------------
export PYTHON_VERSION=3.10

# ------------------------------------------------------------------------------
# HOW COCOA CORE PACKAGES SHOULD BE INSTALLED? ---------------------------------
# ------------------------------------------------------------------------------
export MINICONDA_INSTALLATION=1
#export MANUAL_INSTALLATION=1

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# DERIVED & RARELY CHANGED FLAGS (DO NOT CHANGE) -------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
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
# -------------------- COMPATIBILITY/DEPENDENCIES ------------------------------
# ------------------------------------------------------------------------------
if [ -z "${IGNORE_ACTDR6_CODE}" ]; then
  unset -v IGNORE_FGSPECTRA_CODE
  unset -v IGNORE_SIMONS_OBSERVATORY_LIKELIHOOD_CODE
  unset -v IGNORE_COSMOREC_CODE
fi

if [ -z "${IGNORE_SIMONS_OBSERVATORY_LIKELIHOOD_CODE}" ]; then
  unset -v IGNORE_FGSPECTRA_CODE
fi

if [[ -z "${IGNORE_COSMOPOWER_CODE}" || -z "${IGNORE_EMULTRF_CODE}" ]]; then
  unset -v IGNORE_EMULATOR_GPU_PIP_PACKAGES
fi

#if [ -z "${IGNORE_EUCLID_EMULATOR_V2_CODE}" ]; then
#  unset -v IGNORE_CLASS_CODE
#fi

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# If OVERWRITE_EXISTING_XXX_CODE=1, the setup_cocoa overwrites existing PACKAGES
# overwrite means: delete existing PACKAGE folder and install it again ---------
# these keys are only relevant if you run setup_cocoa multiple times -----------
# ------------------------------------------------------------------------------
if [ -n "${OVERWRITE_EXISTING_ALL_PACKAGES}" ]; then
  export OVERWRITE_EXISTING_COCOA_PRIVATE_PYTHON_ENV=1
  export OVERWRITE_EXISTING_CORE_PACKAGES=1
  export OVERWRITE_EXISTING_COBAYA_CODE=1
  export OVERWRITE_EXISTING_CAMB_CODE=1
  export OVERWRITE_EXISTING_MGCAMB_CODE=1
  export OVERWRITE_EXISTING_CLASS_CODE=1
  export OVERWRITE_EXISTING_HYREC_CODE=1
  export OVERWRITE_EXISTING_COSMOREC_CODE=1
  export OVERWRITE_EXISTING_POLYCHORD_CODE=1
  export OVERWRITE_EXISTING_VELOCILEPTORS_CODE=1
  export OVERWRITE_EXISTING_EE2_CODE=1
  export OVERWRITE_EXISTING_BAO_DATA=1
  export OVERWRITE_EXISTING_SIMONS_OBSERVATORY_CODE=1
  export OVERWRITE_EXISTING_FGSPECTRA_DATA=1
  export OVERWRITE_EXISTING_LIPOP_CMB_CODE=1
  export OVERWRITE_EXISTING_LIPOP_CMB_DATA=1
  export OVERWRITE_EXISTING_ACTDR4_CMB_CODE=1
  export OVERWRITE_EXISTING_ACTDR4_CMB_DATA=1
  export OVERWRITE_EXISTING_ACTDR6_CMB_CODE=1
  export OVERWRITE_EXISTING_ACTDR6_CMB_DATA=1
  export OVERWRITE_EXISTING_BICEP_CMB_DATA=1
  export OVERWRITE_EXISTING_CAMPSPEC_CMB_DATA=1
  export OVERWRITE_EXISTING_SPT3G_CMB_DATA=1
  export OVERWRITE_EXISTING_PLANCK_CMB_DATA=1
  export OVERWRITE_EXISTING_SIMONS_OBSERVATORY_CMB_DATA=1
  export OVERWRITE_EXISTING_SN_DATA=1
  export OVERWRITE_EXISTING_HOLICOW_DATA=1
  export OVERWRITE_EXISTING_COSMOPOWER_CODE=1
  export OVERWRITE_EXISTING_EMULTRF_CODE=1
  export OVERWRITE_EXISTING_DARK_EMULATOR_CODE=1
fi

if [ -n "${REDOWNLOAD_EXISTING_ALL_DATA}" ]; then
  export REDOWNLOAD_EXISTING_CORE_PACKAGES=1
  export REDOWNLOAD_EXISTING_ACTDR6_CMB_DATA=1
  export REDOWNLOAD_EXISTING_LIPOP_CMB_DATA=1
  export REDOWNLOAD_EXISTING_SIMONS_OBSERVATORY_CMB_DATA=1
  export REDOWNLOAD_EXISTING_CAMPSPEC_CMB_DATA=1
  export OVERWRITE_EXISTING_EMULTRF_DATA=1
  export OVERWRITE_EXISTING_COSMOPOWER_DATA=1
fi
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# PACKAGE URL AND VERSIONS. CHANGES IN THE COMMIT ID MAY BREAK COCOA -----------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# This flag saves a lot of time when running setup_cocoa.py 
# Why? Some git repos can be hundreds of MegaBytes (Class is 500 MegaBytes) 
# But, this can create problem if GIT_COMMIT < LAST COMMIT - GIT_MAXIMUM_DEPTH
export SPDLOG_VERSION=v1.15.3
export GIT_CLONE_MAXIMUM_DEPTH=40

export COBAYA_URL="https://github.com/CobayaSampler/cobaya.git"
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
export CAMB_GIT_COMMIT="886b17cbc23137737b7ef4318d165aa7bdbbbbed"
export CAMB_NAME='CAMB'

export CLASS_URL="https://github.com/lesgourg/class_public.git"
export CLASS_GIT_COMMIT="0ceb7a9a4c1e444ef5d5d56a8328a0640be91b18"
export CLASS_NAME="class_public"

export ACTDR4_URL="https://github.com/ACTCollaboration/pyactlike"
export ACTDR4_GIT_COMMIT="1cac8c5d047bc2cad991890f2ebf1d8e3fb483b3"
export ACTDR4_NAME="pyactlike"

export ACTDR6_CMBONLY_URL="https://github.com/ACTCollaboration/DR6-ACT-lite.git"
export ACTDR6_CMBONLY_GIT_COMMIT="627aeafb88ae5ad1aa66b406bea2d65cfa66a27d"

export ACTDR6_MFLIKE_URL="https://github.com/ACTCollaboration/act_dr6_mflike.git"
export ACTDR6_MFLIKE_GIT_COMMIT=4249ff9e2d92f01d6a38d4adbe78ad34f83a33f7

export ACTDR6_LENSING_DATA_URL="https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr6/likelihood/data"
export ACTDR6_LENSING_DATA_FILE="ACT_dr6_likelihood_v1.2.tgz"

export ACTDR6_CMBONLY_DATA_URL="https://lambda.gsfc.nasa.gov/data/act/pspipe/sacc_files/"
export ACTDR6_CMBONLY_DATA_FILE="dr6_data_cmbonly.tar.gz"

export ACTDR6_MFLIKE_DATA_URL=https://lambda.gsfc.nasa.gov/data/act/pspipe/sacc_files/
export ACTDR6_MFLIKE_DATA_FILE=dr6_data.tar.gz

export SO_DATA_URL="https://portal.nersc.gov/cfs/sobs/users/MFLike_data"
export SO_DATA_VERSION="v0.8"

export VELOCILEPTORS_URL="https://github.com/sfschen/velocileptors.git"
export VELOCILEPTORS_GIT_COMMIT="889a0c98895831eb23b250a26162cfb8a93237bd"
export VELOCILEPTORS_NAME="velocileptors"

export FGSPECTRA_URL="https://github.com/simonsobs/fgspectra.git"
export FGSPECTRA_GIT_COMMIT="cd78a3a72274dac2b9e05bef1943370807c46146"
export FGSPECTRA_NAME="fgspectra"

export SO_MFLIKE_URL="https://github.com/simonsobs/LAT_MFLike.git"
export SO_MFLIKE_GIT_COMMIT="531267e8f046f72ddc5fe4ff88432871ad8c9cfd"

export SO_SYSLIB_URL="https://github.com/simonsobs/syslibrary.git"
export SO_SYSLIB_GIT_COMMIT="2471df981053a7526a441d2547eb8dde10d92f70"

export EE2_URL="https://github.com/vivianmiranda/EuclidEmulator2.git"
export EE2_GIT_COMMIT="df250250e31fcd15841b622b42142cb12f9fd68b"
# THIS GIT IS FOR THE ORIGNIAL EE2 BEFORE VM MODS
#export EE2_GIT_COMMIT="ff59f6683069417f6b4d2fb5d59197044d424445"
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

export PLANCK2018_SROLL2_URL="https://web.fe.infn.it/~pagano/low_ell_datasets/sroll2/"
export PLANCK2018_SROLL2_FILE="simall_100x143_sroll2_v3_EE_Aplanck.tgz"

export COSMOPOWER_SOLIKET_URL="https://github.com/simonsobs/SOLikeT.git"
export COSMOPOWER_SOLIKET_GIT_COMMIT="1d8333ea0007c88e7c2de192de39301884093cd8"

export COSMOPOWER_URL="https://github.com/alessiospuriomancini/cosmopower.git"
export COSMOPOWER_GIT_COMMIT="7cac5e71c975c06257b2f95f0dcea5dd09b0f45f"

export COSMOPOWER_URL_DATA="https://github.com/cosmopower-organization/jense_2024_emulators.git"
export COSMOPOWER_URL_DATA_COMMIT="4317635eed70289ee1ec6b3df828027173071e36"

export EMULTRF_URL="git@github.com:SBU-COSMOLIKE/emulators_code.git"
#export EMULTRF_GIT_COMMIT="0a0cba998643670ca90316eada46fecac1a92e85"

export EMULTRF_DATA_URL="git@github.com:SBU-COSMOLIKE/emulators_data_lcdm.git"
#export EMULTRF_DATA_GIT_COMMIT="208fc465a2beb219fabf92e0802c2815d7556588"

export DARKEMULATOR_URL="https://github.com/DarkQuestCosmology/dark_emulator_public.git"
export DARKEMULATOR_GIT_COMMIT="46df5972509624e2eeadc2bf3ac528b02333a7e2"

export FASTPT_URL="https://github.com/jablazek/FAST-PT.git"
export FASTPT_GIT_COMMIT="5e65ad23becaaae5b18aedcaacab99411df92b0f"

# --------------------------------------------------------------------
# --------------------------------------------------------------------
# Cosmolike projects below -------------------------------------------
# --------------------------------------------------------------------
# --------------------------------------------------------------------
export LSST_Y1_URL="git@github.com:SBU-COSMOLIKE/cocoa_lsst_y1.git"
export LSST_Y1_NAME="lsst_y1"
#BRANCH: if unset, load the latest commit on the specified branch
#export LSST_Y1_BRANCH="dev"
#COMMIT: if unset, load the specified commit
#export LSST_Y1_COMMIT="df96af9558c97b07d355df4bfc56f1677e71b201"
#BRANCH: if unset, load the specified TAG
#export LSST_Y1_TAG="v4.0-beta17"

export DES_Y3_URL="git@github.com:SBU-COSMOLIKE/cocoa_des_y3.git"
export DES_Y3_NAME="des_y3"
#BRANCH: if unset, load the latest commit on the specified branch
#export DES_Y3_BRANCH="main"
#COMMIT: if unset, load the specified commit
#export DES_Y3_COMMIT="1a46582b5539c177bd68f8863c054f79a15f8538"
#BRANCH: if unset, load the specified TAG
#export DES_Y3_TAG="v4.0-beta17"

export ROMAN_FOURIER_URL="git@github.com:SBU-COSMOLIKE/cocoa_roman_fourier.git"
export ROMAN_FOURIER_NAME="roman_fourier"
#BRANCH: if unset, load the latest commit on the specified branch
#export ROMAN_FOURIER_BRANCH="main"
#COMMIT: if unset, load the specified commit
#export ROMAN_FOURIER_COMMIT="407a35a15b2a1d96d96cb5f0276cf772c2c60e6d"
#BRANCH: if unset, load the specified TAG
#export ROMAN_FOURIER_TAG="v4.0-beta17"

export ROMAN_REAL_URL="git@github.com:SBU-COSMOLIKE/cocoa_roman_real.git"
export ROMAN_REAL_NAME="roman_real"
#BRANCH: if unset, load the latest commit on the specified branch
#export ROMAN_REAL_BRANCH="main"
#COMMIT: if unset, load the specified commit
export ROMAN_REAL_COMMIT="8a13be52849fc7965b99f41bd173b7dda05fba67"
#BRANCH: if unset, load the specified TAG
#export ROMAN_REAL_TAG="v4.0-beta17"

# --------------------------------------------------------------------
# --------------------------------------------------------------------
# private projects below ---------------------------------------------
# --------------------------------------------------------------------
# --------------------------------------------------------------------
export AXIONS_PROJECT_URL="git@github.com:SBU-COSMOLIKE/cocoa_axions.git"
export AXIONS_PROJECT_NAME="axions"
#export AXIONS_PROJECT_COMMIT="df96af9558c97b07d355df4bfc56f1677e71b201"
