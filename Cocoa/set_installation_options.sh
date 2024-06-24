#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
source ./installation_scripts/impl_unset_keys.sh
export ROOTDIR=$(pwd -P)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ----------------------- HOW COCOA SHOULD BE INSTALLED? -----------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
export MINICONDA_INSTALLATION=1
#export MANUAL_INSTALLATION=1

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# --------------------------- VERBOSE AS DEBUG TOOL ----------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#export COCOA_OUTPUT_VERBOSE=1

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# --- SKIP DOWNLOAD/DECOMPRESSING DATASET (SAVE TIME WHEN INSTALLING COCOA) ----
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# export SKIP_DECOMM_ACT
# export SKIP_DECOMM_SPT
# export SKIP_DECOMM_PLANCK
# export SKIP_DECOMM_BICEP
# export SKIP_DECOMM_STRONG_LENSING
# export SKIP_DECOMM_SN
# export SKIP_DECOMM_BAO
export SKIP_DECOMM_SIMONS_OBSERVATORY=1
export SKIP_DECOMM_CAMSPEC=1
export SKIP_DECOMM_LIPOP=1

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ---- IF DEFINED, COSMOLIKE WILL BE COMPILED WITH DEBUG FLAG ------------------
# ---- DEBUG FLAG = ALL COMPILER WARNINGS + NO MATH OPTIMIZATION + NO OPENMP ---
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
export COSMOLIKE_DEBUG_MODE=1

# ------------------------------------------------------------------------------
# ----- IF TRUE, THEN COCOA USES CLIK FROM https://github.com/benabed/clik -----
# ------------------------------------------------------------------------------
export USE_SPT_CLIK_PLANCK=1

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# --------------- CONTROL OVER THE COMPILATION OF EXTERNAL CODES ---------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#export IGNORE_CAMB_COMPILATION=1
#export IGNORE_CLASS_COMPILATION=1
#export IGNORE_COSMOLIKE_COMPILATION=1
#export IGNORE_POLYCHORD_COMPILATION=1
#export IGNORE_PLANCK_COMPILATION=1
#export IGNORE_ACT_COMPILATION=1
#export IGNORE_COBAYA_INSTALLATION=1
#export IGNORE_CAMSPEC_INSTALLATION=1
#export IGNORE_LIPOP_INSTALLATION=1
#export IGNORE_SO_INSTALLATION=1

# ------------------------------------------------------------------------------
# THREADING COMPILATION/INSTALLATION OF LIBRARIES ------------------------------
# ------------------------------------------------------------------------------
export MAKE_NUM_THREADS=4

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# IF NOT SET, COCOA WILL INSTALL TENSORFLOW, KERAS, AND PYTORCH 
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
export IGNORE_EMULATOR_CPU_PIP_PACKAGES=1
export IGNORE_EMULATOR_GPU_PIP_PACKAGES=1

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ---------------------- DERIVED & RARELY USED FLAGS ---------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

if [ -n "${MANUAL_INSTALLATION}" ]; then      
  
  # NEED TO ADJUST THE FLAGS IN THE FILE BELOW IF MANUAL INSTALLATION  
  source .flags_manual_installation.sh || return 1

elif [ -n "${MINICONDA_INSTALLATION}" ]; then

  source .flags_miniconda_installation.sh || return 1

fi

source .flags_derived.sh || return 1

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ----- PACKAGE URL AND VERSIONS. CHANGES IN THE COMMIT ID MAY BREAK COCOA -----
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

export COBAYA_URL="https://github.com/CobayaSampler/cobaya.git"
export COBAYA_GIT_COMMIT="2636ea9ed399c35c5d276de1acb15aaafbcab10c"

export HILLIPOP_URL="https://github.com/planck-npipe/hillipop.git"
export HILLIPOP_GIT_COMMIT="cc9cbe31991d4662522241543a46d44d2cdec251"

export LOLLIPOP_URL="https://github.com/planck-npipe/lollipop.git"
export LOLLIPOP_GIT_COMMIT="280a9c93d33bc6a058d6bf769ec82d9f7fdbd2b6"

export LIPOP_DATA_VERSION=4.2

export SPT3G_DATA_GIT_COMMIT="66da8e9e2f325024566fe13245788bf8ede897bc"

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

export ACT_NAME="pyactlike"

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------