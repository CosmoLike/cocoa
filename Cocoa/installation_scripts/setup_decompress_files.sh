# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
echo -e '\033[1;44m''SETUP_DECOMPRESS_FILES''\033[0m'

if [ -z "${ROOTDIR}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE ROOTDIR IS NOT DEFINED''\033[0m'
    return 1
fi
if [ -z "${CXX_COMPILER}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE CXX_COMPILER IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
fi
if [ -z "${C_COMPILER}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE C_COMPILER IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
fi
if [ -z "${FORTRAN_COMPILER}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE FORTRAN_COMPILER IS NOT DEFINED''\033[0m'
    cd $ROOTDIR
    return 1
fi
if [ -z "${CMAKE}" ]; then
    echo -e '\033[0;31m''ERROR ENV VARIABLE MAKE IS NOT DEFINED''\033[0m'
    return 1
fi

# ----------------------------------------------------------------------------
# ----------------- COCOA_INSTALLATION_LIBRARIES -----------------------------
# ----------------------------------------------------------------------------
if [ -z "${NO_UNXZ_COCOA_INSTALLATION_LIBRARIES}" ]; then
    echo -e '\033[1;34m''\t DECOMPRESSING FILES ON /COCOA_INSTALLATION_LIBRARIES''\033[0m'

    cd $ROOTDIR/../cocoa_installation_libraries/

    sh unxv_all.sh

    cd $ROOTDIR

    echo -e '\033[1;34m''\t \e[4mDECOMPRESSING FILES ON /COCOA_INSTALLATION_LIBRARIES DONE''\033[0m'
fi
# ----------------------------------------------------------------------------
# ----------------------- /EXTERNAL_MODULES/CODE -----------------------------
# ----------------------------------------------------------------------------
#if [ -z "${NO_UNXZ_EXTERNAL_MODULES_CODE}" ]; then
#    # no xz files under code so far
#fi

# ----------------------------------------------------------------------------
# ----------------------- /EXTERNAL_MODULES/DATA -----------------------------
# ----------------------------------------------------------------------------
if [ -z "${NO_UNXZ_EXTERNAL_MODULES_DATA}" ]; then
    echo -e '\033[1;34m''\t DECOMPRESSING FILES ON /EXTERNAL_MODULES/DATA''\033[0m'

    cd $ROOTDIR/external_modules/data/

    sh unxv_all.sh

    cd $ROOTDIR

    echo -e '\033[1;34m''\t \e[4mDECOMPRESSING FILES ON /EXTERNAL_MODULES/DATA DONE''\033[0m'
fi

echo -e '\033[1;44m''\e[4mSETUP_DECOMPRESS_FILES DONE''\033[0m'
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------