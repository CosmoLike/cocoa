#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${ROOTDIR}" ]; then
    echo 'ERROR ROOTDIR not defined'
    return
fi

( sh "${ROOTDIR:?}"/external_modules/data/clean_all.sh )

cd "${ROOTDIR:?}/external_modules/data"

( sh unxv_sn.sh ) || { cd "${ROOTDIR:?}"; return 1; }

( sh unxv_bao.sh ) || { cd "${ROOTDIR:?}"; return 1; }

( sh unxv_h0licow.sh ) || { cd "${ROOTDIR:?}"; return 1; }

( sh unxv_act_dr6.sh ) || { cd "${ROOTDIR:?}"; return 1; }

( sh unxv_simons_observatory.sh ) || { cd "${ROOTDIR:?}"; return 1; }

( sh unxv_bicep.sh ) || { cd "${ROOTDIR:?}"; return 1; }

( sh unxv_spt.sh ) || { cd "${ROOTDIR:?}"; return 1; }

( sh unxv_planck2018_basic.sh ) || { cd "${ROOTDIR:?}"; return 1; }

( sh unxv_camspec.sh ) || { cd "${ROOTDIR:?}"; return 1; }

( sh unxv_lipop.sh ) || { cd "${ROOTDIR:?}"; return 1; }

cd "${ROOTDIR:?}"

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------