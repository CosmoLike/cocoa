#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -z "${ROOTDIR}" ]; then
  echo 'ERROR ROOTDIR not defined'
  return
fi

unset_env_vars () {
  cdroot || return 1;
}

unset_env_funcs () {
  unset -f cdfolder cpfolder error
  unset -f unset_env_funcs
  cdroot || return 1;
}

unset_all () {
  unset_env_vars
  unset_env_funcs
  unset -f unset_all
  cdroot || return 1;
}

error () {
  fail_script_msg "$(basename ${BASH_SOURCE[0]})" "${1}"
  unset_all || return 1
}

cdfolder() {
  cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

cdfolder || { error "CD ROOTDIR"; return 1; }

( sh unxv_sn.sh ) || { cdroot; return 1; }

( sh unxv_bao.sh ) || { cdroot; return 1; }

( sh unxv_h0licow.sh ) || { cdroot; return 1; }

( sh unxv_act_dr6.sh ) || { cdroot; return 1; }

( sh unxv_simons_observatory.sh ) || { cdroot; return 1; }

( sh unxv_bicep.sh ) || { cdroot; return 1; }

( sh unxv_spt.sh ) || { cdroot; return 1; }

( sh unxv_planck2018_basic.sh ) || { cdroot; return 1; }

( sh unxv_camspec.sh ) || { cdroot; return 1; }

( sh unxv_lipop.sh ) || { cdroot; return 1; }

# ---------------------------------------------------------------------------

cdfolder "${ROOTDIR}" || return 1

unset_all || return 1

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------