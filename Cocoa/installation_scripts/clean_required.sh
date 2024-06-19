#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
pfail() {
  echo -e "\033[0;31m ERROR ENV VARIABLE ${1} NOT DEFINED \033[0m"
  unset pfail
}
if [ -z "${ROOTDIR}" ]; then
  pfail 'ROOTDIR'
  return 1
fi
cdroot() {
  cd "${ROOTDIR}" 2>"/dev/null" || { echo -e \
    "\033[0;31m\t\t CD ROOTDIR (${ROOTDIR}) FAILED \033[0m"; return 1; }
  unset cdroot
}
unset_env_vars_clean_req () {
  unset SETUP_PREREQUISITE_DONE
  unset pfail
  unset unset_env_vars_clean_req
  cdroot || return 1;
}
fail_clean_req () {
  local MSG="\033[0;31m\t\t (clean_required.sh) WE CANNOT RUN \e[3m"
  local MSG2="\033[0m"
  echo -e "${MSG} ${1} ${MSG2}"
  unset fail_clean_req
  unset_env_vars_clean_req
  return 1
}
cdfolder() {
  cd "${1}" 2>"/dev/null" || { fail_clean_req "CD FOLDER: ${1}"; return 1; }
}

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
ptop 'CLEANING COCOA REQUIRED LIBRARIES'

cdfolder "${ROOTDIR}"/../cocoa_installation_libraries || return 1

sh clean_all || { fail_clean_req "SCRIPT clean_all.sh"; return 1; }

unset_env_vars_clean_req || return 1
pbottom 'CLEANING COCOA REQUIRED LIBRARIES'
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------