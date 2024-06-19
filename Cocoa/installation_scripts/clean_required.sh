#!/bin/bash
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
pfail() {
  echo -e "\033[0;31m ERROR ENV VARIABLE ${1} IS NOT DEFINED \033[0m"
  unset pfail
}
if [ -z "${ROOTDIR}" ]; then
  pfail 'ROOTDIR'
  return 1
fi
unset_env_vars_clean_req () {
  unset SETUP_PREREQUISITE_DONE
  unset pfail
  unset unset_env_vars_clean_req
  cd "${ROOTDIR}" || return 1
}
fail_clean_clean_req () {
  local MSG="\033[0;31m (clean_required.sh) WE CANNOT RUN \e[3m"
  local MSG2="\033[0m"
  echo -e "${MSG} ${1} ${MSG2}"
  unset fail_clean_clean_req
  unset_env_vars_clean_req
  return 1
}

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
ptop 'CLEANING COCOA REQUIRED LIBRARIES'

cd "${ROOTDIR}"/../cocoa_installation_libraries 2> ${OUT2}
if [ $? -ne 0 ]; then
  fail_clean_clean_req "CD COCOA_INSTALLATION_LIBRARIES FOLDER"
  return 1
fi

sh clean_all
if [ $? -ne 0 ]; then
  fail_clean_clean_req "SCRIPT clean_all.sh"
fi

unset_env_vars_clean_req
pbottom 'CLEANING COCOA REQUIRED LIBRARIES'
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------