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
  cd $ROOTDIR
  unset SETUP_PREREQUISITE_DONE
  unset pfail
  unset unset_env_vars_clean_req
}
fail_clean_clean_req () {
  export MSG="\033[0;31m (clean_required.sh) WE CANNOT RUN \e[3m"
  export MSG2="\033[0m"
  echo -e "${MSG} ${1} ${MSG2}"
  unset_env_vars_clean_req
  unset MSG
  unset MSG2
  unset fail_clean_clean_req
}
# ---------------------------------------------------------------------------
ptop2 'COCOA REQUIRED LIBRARIES'

cd $ROOTDIR/../cocoa_installation_libraries

sh clean_all
if [ $? -ne 0 ]; then
  fail_clean_clean_req "clean_all.sh"
  return 1
fi

unset_env_vars_clean_req
pbottom2 'COCOA REQUIRED LIBRARIES'
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------