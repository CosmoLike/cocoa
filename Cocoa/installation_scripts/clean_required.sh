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
if [ -z "${ptop}" || -z "${pbottom}" ]; then
  pfail "PTOP/PBOTTOM"
  cd $ROOTDIR
  return 1
fi
unset_env_vars () {
  cd $ROOTDIR
  unset SETUP_PREREQUISITE_DONE
  unset unset_env_vars
}

# ---------------------------------------------------------------------------
ptop 'COCOA REQUIRED LIBRARIES'

cd $ROOTDIR/../cocoa_installation_libraries

sh clean_all
if [ $? -ne 0 ]; then
  unset_env_vars
  return 1
fi

unset_env_vars
pbottom 'COCOA REQUIRED LIBRARIES'
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------