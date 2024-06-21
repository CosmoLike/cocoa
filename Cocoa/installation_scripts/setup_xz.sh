#!/bin/bash
# ----------------------------------------------------------------------------
# -------------------------- XZ COMPRESSION Library --------------------------
# ----------------------------------------------------------------------------
if [ -z "${IGNORE_XZ_INSTALLATION}" ]; then
  
  pfail() {
    echo -e "\033[0;31m ERROR ENV VARIABLE ${1} IS NOT DEFINED \033[0m"
    unset pfail
  }
  
  if [ -z "${ROOTDIR}" ]; then
    pfail 'ROOTDIR'; return 1
  fi

  cdroot() {
    cd "${ROOTDIR}" 2>"/dev/null" || { echo -e \
      "\033[0;31m\t\t CD ROOTDIR (${ROOTDIR}) FAILED \033[0m"; return 1; }
    unset cdroot
  }

  if [ -z "${C_COMPILER}" ]; then
    pfail 'C_COMPILER'; cdroot; return 1;
  fi

  unset_env_vars_sxz () {
    unset OUT1
    unset OUT2
    unset CCIL
    unset pfail
    unset unset_env_vars_sxz
    cdroot || return 1;
  }
  
  fail_sxz () {
    local MSG="\033[0;31m (setup_xz.sh) WE CANNOT RUN \e[3m"
    local MSG2="\033[0m"
    echo -e "${MSG} ${1} ${MSG2}"  
    unset fail_sxz
    unset_env_vars_sxz
  }

  cdfolder() {
    cd "${1}" 2>"/dev/null" || { fail_sxz "CD FOLDER: ${1}"; return 1; }
  }

  if [ -z "${DEBUG_XZ_PACKAGE}" ]; then
    if [ -z "${MAKE_NUM_THREADS}" ]; then
      pfail 'MAKE_NUM_THREADS'; cdroot; return 1;
    fi
    export OUT1="/dev/null"; export OUT2="/dev/null"
    export XZMNT="${MAKE_NUM_THREADS}"
  else
    export OUT1="/dev/tty"; export OUT2="/dev/tty"
    export XZMNT=1
  fi

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  
  ptop2 "SETUP_XZ"
  
  export CCIL="${ROOTDIR}/../cocoa_installation_libraries"

  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  ptop "INSTALLING XZ LIBRARY"

  cdfolder "${CCIL}" || return 1;

  #False xz file: just to trigger GIT LFS
  cp xz-5.2.5.tar.gz.xz xz-5.2.5.tar.gz 2>${OUT2} ||  
    { fail_sxz "CP XZ TAR"; return 1; }

  tar -xf xz-5.2.5.tar.gz.xz >${OUT1} 2>${OUT2} ||  
    { fail_sxz "TAR XZ TAR"; return 1; }

  cdfolder "${CCIL}/xz-5.2.5/" || return 1;

  CC=$C_COMPILER ./configure --prefix="${ROOTDIR}/.local" >${OUT1} 2>${OUT2} || 
    { fail_sxz "CONFIGURE"; return 1; }

  make -j $XZMNT all >${OUT1} 2>${OUT2} || { fail_sxz "MAKE"; return 1; }

  make install >${OUT1} 2>${OUT2} || { fail_sxz "MAKE INSTALL"; return 1; }

  cdfolder "${ROOTDIR}" || return 1;

  pbottom "INSTALLING XZ LIBRARY"
  
  # ----------------------------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  unset_env_vars_sxz || return 1;

  pbottom2 "SETUP_XZ"
fi

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------