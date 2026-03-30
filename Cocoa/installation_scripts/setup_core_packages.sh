#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if [ -n "${IGNORE_CORE_INSTALLATION:-}" ]; then
  return 99
fi

if [ -z "${ROOTDIR:-}" ]; then
  pfail 'ROOTDIR'; return 1
fi

# parenthesis = run in a subshell
( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" )  || return 1;

unset_env_vars () {
  unset -v URL_BASE URL FOLDER VER XZF CCIL CNAME PACKDIR PACKAGE_VERSION 
  cdroot || return 1;
}

unset_env_funcs () {
  unset -f cdfolder cpfolder cpfile error wgetact gitact gitact1 gitact2
  unset -f unset_env_funcs wgetact1 wgetact2
  cdroot || return 1;
}

unset_all () {
  unset_env_vars
  unset_env_funcs
  unset -f unset_all
  cdroot || return 1;
}

error () {
  fail_script_msg "$(basename "${BASH_SOURCE[0]}")" "${1}"
  unset_all || return 1
}

cdfolder() {
  cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
}

cpfolder() {
  cp -r "${1:?}" "${2:?}"  \
    2>"/dev/null" || { error "CP FOLDER ${1} on ${2}"; return 1; }
}

cpfile() {
  cp "${1:?}" "${2:?}" \
    2>"/dev/null" || { error "CP FILE ${1} on ${2}"; return 1; }
}

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Logic of these functions: Cocoa processes core packages using three scripts
#
# (1) setup_core_packages: download core packages and save them in the appropriate
#                          compressed format. What is the appropriate compressed format? .xz
#
# (2) unxv_core_packages: uncompress the compressed core packages (assuming
#                         appropriate format, i.e., “.xz”).
#
# (3) compile_core_packages: compile core packages
#-----------------------------------------------------------------------------
#Given that,
#
#wgetact1: This function downloads core packages and saves them in the
#          appropriate folder (uncompressed) format. The name of the uncompressed 
#          folder is chosen by the user via a set of keys. Example of
#          such keys: COCOA_CMAKE_DIR and COCOA_WGET_DIR.
#          These keys all have reasonable default values, so setting them is optional.
#
#wgetact2: compress the uncompressed folder to .xz to be later uncompressed again
#          by unxv_core_packages (we know is pedantic, but for historical reasons
#          we wanted the setup / unxv / compile scripts to have well defined roles
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
wgetact1() {
  #ARGUMENTS: ORIGINAL COMPRESSED FILE ROOT (1)
  #           ORIGINAL COMPRESSED FILE TYPE (2) 
  #           URL WHERE TO DOWNLOAD THE PACKAGE (3) 
  #           FINAL UNCOMPRESS DIRECTORY WITH APPROPRIATE NAME CHOSEN BY USER(4) 

  local CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"
  
  local FILE="${1:?}.${2:?}" # Full archive filename, e.g. "fftw.tar.gz"

  # --------------------------------------------------------------------------
  # In case this script runs twice -------------------------------------------
  # --------------------------------------------------------------------------
  if [ -n "${OVERWRITE_EXISTING_CORE_PACKAGES:-}" ]; then
  
    rm -rf "${CCIL:?}/${1:?}" # INITIAL DIRECTORY AFTER DECOMPRESSION
    rm -rf "${CCIL:?}/${4:?}" # FINAL DIRECTORY WHERE THE LIB WILL BE INSTALLED
    
    if [ -n "${REDOWNLOAD_EXISTING_CORE_PACKAGES:-}" ]; then
    
      rm -f  "${CCIL:?}/${FILE:?}" # COMPRESSED FILE
    
    fi 
  
  fi

  cdfolder "${CCIL:?}" || { unset_all; return 1; }
  
  if [[ ! -d "${CCIL:?}/${4:?}" && ! -d "${CCIL:?}/${1:?}" ]]; then

    if [[ ! -e "${CCIL:?}/${FILE:?}" ]]; then
    
      "${WGET:?}" "${3:?}" -q --show-progress --no-check-certificate \
        --progress=bar:force:noscroll --timeout=30 --tries=2 --waitretry=0 \
        --retry-connrefused --read-timeout=30 \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC24:?}"; return 1; }

    fi

    # ------------------------------------------------------------------------
    # Check the archive layout before uncompressing the downloaded file
    # ------------------------------------------------------------------------
    local EXPECTED="${1:?}" # Expected top-level directory after extraction
    local FOUND=0           # Will become 1 if we find at least one valid entry
    if [ "${2:?}" = "tar.gz" ]; then
    
      # List the contents of the .tar.gz archive
      while IFS= read -r entry; do
        entry="${entry#./}" # Remove a possible leading "./"
        [ -z "${entry}" ] && continue

        case "${entry}" in
          # Accept: EXPECTED, EXPECTED/, EXPECTED/anything/else
          "${EXPECTED}"|"${EXPECTED}/"|"${EXPECTED}/"*)
            FOUND=1
            ;;
          *)
            error "${EC25:?} (unexpected archive layout)"
            return 1
            ;;
        esac
      done < <(tar tzf "${FILE}" 2>>${OUT2:?})
      # explanation of the < <(tar tzf "${FILE}" 2>>${OUT2:?})
      # The left  < operator: is the normal input redirection operator 
      #                       i.e., take stdin for this loop from the thing on my right)
      #                       i.e., < [some temporary file-like stream created by <(...)]
      # The right <(...) operator: process substitution that allows you to run 
      #                            the bash program tar tzf "${FILE}" 2>>${OUT2:?}
      #                            and redirect its output to here.

    elif [ "${2:?}" = "tar.xz" ]; then

      while IFS= read -r entry; do
        entry="${entry#./}"
        [ -z "${entry}" ] && continue
        case "${entry}" in
          "${EXPECTED}"|"${EXPECTED}/"|"${EXPECTED}/"*)
            FOUND=1
            ;;
          *)
            error "${EC25:?} (unexpected archive layout)"
            return 1
            ;;
        esac
      done < <(tar tf "${FILE}" 2>>${OUT2:?})
    
    else
    
      error "UNKNOWN FILE EXTENSION (validation)"; return 1
    
    fi

    if [ "${FOUND}" -eq 0 ]; then
    
      error "${EC25:?} (empty archive)"; return 1
    
    fi

    # ------------------------------------------------------------------------
    # finally: uncompress the file
    # ------------------------------------------------------------------------

    if [ "${2:?}" == "tar.gz" ]; then
    
      tar zxf "${FILE}" \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC25:?} (gz)"; return 1; }
    
    elif [ "${2:?}" == "tar.xz" ]; then
    
      tar xf "${FILE}" \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC25:?} (xz)"; return 1; }
    
    else
    
      error "UNKNOWN FILE EXTENSION"; return 1;
    
    fi

    if [[ "${1:?}/" != "${4:?}" && \
          "${1:?}"  != "${4:?}" && \
          "${1:?}"  != "${4:?}/" ]]; then
      
      # In case this script runs twice (e.g., after being killed by CTRL-D)
      rm -rf "${CCIL:?}/${4:?}" 
      
      mv "${1:?}/" "${4:?}" \
        2>>${OUT2:?} || { error "MV FOLDER"; return 1; }
    
    fi
  
  fi

  cdfolder "${ROOTDIR:?}" || { unset_all; return 1; }
}

wgetact2() {
  #ARGUMENTS: FINAL UNCOMPRESS DIRECTORY WITH APPROPRIATE NAME CHOSEN BY USER(4) 
  #           FINAL COMPRESSED DIRECTORY IN THE APPROPRIATE XZ FORMAT (5)
  
  local CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

  cdfolder "${CCIL:?}" || { unset_all; return 1; }

  # --------------------------------------------------------------------------
  # In case this script runs twice -------------------------------------------
  # --------------------------------------------------------------------------
  if [ -n "${OVERWRITE_EXISTING_CORE_PACKAGES}" ]; then
  
    rm -f  "${CCIL:?}/${2:?}"
  
  fi

  if [ ! -e "${CCIL:?}/${2:?}" ]; then
  
    # Why this compress? In the old Cocoa we saved the libraries in the github repo 
    #                    using git lfs. So this way preserves the old scripts
    #                    We set "-k 1" to be the minimum compression. As we said, 
    #                    xz is used just for historical compatibility)
    (
      set -o pipefail
      tar -cf - "${1:?}" 2>>"${OUT2:?}" | 
      xz -k -1 --threads="${MNT:-1}" -c - 2>>"${OUT2:?}" > "${2:?}" \
        2>>${OUT2:?} 
    ) || { error "TAR (COMPRESS)"; return 1; }
  
  fi

  cdfolder "${ROOTDIR:?}" || { unset_all; return 1; }
}

wgetact() {
  #ARGUMENTS: ORIGINAL COMPRESSED FILE ROOT (1)
  #           ORIGINAL COMPRESSED FILE TYPE (2) 
  #           URL WHERE TO DOWNLOAD THE PACKAGE (3) 
  #           FINAL UNCOMPRESS DIRECTORY WITH APPROPRIATE NAME CHOSEN BY USER(4) 
  #           FINAL COMPRESSED DIRECTORY IN THE APPROPRIATE XZ FORMAT (5)

  wgetact1 "${1:?}" "${2:?}" "${3:?}" "${4:?}" || { unset_all; return 1; }

  wgetact2 "${4:?}" "${5:?}" || { unset_all; return 1; }

}

gitact1() {
  #ARGUMENTS: UNCOMPRESS DIRECTORY CONTAINING THE GIT REPO (1)
  #           GIT VERSION (TAG) (2) 
  #           URL WHERE TO DOWNLOAD THE PACKAGE (3) 
  
  local CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

  cdfolder "${CCIL:?}" || { unset_all; return 1; }

  # --------------------------------------------------------------------------
  # In case this script runs twice -------------------------------------------
  # --------------------------------------------------------------------------
  if [ -n "${OVERWRITE_EXISTING_CORE_PACKAGES:-}" ]; then
  
    rm -rf "${CCIL:?}/${1:?}"
  
  fi

  if [ ! -d "${CCIL:?}/${1:?}" ]; then
  
    "${GIT:?}" clone "${3:?}" "${1:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "GIT CLONE"; return 1; }
    
    cdfolder "${CCIL:?}/${1:?}" || return 1;
    
    "${GIT:?}" checkout "${2:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "GIT CHECKOUT"; return 1; }
    
    rm -rf "${CCIL:?}/${1:?}/.git/"
  
  fi

  cdfolder "${ROOTDIR:?}" || { unset_all; return 1; }
}

gitact2() {
  #ARGUMENTS: UNCOMPRESS DIRECTORY (1)
  #           FINAL COMPRESSED DIRECTORY IN THE APPROPRIATE XZ FORMAT (2)
  
  local CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

  # --------------------------------------------------------------------------
  # In case this script runs twice -------------------------------------------
  # --------------------------------------------------------------------------
  if [ -n "${OVERWRITE_EXISTING_CORE_PACKAGES}" ]; then
  
    rm -f  "${CCIL:?}/${2:?}"
  
  fi

  if [ ! -f "${CCIL:?}/${2:?}" ]; then
  
    cdfolder "${CCIL:?}" || return 1;

    # Why this compress? In the old Cocoa we saved the libraries in the github repo 
    #                    using git lfs. So this way preserves the old scripts
    #                    We set "-k 1" to be the minimum compression. As we said, 
    #                    xz is used just for historical compatibility)
    (
      set -o pipefail
      tar -cf - "${1:?}" 2>>${OUT2:?}  | 
      xz -k -1 --threads="${MNT:-1}" -c - 2>>${OUT2:?}  > "${2:?}"
    ) || { error "TAR (COMPRESS)"; return 1; }


  fi

  cdfolder "${ROOTDIR:?}" || { unset_all; return 1; }
}

gitact() {
  #ARGUMENTS: UNCOMPRESS DIRECTORY CONTAINING THE GIT REPO (1)
  #           GIT VERSION (TAG) (2) 
  #           URL WHERE TO DOWNLOAD THE PACKAGE (3) 
  #           FINAL COMPRESSED DIRECTORY IN THE APPROPRIATE XZ FORMAT (4)
  gitact1 "${1:?}" "${2:?}" "${3:?}" || return 1;
  
  gitact2 "${1:?}" "${4:?}" || return 1;
}

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
 
unset_env_vars || return 1;

CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

# ----------------------------------------------------------------------------
# ------------------------------ XZ LIBRARY ----------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_XZ_INSTALLATION}" ]; then

  ptop "GETTING AND COMPILING XZ LIBRARY (CORE LIBS)" || { unset_all; return 1; }

  cdfolder "${CCIL:?}" || { unset_all; return 1; }

  #False xz file (if I am installing xz - this has to be a fake xz file): just to trigger GIT LFS
  cp xz-5.2.5.tar.gz.xz xz-5.2.5.tar.gz \
    2>>${OUT2} ||  { error "CP XZ TAR"; return 1; }

  tar -xf xz-5.2.5.tar.gz.xz \
    >>${OUT1:?} 2>>${OUT2:?} ||  { error "TAR XZ TAR"; return 1; }

  cdfolder "${CCIL:?}/xz-5.2.5/" || { unset_all; return 1; }

  CC="${C_COMPILER:?}" ./configure --prefix="${ROOTDIR:?}/.local" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC11:?}"; return 1; }

  make -j "${MNT:-1}" all >>${OUT1:?} 2>>${OUT2:?} || { error "${EC8:?}"; return 1; }

  make install >>${OUT1:?} 2>>${OUT2:?} || { error "${EC10:?}"; return 1; }

  cdfolder "${ROOTDIR}" || return 1;

  pbottom "GETTING AND COMPILING XZ LIBRARY (CORE LIBS)" || { unset_all; return 1; }

fi

# ----------------------------------------------------------------------------
# ------------------------------ CMAKE LIBRARY --------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_CMAKE_INSTALLATION}" ]; then
  
  ptop "GETTING CMAKE LIBRARY (CORE LIBS)" || { unset_all; return 1; }

  URL='https://github.com/Kitware/CMake.git'
  
  PACKAGE_VERSION="${CMAKE_VERSION:-"3.26.4"}"

  FOLDER="cmake-${PACKAGE_VERSION:?}"
  
  VER=v${PACKAGE_VERSION:?}

  # note: Do not change XFZ filename. 
  # note: Otherwise the script unxv_core_packages.sh will stop working
  XZF="cmake.xz"

  gitact "${COCOA_CMAKE_DIR:-"${FOLDER:?}"}" \
         "${VER:?}" \
         "${URL:?}" \
         "${XZF:?}" || { unset_all; return 1; }

  unset -v URL FOLDER VER XZF PACKAGE_VERSION

  # COMPILING CMAKE

  pbottom "GETTING CMAKE LIBRARY (CORE LIBS)" || { unset_all; return 1; }

fi

# ----------------------------------------------------------------------------
# ------------------------------ WGET LIBRARY --------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_WGET_INSTALLATION}" ]; then
  
  ptop "GETTING WGET LIBRARY (CORE LIBS)" || { unset_all; return 1; }

  PACKAGE_VERSION="${WGET_VERSION:-"1.24.5"}"

  FILE_ROOT="wget-${PACKAGE_VERSION:?}"

  COMPRESSED_FILE_EXT="tar.gz"

  URL="https://ftp.gnu.org/gnu/wget/${FILE_ROOT:?}.${COMPRESSED_FILE_EXT}"

  # note: Do not change XFZ filename. 
  # note: Otherwise the script unxv_core_packages.sh will stop working
  XZF="wget.xz"

  # note: Circular dependence (wget depends on wget!). Why? 
  # note: The wget commands show progress bar when downloading large likelihoods
  # note: The specific progress bar option only exists on recent wget versions  
  wgetact "${FILE_ROOT:?}" \
          "${COMPRESSED_FILE_EXT:?}" \
          "${URL:?}" \
          "${COCOA_WGET_DIR:-"${FILE_ROOT:?}"}" \
          "${XZF:?}" || { unset_all; return 1; }

  unset -v URL FILE_ROOT COMPRESSED_FILE_EXT XZF PACKAGE_VERSION

  pbottom "GETTING WGET LIBRARY (CORE LIBS)" || { unset_all; return 1; }

fi

# ----------------------------------------------------------------------------
# ------------------------------ BinUtils LIBRARY ----------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_DISTUTILS_INSTALLATION}" ]; then

  ptop "GETTING BINUTILS LIBRARY (CORE LIBS)" || { unset_all; return 1; }

  PACKAGE_VERSION="${BINUTILS_VERSION:-"2.37"}"

  FILE_ROOT="binutils-${PACKAGE_VERSION:?}"
  
  COMPRESSED_FILE_EXT="tar.gz"
  
  URL="https://ftp.gnu.org/gnu/binutils/${FILE_ROOT:?}.${COMPRESSED_FILE_EXT:?}"
  
  # note: Do not change XFZ filename. 
  # note: Otherwise the script unxv_core_packages.sh will stop working
  XZF="binutils.xz"

  wgetact "${FILE_ROOT:?}" \
          "${COMPRESSED_FILE_EXT:?}" "${URL:?}" \
          "${COCOA_BINUTILS_DIR:-"${FILE_ROOT:?}"}" \
          "${XZF:?}" || { unset_all; return 1; }

  unset -v URL FILE_ROOT COMPRESSED_FILE_EXT XZF PACKAGE_VERSION

  pbottom "GETTING BINUTILS LIBRARY (CORE LIBS)" || { unset_all; return 1; }

fi

# ----------------------------------------------------------------------------
# ------------------------------ TEXTINFO LIBRARY ----------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_DISTUTILS_INSTALLATION}" ]; then
  
  ptop "GETTING TEXINFO LIBRARY (CORE LIBS)" || { unset_all; return 1; }

  PACKAGE_VERSION="${TEXINFO_VERSION:-"7.0.3"}"

  FILE_ROOT="texinfo-${PACKAGE_VERSION:?}"
  
  COMPRESSED_FILE_EXT="tar.xz"
  
  URL="https://ftp.gnu.org/gnu/texinfo/${FILE_ROOT}.${COMPRESSED_FILE_EXT}"
 
  # note: Do not change XFZ filename. 
  # note: Otherwise the script unxv_core_packages.sh will stop working 
  XZF="texinfo.xz"

  wgetact "${FILE_ROOT:?}" \
          "${COMPRESSED_FILE_EXT:?}" \
          "${URL:?}" \
          "${COCOA_TEXINFO_DIR:-"${FILE_ROOT:?}"}" \
          "${XZF:?}" || { unset_all; return 1; }

  unset -v URL FILE_ROOT COMPRESSED_FILE_EXT XZF PACKAGE_VERSION

  pbottom "GETTING TEXINFO LIBRARY (CORE LIBS)" || { unset_all; return 1; }

fi

# ----------------------------------------------------------------------------
# ------------------------------ OpenBLAS LIBRARY ----------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_OPENBLAS_INSTALLATION}" ]; then

  ptop "GETTING OPENBLAS LIBRARY (CORE LIBS)" || { unset_all; return 1; }

  PACKAGE_VERSION="${OPENBLAS_VERSION:-"0.3.23"}" 

  FOLDER="OpenBLAS-${PACKAGE_VERSION:?}"

  URL='https://github.com/OpenMathLib/OpenBLAS.git'

  VER="v${PACKAGE_VERSION:?}"

  # note: Do not change XFZ filename. 
  # note: Otherwise the script unxv_core_packages.sh will stop working
  XZF="OpenBLAS.xz"

  gitact "${COCOA_OPENBLAS_DIR:-"${FOLDER:?}"}" \
         "${VER:?}" \
         "${URL:?}" \
         "${XZF:?}" || { unset_all; return 1; }

  unset -v URL FOLDER VER XZF PACKAGE_VERSION

  pbottom "GETTING OPENBLAS LIBRARY (CORE LIBS)" || { unset_all; return 1; }

fi

# ----------------------------------------------------------------------------
# ------------------------------ LAPACK LIBRARY ------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_FORTRAN_LAPACK_INSTALLATION}" ]; then

  ptop "GETTING LAPACK LIBRARY (CORE LIBS)" || { unset_all; return 1; }

  PACKAGE_VERSION="${LAPACK_VERSION:-"3.11"}" 

  FOLDER="lapack-${PACKAGE_VERSION:?}"

  URL='https://github.com/Reference-LAPACK/lapack.git'

  VER=v"${PACKAGE_VERSION:?}"

  # note: Do not change XFZ filename. 
  # note: Otherwise the script unxv_core_packages.sh will stop working
  XZF="lapack.xz"

  gitact "${COCOA_LAPACK_DIR:-"${FOLDER:?}"}" \
         "${VER:?}" \
         "${URL:?}" \
         "${XZF:?}" || { unset_all; return 1; }

  unset -v URL FOLDER VER XZF PACKAGE_VERSION

  pbottom "GETTING LAPACK LIBRARY DONE (CORE LIBS)" || { unset_all; return 1; }

fi

# ----------------------------------------------------------------------------
# ------------------------------ HDF5 LIBRARY --------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_HDF5_INSTALLATION}" ]; then

  ptop "GETTING HDF5 LIBRARY (CORE LIBS)" || { unset_all; return 1; }

  FILE_ROOT="hdf5-1.12.3"

  COMPRESSED_FILE_EXT="tar.gz"

  URL_BASE="https://hdf-wordpress-1.s3.amazonaws.com/wp-content/uploads"
  
  URL="${URL_BASE:?}/manual/HDF5/HDF5_1_12_3/src/${FILE_ROOT:?}.${COMPRESSED_FILE_EXT:?}"

  # note: Do not change XFZ filename. 
  # note: Otherwise the script unxv_core_packages.sh will stop working
  XZF="hdf5.xz"

  wgetact "${FILE_ROOT:?}" \
          "${COMPRESSED_FILE_EXT:?}" \
          "${URL:?}" \
          "${COCOA_HDF5_DIR:-"${FILE_ROOT:?}"}" \
          "${XZF:?}" || { unset_all; return 1; }

  unset -v URL_BASE URL FILE_ROOT COMPRESSED_FILE_EXT XZF

  pbottom "GETTING HDF5 LIBRARY DONE (CORE LIBS)" || { unset_all; return 1; }

fi

# ----------------------------------------------------------------------------
# ------------------------------ CFITSIO LIBRARY -----------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_C_CFITSIO_INSTALLATION}" ]; then
  
  ptop "GETTING CFITSIO LIBRARY (CORE LIBS)" || { unset_all; return 1; }

  PACKAGE_VERSION="${CFITSIO_VERSION:-"4.0.0"}" 

  FILE_ROOT="cfitsio-${PACKAGE_VERSION:?}"
  
  COMPRESSED_FILE_EXT="tar.gz"
  
  URL_BASE="http://heasarc.gsfc.nasa.gov"
  
  URL="${URL_BASE}/FTP/software/fitsio/c/${FILE_ROOT:?}.${COMPRESSED_FILE_EXT:?}"
 
  # note: Do not change XFZ filename. 
  # note: Otherwise the script unxv_core_packages.sh will stop working 
  XZF="cfitsio.xz"

  wgetact "${FILE_ROOT:?}" \
          "${COMPRESSED_FILE_EXT:?}" \
          "${URL:?}" \
          "${COCOA_CFITSIO_DIR:-"${FILE_ROOT:?}"}" \
          "${XZF:?}" || { unset_all; return 1; }

  unset -v URL_BASE URL FILE_ROOT COMPRESSED_FILE_EXT XZF PACKAGE_VERSION

  pbottom "GETTING CFITSIO LIBRARY DONE (CORE LIBS)" || { unset_all; return 1; }

fi

# ----------------------------------------------------------------------------
# ------------------------------ FFTW LIBRARY --------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_C_FFTW_INSTALLATION}" ]; then

  ptop "GETTING FFTW LIBRARY (CORE LIBS)" || { unset_all; return 1; }

  PACKAGE_VERSION="${FFTW_VERSION:-"3.3.10"}" 

  FILE_ROOT="fftw-${PACKAGE_VERSION:?}"

  COMPRESSED_FILE_EXT="tar.gz"

  URL="http://www.fftw.org/${FILE_ROOT:?}.${COMPRESSED_FILE_EXT:?}"

  # note: Do not change XFZ filename. 
  # note: Otherwise the script unxv_core_packages.sh will stop working
  XZF="fftw.xz"

  wgetact "${FILE_ROOT:?}" \
          "${COMPRESSED_FILE_EXT:?}" \
          "${URL:?}" \
          "${COCOA_FFTW_DIR:-"${FILE_ROOT:?}"}" \
          "${XZF:?}" || { unset_all; return 1; }

  unset -v URL FILE_ROOT COMPRESSED_FILE_EXT XZF PACKAGE_VERSION

  pbottom "GETTING FFTW LIBRARY DONE (CORE LIBS)" || { unset_all; return 1; }

fi

# ----------------------------------------------------------------------------
# ------------------------------ GSL LIBRARY ---------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_C_GSL_INSTALLATION}" ]; then

  ptop "GETTING GSL LIBRARY (CORE LIBS)" || { unset_all; return 1; }
  
  PACKAGE_VERSION="${GSL_VERSION:-"2.7"}" 

  FILE_ROOT="gsl-${PACKAGE_VERSION:?}"

  COMPRESSED_FILE_EXT="tar.gz"

  URL="http://ftp.wayne.edu/gnu/gsl/${FILE_ROOT:?}.${COMPRESSED_FILE_EXT:?}"

  # note: Do not change XFZ filename. 
  # note: Otherwise the script unxv_core_packages.sh will stop working
  XZF="gsl.xz"

  wgetact "${FILE_ROOT:?}" \
          "${COMPRESSED_FILE_EXT:?}" \
          "${URL:?}" \
          "${COCOA_GSL_DIR:-"${FILE_ROOT:?}"}" \
          "${XZF:?}" || { unset_all; return 1; }

  unset -v URL FILE_ROOT COMPRESSED_FILE_EXT XZF PACKAGE_VERSION

  pbottom "GETTING GSL LIBRARY DONE (CORE LIBS)" || { unset_all; return 1; }

fi

# ----------------------------------------------------------------------------
# ------------------------------ SPDLOG LIBRARY ------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_CPP_SPDLOG_INSTALLATION}" ]; then

  ptop "GETTING SPDLOG LIBRARY (CORE LIBS)" || { unset_all; return 1; }

  PACKAGE_VERSION="${SPDLOG_VERSION:-"v1.15.3"}"  

  URL='https://github.com/gabime/spdlog.git'

  FOLDER='spdlog'

  VER="${PACKAGE_VERSION:?}"

  # note: Do not change XFZ filename. 
  # note: Otherwise the script unxv_core_packages.sh will stop working
  XZF="spdlog.xz"

  gitact "${COCOA_SPDLOG_DIR:-"spdlog"}" \
         "${VER:?}" \
         "${URL:?}" \
         "${XZF:?}" || { unset_all; return 1; }

  unset -v URL FOLDER VER XZF

  pbottom "GETTING SPDLOG LIBRARY (CORE LIBS)" || { unset_all; return 1; }

fi

# ----------------------------------------------------------------------------
# ------------------------------ ARMA LIBRARY --------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_CPP_ARMA_INSTALLATION}" ]; then
  
  ptop "GETTING ARMA LIBRARY (CORE LIBS)" || { unset_all; return 1; }

  PACKAGE_VERSION="${ARMA_VERSION:-"12.8.4"}" 
  
  FILE_ROOT="armadillo-${PACKAGE_VERSION:?}"
  
  COMPRESSED_FILE_EXT="tar.xz"
  
  URL="https://sourceforge.net/projects/arma/files/${FILE_ROOT:?}.${COMPRESSED_FILE_EXT:?}"

  # note: Do not change XFZ filename. 
  # note: Otherwise the script unxv_core_packages.sh will stop working
  XZF="armadillo.xz"

  wgetact "${FILE_ROOT:?}" \
          "${COMPRESSED_FILE_EXT:?}" \
          "${URL:?}" \
          "${COCOA_ARMADILLO_DIR:-"${FILE_ROOT:?}"}" \
          "${XZF:?}" || { unset_all; return 1; }

  unset -v URL FILE_ROOT COMPRESSED_FILE_EXT XZF

  pbottom "GETTING ARMA LIBRARY DONE (CORE LIBS)" || { unset_all; return 1; }

fi

# ----------------------------------------------------------------------------
# ------------------------------ BOOST LIBRARY -------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_CPP_BOOST_INSTALLATION}" ]; then

  ptop "GETTING BOOST LIBRARY (CORE LIBS)" || { unset_all; return 1; }

  PACKAGE_VERSION="${CPP_BOOST_VERSION:-"81"}"

  FILE_ROOT="boost_1_${PACKAGE_VERSION:?}_0"

  COMPRESSED_FILE_EXT="tar.gz"

  URL_BASE="https://boostorg.jfrog.io/artifactory/main/release"
  
  URL="${URL_BASE}/1.${PACKAGE_VERSION:?}.0/source/${FILE_ROOT:?}.${COMPRESSED_FILE_EXT:?}"

  # note: Do not change XFZ filename. 
  # note: Otherwise the script unxv_core_packages.sh will stop working
  XZF="boost.xz"

  wgetact "${FILE_ROOT:?}" \
          "${COMPRESSED_FILE_EXT:?}" \
          "${URL:?}" \
          "${COCOA_BOOST_DIR:-"${FILE_ROOT:?}"}" \
          "${XZF:?}" || { unset_all; return 1; }

  unset -v URL_BASE URL FILE_ROOT COMPRESSED_FILE_EXT XZF PACKAGE_VERSION

  pbottom "GETTING BOOST LIBRARY (CORE LIBS)" || { unset_all; return 1; }

fi

# ----------------------------------------------------------------------------
# ------------------------------ CUBA LIBRARY --------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_CPP_CUBA_INSTALLATION}" ]; then

  ptop "GETTING CUBA LIBRARY (CORE LIBS)" || { unset_all; return 1; }

  PACKAGE_VERSION="${CPP_CUBA_VERSION:-"4.2.2"}"

  FILE_ROOT="Cuba-${PACKAGE_VERSION:?}"

  COMPRESSED_FILE_EXT="tar.gz"
  
  URL="https://feynarts.de/cuba/${FILE_ROOT:?}.${COMPRESSED_FILE_EXT:?}"

  # note: Do not change XFZ filename. 
  # note: Otherwise the script unxv_core_packages.sh will stop working
  XZF="cuba.xz"

  wgetact "${FILE_ROOT:?}" \
          "${COMPRESSED_FILE_EXT:?}" \
          "${URL:?}" \
          "${COCOA_CUBA_DIR:-"${FILE_ROOT:?}"}" \
          "${XZF:?}" || { unset_all; return 1; }

  unset -v URL FILE_ROOT COMPRESSED_FILE_EXT XZF PACKAGE_VERSION PACKDIR

  pbottom "GETTING CUBA LIBRARY (CORE LIBS)" || { unset_all; return 1; }

fi

# ----------------------------------------------------------------------------
# ------------------------------ CARMA LIBRARY -------------------------------
# ----------------------------------------------------------------------------

if [ -z "${IGNORE_CPP_CARMA_INSTALLATION}" ]; then

  ptop "GETTING CARMA LIBRARY (CORE LIBS)" || { unset_all; return 1; }

  CNAME="${COCOA_CARMA_DIR:-"carma"}"
  
  PACKDIR="${CCIL:?}/${CNAME:?}"

  FOLDER="carma_tmp"

  URL="https://github.com/RUrlus/carma.git"

  VER=v0.7.0

  XZF="carma.xz"

  # --------------------------------------------------------------------------
  # In case this script runs twice -------------------------------------------
  # --------------------------------------------------------------------------
  if [[ -n "${OVERWRITE_EXISTING_CORE_PACKAGES:-}" ]]; then
    
    rm -rf "${CCIL:?}/${XZF:?}"
  
  fi
  # splitting also catches partial states
  if [[ ! -f "${CCIL:?}/${XZF:?}" ]]; then 
    
    rm -rf "${CCIL:?}/${FOLDER:?}"
    rm -rf "${PACKDIR:?}"
    rm -rf "${CCIL:?}/include"
  
  fi
 
  if [ ! -f "${CCIL:?}/${XZF:?}" ]; then
    
    gitact1 "${FOLDER:?}" "${VER:?}" "${URL:?}" || { unset_all; return 1; }

    # ------------------------------------------------------------------------
    # move/rename include file and carma.h folder
    # ------------------------------------------------------------------------
    [ -d "${CCIL:?}/${FOLDER:?}/include" ] || { error "CARMA include missing"; return 1; }
    [ ! -e "${PACKDIR:?}" ] || { error "destination already exists"; return 1; }

    mv "${CCIL:?}/${FOLDER:?}/include" "${CCIL:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "MV CARMA INCLUDE FOLDER"; return 1; }

    mv "${CCIL:?}/include" "${PACKDIR:?}" \
      2>>${OUT2:?} || { error "RENAME CARMA INCLUDE FOLDER"; return 1; }

    mv "${PACKDIR:?}/carma" "${PACKDIR:?}/carma.h" \
      2>>${OUT2:?} || { error "RENANE CARMA HEADER"; return 1; }

    # --------------------------------------------------------------------------
    
    gitact2 "${CNAME:?}" "${XZF:?}" || { unset_all; return 1; }

  fi

  unset -v CNAME PACKDIR URL FOLDER VER XZF

  pbottom "GETTING CARMA LIBRARY DONE (CORE LIBS)" || { unset_all; return 1; }

fi

# ----------------------------------------------------------------------------

unset_all || return 1;

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------