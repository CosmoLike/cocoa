#!/bin/bash
# ------------------------------------------------------------------------------
# setup_bcemu.sh: download BCemu (baryonic correction model emulators).
#
# Sourced, never executed (only `source` keeps its environment changes).
# setup_cocoa.sh runs it during a full installation, and sourcing it on its
# own redoes just this package (see the README FAQ on installing a single
# package). Setup scripts MAY use the internet (git clone, wget, pip index);
# the matching compile step never does, so anything that must be downloaded
# happens here.
#
# Skipped (return 99) when IGNORE_BCEMU_CODE is set. Reruns are safe: existing
# results are kept unless OVERWRITE_EXISTING_BCEMU_CODE is set. On failure
# error() names the failing step and unsets everything defined here; on
# success the script returns 55, the value the runner caches so a finished
# step is not repeated.
# ------------------------------------------------------------------------------
if [ -n "${IGNORE_BCEMU_CODE:-}" ]; then
  return 99
fi

if [ -z "${ROOTDIR:-}" ]; then
  source start_cocoa.sh || { pfail 'ROOTDIR'; return 1; }
fi

# parenthesis = run in a subshell
( source "${ROOTDIR:?}/installation_scripts/flags_check.sh" ) || return 1;

# unset_env_vars: forget every variable the body defines. Sourced scripts
# share the caller's shell, so anything not unset here would leak into the
# user's environment. cdroot returns to the directory the user launched from.
unset_env_vars () {
  unset -v URL CCIL ECODEF FOLDER PACKDIR PLIB
  cdroot || return 1;
}

# unset_env_funcs: forget every helper function defined below (functions leak
# into the user's shell exactly like variables).
unset_env_funcs () {
  unset -f cdfolder cpfolder error cpfile
  unset -f unset_env_funcs
  cdroot || return 1;
}

# unset_all: the single cleanup entry point: variables, then functions, then
# itself.
unset_all () {
  unset_env_vars
  unset_env_funcs
  unset -f unset_all
  cdroot || return 1;
}

# error: print which step failed (fail_script_msg names this script and the
# message) and run the FULL cleanup. A call site whose helper already routes
# through error() therefore uses plain `|| return 1`; adding another unset_all
# there would run the cleanup twice.
error () {
  fail_script_msg "$(basename "${BASH_SOURCE[0]}")" "${1}"
  unset_all || return 1
}

# cdfolder: cd that reports the target and cleans up on failure.
cdfolder() {
  cd "${1:?}" 2>"/dev/null" || { error "CD FOLDER: ${1}"; return 1; }
}

# cpfolder: recursive copy that reports source/target and cleans up on
# failure.
cpfolder() {
  cp -r "${1:?}" "${2:?}"  \
    2>"/dev/null" || { error "CP FOLDER ${1} on ${2}"; return 1; }
}

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

unset_env_vars || return 1

CCIL="${ROOTDIR:?}/../cocoa_installation_libraries"

# E = EXTERNAL, CODE, F=FODLER
ECODEF="${ROOTDIR:?}/external_modules/code"

URL="${BCEMU_URL:-"https://github.com/sambit-giri/BCemu.git"}"

FOLDER="${BCEMU_NAME:-"bcemu"}"

PACKDIR="${ECODEF:?}/${FOLDER:?}"
ptop "INSTALLING BARYONIC BCEMU FEEDBACK MODELS" || { unset_all; return 1; }

# ----------------------------------------------------------------------------
# In case this script is called twice ----------------------------------------
# ----------------------------------------------------------------------------
if [ -n "${OVERWRITE_EXISTING_BCEMU_CODE:-}" ]; then
  rm -rf "${PACKDIR:?}"
fi

if [ ! -d "${PACKDIR:?}" ]; then
  echo "${PACKDIR:?}"

  cdfolder "${ECODEF:?}" || return 1;

  "${GIT:?}" clone "${URL:?}" --depth ${GIT_CLONE_MAXIMUM_DEPTH:-1000} \
    --recursive --no-single-branch "${PACKDIR:?}" \
    >>${OUT1:?} 2>>${OUT2:?} || { error "${EC15:?}"; return 1; }

  cdfolder "${PACKDIR:?}" || return 1;

  if [[ -n "${BCEMU_GIT_COMMIT:-}" ||
        -n "${BCEMU_GIT_BRANCH:-}" ||
        -n "${BCEMU_GIT_TAG:-}" ]]; then
    if [ "$("${GIT:?}" rev-parse --is-shallow-repository)" = "true" ]; then
      "${GIT:?}" fetch --unshallow --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    else
      "${GIT:?}" fetch --all --tags --prune \
        >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }
    fi
  fi

  if [ -n "${BCEMU_GIT_COMMIT:-}" ]; then

    "${GIT:?}" checkout "${BCEMU_GIT_COMMIT:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }

  elif [ -n "${BCEMU_GIT_BRANCH:-}" ]; then

    "${GIT:?}" checkout -b "${BCEMU_GIT_BRANCH:?}" "origin/${BCEMU_GIT_BRANCH:?}" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }

  elif [ -n "${BCEMU_GIT_TAG:-}" ]; then

    "${GIT:?}" checkout "tags/${BCEMU_GIT_TAG:?}" -b "${BCEMU_GIT_TAG:?}TMP" \
      >>${OUT1:?} 2>>${OUT2:?} || { error "${EC16:?}"; return 1; }

  fi

fi

# Cocoa patch: BCemu's __init__ imports its spectra module, which imports camb
# at package import. That hijacks cobaya's path-checked CAMB load (bfmt only
# needs the BCM_* emulator, not BCemu's CAMB helpers). Drop the import.
sed --in-place --regexp-extended \
  's@^from \.spectra import@#from .spectra import@' \
  "${PACKDIR:?}/src/BCemu/__init__.py" \
  2>>${OUT2:?} || { error "SED BCEMU SPECTRA PATCH"; return 1; }

# ------------------------------------------------------------------------------
# Install BCemu on .local and download its trained emulator files. Both steps
# may need internet and therefore live in setup (not compile): compile_*.sh
# scripts must stay offline-safe (they may run on compute nodes)
# ------------------------------------------------------------------------------

PLIB="${ROOTDIR:?}/.local/lib/python${PYTHON_VERSION:?}/site-packages"
rm -rf "${PLIB:?}/BCemu"
rm -rf "${PLIB:?}"/BCemu-*
rm -rf "${PLIB:?}"/bcemu-*

(
  env CXX="${CXX_COMPILER:?}" CC="${C_COMPILER:?}" \
    ${PIP3:?} install "${PACKDIR:?}" \
      --no-dependencies \
      --prefix="${ROOTDIR:?}/.local" \
      --no-index \
      --no-build-isolation
) >>${OUT1:?} 2>>${OUT2:?} || { error "${EC3:?}"; return 1; }

# Existing emulator files are skipped; without this step, BCemu would download
# them at the FIRST MCMC evaluation, which fails on offline compute nodes
(
  env PYTHONPATH="${PLIB:?}:${PYTHONPATH:-}" ${PYTHON3:?} -c \
    "from BCemu.download import download_emulators; \
     download_emulators(model_name='BCemu2021'); \
     download_emulators(model_name='BCemu2025')"
) >>${OUT1:?} 2>>${OUT2:?} || { error "BCEMU EMULATOR FILES DOWNLOAD"; return 1; }

cdfolder "${ROOTDIR:?}" || return 1;

pbottom "INSTALLING BARYONIC BCEMU FEEDBACK MODELS" || { unset_all; return 1; }

# ---------------------------------------------------------------------------

unset_all || return 1

#-------------------------------------------------------------------------------

return 55; # why this odd number? Setup_cocoa will cache this installation only
           #   if this script runs entirely. What if the user close the terminal
           #   or the system shuts down in the middle of a git clone?
           #   In this case, PACKDIR would exists, but it is corrupted

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
