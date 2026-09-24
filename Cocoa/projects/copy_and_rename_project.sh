#!/bin/bash
# ------------------------------------------------------------------------------
# copy_and_rename_project.sh: create a new Cosmolike project from a donor
# ------------------------------------------------------------------------------
# Copies projects/OLD_PROJECT into projects/NEW_PROJECT and renames the
# project and survey names in file names and file contents: code, yamls,
# datasets, the parameter names on the covmat header lines, the Git LFS
# patterns in .gitattributes, and the unit-test suite. It then deletes
# what does not transfer to a new survey (chains, caches, the data-vector
# emulators and the EXAMPLE_EMUL_* examples that load them, the tests'
# frozen references) and replaces the donor README with a stub.
#
# Usage (from the Cocoa/ folder, cocoa environment active):
#
#     1. edit NEW_PROJECT / NEW_SURVEY below
#     2. bash ./projects/copy_and_rename_project.sh
#     3. source start_cocoa.sh   # recreates the project symlinks
#     4. source ./projects/<new>/scripts/compile_<new>.sh
#
# Runs on Linux and macOS: bash 3.2 suffices, and GNU sed comes from the
# cocoa environment (the script aborts when GNU sed is missing). The
# script only renames: data vectors, covariances, masks, and n(z) stay
# the donor's, and the tests' frozen references must be regenerated with
# tests/generate_frozen_reference.py --overwrite. The FAQ "How do we
# create a new Cosmolike project?" in projects/README.md documents the
# remaining manual steps.
# ------------------------------------------------------------------------------

OLD_PROJECT="lsst_y1"
OLD_SURVEY="LSST"

NEW_PROJECT="xxx"
NEW_SURVEY="XXX"

# Case variants precomputed with tr so the script runs on bash 3.2
# (macOS /bin/bash): the ${var,,} and ${var^^} expansions need bash 4.
OLD_PROJECT_U=$(echo "${OLD_PROJECT:?}" | tr '[:lower:]' '[:upper:]')
OLD_PROJECT_L=$(echo "${OLD_PROJECT:?}" | tr '[:upper:]' '[:lower:]')
NEW_PROJECT_L=$(echo "${NEW_PROJECT:?}" | tr '[:upper:]' '[:lower:]')
OLD_SURVEY_U=$(echo "${OLD_SURVEY:?}" | tr '[:lower:]' '[:upper:]')
OLD_SURVEY_L=$(echo "${OLD_SURVEY:?}" | tr '[:upper:]' '[:lower:]')
NEW_SURVEY_U=$(echo "${NEW_SURVEY:?}" | tr '[:lower:]' '[:upper:]')
NEW_SURVEY_L=$(echo "${NEW_SURVEY:?}" | tr '[:upper:]' '[:lower:]')

# GNU sed ships with the cocoa environment on Linux and macOS; BSD sed
# (the macOS default) lacks --in-place --regexp-extended and would fail
# silently below.
if ! sed --version >/dev/null 2>&1; then
  echo "GNU sed not found: activate the cocoa environment first" \
       "(conda activate cocoa; source start_cocoa.sh)" >&2
  return 1 2>/dev/null || exit 1
fi

PRJ="${ROOTDIR:?}/projects/${NEW_PROJECT:?}"

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

# Rename every file and folder under ${1} whose name contains ${2},
# replacing ${2} with ${3}. Line-by-line notes below, since this is
# advanced bash: ${1} ${2} ${3} are the function's positional
# arguments, and ${1:?} means "expand ${1}, but abort with an error
# when it is empty or unset" (the same guard every rm in this script
# uses).
rename_tree () {
  # Work from inside the target folder, so find prints short relative
  # paths (./some_file.py). If cd fails the function must stop at
  # once: find would otherwise run in whatever folder the shell
  # happened to be in and rename files THERE. "cmd || action" runs
  # action only when cmd fails, and >&2 sends the message to stderr.
  cd "${1:?}" || { echo "missing folder: ${1}" >&2; return 1; }

  # find lists everything (files AND folders) whose name contains ${2}:
  #   -iname "*${2}*"  matches the name case-insensitively;
  #   -depth           lists a folder's CONTENTS before the folder
  #                    itself, so renaming the folder never invalidates
  #                    paths still waiting on the list;
  #   -print0          separates results with a \0 byte instead of a
  #                    newline. \0 is the one character a file name can
  #                    never contain, so names with spaces or newlines
  #                    survive the pipe to read intact.
  # Example: with ${2} = lsst_y1, run inside the copied data folder,
  # the output is one stream of \0-separated relative names:
  #     ./lsst_y1_cov\0./lsst_y1_source.nz\0./lsst_y1_M1_GGL0.05.dataset\0...
  find . -depth -iname "*${2:?}*" -print0 |
    # read pulls the \0-separated names off the pipe one at a time:
    #   IFS=   empty field separator, so surrounding spaces are kept;
    #   -r     raw mode: backslashes in names are not special;
    #   -d ''  read up to the next \0, matching -print0 above.
    while IFS= read -r -d '' f; do
      # ${f//old/new} is bash pattern substitution: ${f} with EVERY
      # occurrence of ${2} replaced by ${3} (a single slash,
      # ${f/old/new}, would replace only the first occurrence).
      # Example: f = ./lsst_y1_cov, ${2} = lsst_y1, ${3} = xxx
      #          ->  g = ./xxx_cov
      g="${f//${2}/${3:?}}"

      # -iname matched case-insensitively but the substitution above
      # is case-sensitive, so g can come out equal to f (a file named
      # LSST_Y1_notes when ${2} is lsst_y1). "mv name name" would be
      # an error, so rename only when the name actually changed:
      # "[ a = b ] || cmd" runs cmd only when the test fails.
      # Example: f = ./lsst_y1_cov    g = ./xxx_cov        -> mv runs
      #          f = ./LSST_Y1_notes  g = ./LSST_Y1_notes  -> skipped
      [ "$f" = "$g" ] || mv "$f" "$g"
    done
}

# Replace every case variant of the project and survey names inside the
# file ${1}. Order matters: the project expressions run first, so an
# uppercase project name (LSST_Y1) is consumed as a project name before
# the survey substitution (LSST) can touch it.
rename_in_file () {
  sed --in-place --regexp-extended \
    -e "s@${OLD_PROJECT}@${NEW_PROJECT_L}@g" \
    -e "s@${OLD_PROJECT_U}@${NEW_PROJECT_L}@g" \
    -e "s@${OLD_PROJECT_L}@${NEW_PROJECT_L}@g" \
    -e "s@${OLD_SURVEY}@${NEW_SURVEY_U}@g" \
    -e "s@${OLD_SURVEY_U}@${NEW_SURVEY_U}@g" \
    -e "s@${OLD_SURVEY_L}@${NEW_SURVEY_U}@g" \
    "${1:?}" 2>/dev/null
}

# ------------------------------------------------------------------------------
# Copy the donor project
# ------------------------------------------------------------------------------

rm -rf "${PRJ:?}"
cp -r "${ROOTDIR:?}/projects/${OLD_PROJECT:?}" "${PRJ:?}"

# ------------------------------------------------------------------------------
# Rename file and folder NAMES
# ------------------------------------------------------------------------------

declare -a MVDIRS=("data"
                   "likelihood"
                   "scripts"
                   "interface"
                  )

for (( i=0; i<${#MVDIRS[@]}; i++ ));
do
  rename_tree "${PRJ:?}/${MVDIRS[$i]}" "${OLD_PROJECT}" "${NEW_PROJECT_L}" || \
    { return 1 2>/dev/null || exit 1; }
done

rename_tree "${PRJ:?}/interface" "${OLD_SURVEY}" "${NEW_SURVEY_L}" || \
  { return 1 2>/dev/null || exit 1; }

# ------------------------------------------------------------------------------
# Rename file CONTENTS
# ------------------------------------------------------------------------------

# Named files the extension globs below cannot reach: the Makefile, and
# .gitattributes (its Git LFS pattern names the covariance data file; a
# stale pattern would let the renamed large file escape LFS tracking).
declare -a SEDFILES=("interface/MakefileCosmolike"
                     ".gitattributes"
                    )

for (( i=0; i<${#SEDFILES[@]}; i++ ));
do
  rename_in_file "${PRJ:?}/${SEDFILES[$i]}"
done

# covmat included: the first line of a proposal covariance names the
# sampled parameters, which carry the old survey prefix.
declare -a SEDDIRS=(""
                    "likelihood/"
                    "interface/"
                    "data/"
                    "scripts/"
                   )

for (( i=0; i<${#SEDDIRS[@]}; i++ ));
do
  for f in "${PRJ:?}/${SEDDIRS[$i]}"*.{sh,py,cpp,dataset,yaml,covmat}; do
    [ -e "$f" ] || continue
    rename_in_file "$f"
  done
done

# tests/ ships the project unit-test suite: the code and README carry the
# project name in imports, likelihood references, and parameter prefixes.
for f in "${PRJ:?}"/tests/*.{py,md}; do
  [ -e "$f" ] || continue
  rename_in_file "$f"
done

# ------------------------------------------------------------------------------
# Delete what does not transfer to a new survey
# ------------------------------------------------------------------------------

# The data-vector emulators (and the EXAMPLE_EMUL_* examples that load
# them) were trained on the old survey, so they are untransferable; the
# hybrid EXAMPLE_EMUL2_* examples emulate only Boltzmann outputs and are
# kept. The tests' frozen state (frozen/ + manifest_sha256.json) is a
# SHA-256 pinned snapshot of the old survey's data, and the figures are
# measurements of it: none transfer. Once the new survey's data files
# are in place, regenerate them with
# tests/generate_frozen_reference.py --overwrite (see tests/README.md).
declare -a CLEAN_FILES=("*.txt"
                        "*.sbatch"
                        "*.ipynb"
                        "EXAMPLE_EMUL_*.yaml"
                        "EXAMPLE_EMUL_*.py"
                        "interface/*.so"
                        "interface/*.o"
                        "scripts/EXAMPLE_PLOT_*.py" # they plot emulator chains
                        "scripts/*.sbatch"
                        "chains/*.txt"          # also covers rank-suffixed *.N.txt chains
                        "chains/*.progress"
                        "chains/*.covmat"
                        "chains/*.locked"
                        "chains/*.checkpoint"
                        "chains/*.py."
                        "chains/*.yaml."
                        "chains/*.input.yaml"
                        "chains/*.updated.yaml"
                        "chains/*.pyc"
                        "tests/manifest_sha256.json"
                        "tests/*.png"
                       )

for (( i=0; i<${#CLEAN_FILES[@]}; i++ ));
do
  rm -f "${PRJ:?}"/${CLEAN_FILES[$i]} 2>/dev/null
done

declare -a CLEAN_FOLDERS=("emulators"
                          ".git"
                          "interface/__pycache__"
                          "likelihood/__pycache__"
                          "tests/__pycache__"
                          "tests/frozen"
                          "scripts/random_scripts_used_by_dev"
                         )

for (( i=0; i<${#CLEAN_FOLDERS[@]}; i++ ));
do
  rm -rf "${PRJ:?}/${CLEAN_FOLDERS[$i]}" 2>/dev/null
done

# ------------------------------------------------------------------------------
# Replace the donor README with a stub
# ------------------------------------------------------------------------------

# The top-level README documents the OLD survey (its releases, pinned
# installation keys, and data provenance), so renaming it would only
# fabricate a history the new project never had.
cat > "${PRJ:?}/README.md" <<EOF
# The ${NEW_PROJECT_L} project

Created by projects/copy_and_rename_project.sh from a donor project.
The donor README does not transfer, so this stub replaces it: write
here the new survey documentation (data provenance, scale cuts,
releases). The FAQ "How do we create a new Cosmolike project?" in
projects/README.md lists the remaining manual steps, and
tests/README.md documents the unit-test suite.
EOF

# ------------------------------------------------------------------------------
# Leave no variables or functions behind (the script may be sourced)
# ------------------------------------------------------------------------------

unset -v PRJ OLD_PROJECT OLD_SURVEY NEW_PROJECT NEW_SURVEY
unset -v OLD_PROJECT_U OLD_PROJECT_L NEW_PROJECT_L
unset -v OLD_SURVEY_U OLD_SURVEY_L NEW_SURVEY_U NEW_SURVEY_L
unset -v MVDIRS SEDFILES SEDDIRS CLEAN_FILES CLEAN_FOLDERS i f
unset -f rename_tree rename_in_file

# ------------------------------------------------------------------------------
