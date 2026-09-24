#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
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
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

rm -rf "${PRJ:?}"
cp -r "${ROOTDIR:?}/projects/${OLD_PROJECT:?}" "${PRJ:?}"

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# mv + bash pattern substitution replaces the linux-only rename tool
# (validated byte-identical to it on the full lsst_y1 tree; the rename
# tool is also flavor-dependent: the perl variant shipped by
# Debian-family systems would misread these arguments). -depth lists
# the contents of a folder before the folder itself, so a renamed
# folder never invalidates the paths still on the list.
for d in data likelihood scripts interface; do
  cd "${PRJ:?}/${d}/"
  find . -depth -iname "*${OLD_PROJECT}*" -print0 | \
    while IFS= read -r -d '' f; do
      g="${f//${OLD_PROJECT}/${NEW_PROJECT_L}}"
      [ "$f" = "$g" ] || mv "$f" "$g"
    done
done

cd "${PRJ:?}/interface/"
find . -depth -iname "*${OLD_SURVEY}*" -print0 | \
  while IFS= read -r -d '' f; do
    g="${f//${OLD_SURVEY}/${NEW_SURVEY_L}}"
    [ "$f" = "$g" ] || mv "$f" "$g"
  done

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
declare -a TMP=(
                "interface/MakefileCosmolike"
               ) 

for (( i=0; i<${#TMP[@]}; i++ ));
do
  sed --in-place --regexp-extended "s@${OLD_PROJECT}@${NEW_PROJECT_L}@g"   "${PRJ:?}/${TMP[$i]}" 2>/dev/null
  sed --in-place --regexp-extended "s@${OLD_PROJECT_U}@${NEW_PROJECT_L}@g" "${PRJ:?}/${TMP[$i]}" 2>/dev/null
  sed --in-place --regexp-extended "s@${OLD_PROJECT_L}@${NEW_PROJECT_L}@g" "${PRJ:?}/${TMP[$i]}" 2>/dev/null

  sed --in-place --regexp-extended "s@${OLD_SURVEY}@${NEW_SURVEY_U}@g"   "${PRJ:?}/${TMP[$i]}" 2>/dev/null
  sed --in-place --regexp-extended "s@${OLD_SURVEY_U}@${NEW_SURVEY_U}@g" "${PRJ:?}/${TMP[$i]}" 2>/dev/null
  sed --in-place --regexp-extended "s@${OLD_SURVEY_L}@${NEW_SURVEY_U}@g" "${PRJ:?}/${TMP[$i]}" 2>/dev/null
done


for f in ${PRJ}/{,likelihood/,interface/,data/,scripts/}*.{sh,py,cpp,dataset,yaml}; do
  [ -e "$f" ] || continue

  sed --in-place --regexp-extended "s@${OLD_PROJECT}@${NEW_PROJECT_L}@g" "${f}" 2>/dev/null
  sed --in-place --regexp-extended "s@${OLD_PROJECT_U}@${NEW_PROJECT_L}@g" "${f}" 2>/dev/null
  sed --in-place --regexp-extended "s@${OLD_PROJECT_L}@${NEW_PROJECT_L}@g" "${f}" 2>/dev/null

  sed --in-place --regexp-extended "s@${OLD_SURVEY}@${NEW_SURVEY_U}@g" "${f}" 2>/dev/null
  sed --in-place --regexp-extended "s@${OLD_SURVEY_U}@${NEW_SURVEY_U}@g" "${f}" 2>/dev/null
  sed --in-place --regexp-extended "s@${OLD_SURVEY_L}@${NEW_SURVEY_U}@g" "${f}" 2>/dev/null
done

# tests/ ships the project unit-test suite: the code and README carry the
# project name in imports, likelihood references, and parameter prefixes.
for f in ${PRJ}/tests/*.{py,md}; do
  [ -e "$f" ] || continue

  sed --in-place --regexp-extended "s@${OLD_PROJECT}@${NEW_PROJECT_L}@g" "${f}" 2>/dev/null
  sed --in-place --regexp-extended "s@${OLD_PROJECT_U}@${NEW_PROJECT_L}@g" "${f}" 2>/dev/null
  sed --in-place --regexp-extended "s@${OLD_PROJECT_L}@${NEW_PROJECT_L}@g" "${f}" 2>/dev/null

  sed --in-place --regexp-extended "s@${OLD_SURVEY}@${NEW_SURVEY_U}@g" "${f}" 2>/dev/null
  sed --in-place --regexp-extended "s@${OLD_SURVEY_U}@${NEW_SURVEY_U}@g" "${f}" 2>/dev/null
  sed --in-place --regexp-extended "s@${OLD_SURVEY_L}@${NEW_SURVEY_U}@g" "${f}" 2>/dev/null
done

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

rm -f "${PRJ:?}"/*.txt                   2>/dev/null
rm -f "${PRJ:?}"/*.sbatch                2>/dev/null
rm -f "${PRJ:?}"/*.ipynb                 2>/dev/null
# The data-vector emulators (and the EXAMPLE_EMUL_* examples that load them)
# were trained on the old survey, so they are untransferable. The hybrid
# EXAMPLE_EMUL2_* examples emulate only Boltzmann outputs and are kept.
rm -rf "${PRJ:?}"/emulators              2>/dev/null
rm -f "${PRJ:?}"/EXAMPLE_EMUL_*.yaml     2>/dev/null
rm -f "${PRJ:?}"/EXAMPLE_EMUL_*.py       2>/dev/null
rm -f "${PRJ:?}"/interface/*.so          2>/dev/null
rm -f "${PRJ:?}"/interface/*.o           2>/dev/null
rm -f "${PRJ:?}"/chains/*.txt            2>/dev/null # also covers rank-suffixed *.N.txt chains
rm -f "${PRJ:?}"/chains/*.progress       2>/dev/null
rm -f "${PRJ:?}"/chains/*.covmat         2>/dev/null
rm -f "${PRJ:?}"/chains/*.locked         2>/dev/null
rm -f "${PRJ:?}"/chains/*.checkpoint     2>/dev/null
rm -f "${PRJ:?}"/chains/*.py.            2>/dev/null
rm -f "${PRJ:?}"/chains/*.yaml.          2>/dev/null
rm -f "${PRJ:?}"/chains/*.input.yaml     2>/dev/null
rm -f "${PRJ:?}"/chains/*.updated.yaml   2>/dev/null
rm -f "${PRJ:?}"/chains/*.pyc            2>/dev/null
rm -rf "${PRJ:?}"/.git/                  2>/dev/null
rm -rf "${PRJ:?}"/interface/__pycache__  2>/dev/null
rm -rf "${PRJ:?}"/likelihood/__pycache__ 2>/dev/null
rm -rf "${PRJ:?}"/tests/__pycache__      2>/dev/null
rm -rf "${PRJ:?}"/scripts/random_scripts_used_by_dev 2>/dev/null
# The tests' frozen state (frozen/ + manifest_sha256.json) is a SHA-256
# pinned snapshot of the old survey's data, and the figures are measurements
# of it: none transfer. Once the new survey's data files are in place,
# regenerate them with tests/generate_frozen_reference.py --overwrite
# (see tests/README.md).
rm -rf "${PRJ:?}"/tests/frozen               2>/dev/null
rm -f  "${PRJ:?}"/tests/manifest_sha256.json 2>/dev/null
rm -f  "${PRJ:?}"/tests/*.png                2>/dev/null

unset -v PRJ OLD_PROJECT OLD_SURVEY NEW_PROJECT NEW_SURVEY
unset -v OLD_PROJECT_U OLD_PROJECT_L NEW_PROJECT_L
unset -v OLD_SURVEY_U OLD_SURVEY_L NEW_SURVEY_U NEW_SURVEY_L

# ------------------------------------------------------------------------------------
