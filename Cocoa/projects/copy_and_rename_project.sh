#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

OLD_PROJECT="lsst_y1"
OLD_SURVEY="LSST"

NEW_PROJECT="xxx"
NEW_SURVEY="XXX" 

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

for d in data likelihood scripts interface; do
  cd "${PRJ:?}/${d}/"
  find . -iname "*${OLD_PROJECT}*" -exec rename ${OLD_PROJECT} ${NEW_PROJECT,,} '{}' \;
done

cd "${PRJ:?}/interface/"
find . -iname "*${OLD_SURVEY}*" -exec rename ${OLD_SURVEY} ${NEW_SURVEY,,} '{}' \;
find . -iname "*${OLD_SURVEY,,}*" -exec rename ${OLD_SURVEY} ${NEW_SURVEY,,} '{}' \;

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
declare -a TMP=(
                "interface/MakefileCosmolike"
               ) 

for (( i=0; i<${#TMP[@]}; i++ ));
do
  sed --in-place --regexp-extended "s@${OLD_PROJECT}@${NEW_PROJECT,,}@g"   "${PRJ:?}/${TMP[$i]}" 2>/dev/null
  sed --in-place --regexp-extended "s@${OLD_PROJECT^^}@${NEW_PROJECT,,}@g" "${PRJ:?}/${TMP[$i]}" 2>/dev/null
  sed --in-place --regexp-extended "s@${OLD_PROJECT,,}@${NEW_PROJECT,,}@g" "${PRJ:?}/${TMP[$i]}" 2>/dev/null

  sed --in-place --regexp-extended  "s@"${OLD_SURVEY}"@"${NEW_SURVEY^^}"@g"   "${PRJ:?}/${TMP[$i]}" 2>/dev/null
  sed --in-place --regexp-extended  "s@"${OLD_SURVEY^^}"@"${NEW_SURVEY^^}"@g" "${PRJ:?}/${TMP[$i]}" 2>/dev/null
  sed --in-place --regexp-extended  "s@"${OLD_SURVEY,,}"@"${NEW_SURVEY^^}"@g" "${PRJ:?}/${TMP[$i]}" 2>/dev/null
done


for f in ${PRJ}/{,likelihood/,interface/,data/,scripts/}*.{sh,py,cpp,dataset,yaml}; do
  [ -e "$f" ] || continue

  sed --in-place --regexp-extended "s@${OLD_PROJECT}@${NEW_PROJECT,,}@g" "${f}" 2>/dev/null
  sed --in-place --regexp-extended "s@${OLD_PROJECT^^}@${NEW_PROJECT,,}@g" "${f}" 2>/dev/null
  sed --in-place --regexp-extended "s@${OLD_PROJECT,,}@${NEW_PROJECT,,}@g" "${f}" 2>/dev/null

  sed --in-place --regexp-extended  "s@"${OLD_SURVEY}"@"${NEW_SURVEY^^}"@g" "${f}" 2>/dev/null
  sed --in-place --regexp-extended  "s@"${OLD_SURVEY^^}"@"${NEW_SURVEY^^}"@g" "${f}" 2>/dev/null
  sed --in-place --regexp-extended  "s@"${OLD_SURVEY,,}"@"${NEW_SURVEY^^}"@g" "${f}" 2>/dev/null
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
rm -rf "${PRJ:?}"/scripts/random_scripts_used_by_dev 2>/dev/null

unset -v PRJ OLD_PROJECT OLD_SURVEY NEW_PROJECT NEW_SURVEY

# ------------------------------------------------------------------------------------
