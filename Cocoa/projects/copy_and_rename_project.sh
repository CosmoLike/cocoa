#!/bin/bash
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

error_cem_msg () {
  local FILE="$(basename ${BASH_SOURCE[0]})"
  local MSG="\033[0;31m\t\t (${FILE}) we cannot run "
  local MSG2="\033[0m"
  echo -e "${MSG}${1:?}${MSG2}" 2>"/dev/null"
}

error_cem () {
  error_cem_msg ${1:?}
  unset -v SCRIPTS ERRORCODE
  unset -f error_cem error_cem_msg
  cd $(pwd -P) 2>"/dev/null"
  source stop_cocoa 2>"/dev/null"
  return 1
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

OLD_PROJECT="lsst_y1"
OLD_SURVEY="LSST"

NEW_PROJECT="roman_real_y1"
NEW_SURVEY="ROMAN" 

PRJ="${ROOTDIR:?}/projects/${NEW_PROJECT}"

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

cd "${PRJ:?}/data/"
find . -iname "*${OLD_PROJECT}*" -exec rename ${OLD_PROJECT} ${NEW_PROJECT,,} '{}' \;

cd "${PRJ:?}/likelihood/"
find . -iname "*${OLD_PROJECT}*" -exec rename ${OLD_PROJECT} ${NEW_PROJECT,,} '{}' \;

cd "${PRJ:?}/scripts/"
find . -iname "*${OLD_PROJECT}*" -exec rename ${OLD_PROJECT} ${NEW_PROJECT,,} '{}' \;

cd "${PRJ:?}/interface/"
find . -iname "*${OLD_PROJECT}*" -exec rename ${OLD_PROJECT} ${NEW_PROJECT,,} '{}' \;
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
  sed --in-place --regexp-extended "s@${OLD_PROJECT}@${NEW_PROJECT,,}@g" \
    "${PRJ:?}/${TMP[$i]}" 2>/dev/null
  sed --in-place --regexp-extended "s@${OLD_PROJECT^^}@${NEW_PROJECT,,}@g" \
    "${PRJ:?}/${TMP[$i]}" 2>/dev/null
  sed --in-place --regexp-extended "s@${OLD_PROJECT,,}@${NEW_PROJECT,,}@g" \
    "${PRJ:?}/${TMP[$i]}" 2>/dev/null

  sed --in-place --regexp-extended  "s@"${OLD_SURVEY}"@"${NEW_SURVEY,,}"@g" \
    "${PRJ:?}/${TMP[$i]}" 2>/dev/null
  sed --in-place --regexp-extended  "s@"${OLD_SURVEY^^}"@"${NEW_SURVEY,,}"@g" \
    "${PRJ:?}/${TMP[$i]}" 2>/dev/null
  sed --in-place --regexp-extended  "s@"${OLD_SURVEY,,}"@"${NEW_SURVEY,,}"@g" \
    "${PRJ:?}/${TMP[$i]}" 2>/dev/null
done


for f in ${PRJ}/{,likelihood/,interface/,data/,scripts/}*.{sh,py,cpp,dataset,yaml}; do
  [ -e "$f" ] || continue

  sed --in-place --regexp-extended "s@${OLD_PROJECT}@${NEW_PROJECT,,}@g" "${f}" 2>/dev/null
  sed --in-place --regexp-extended "s@${OLD_PROJECT^^}@${NEW_PROJECT,,}@g" "${f}" 2>/dev/null
  sed --in-place --regexp-extended "s@${OLD_PROJECT,,}@${NEW_PROJECT,,}@g" "${f}" 2>/dev/null

  sed --in-place --regexp-extended  "s@"${OLD_SURVEY}"@"${NEW_SURVEY,,}"@g" "${f}" 2>/dev/null
  sed --in-place --regexp-extended  "s@"${OLD_SURVEY^^}"@"${NEW_SURVEY,,}"@g" "${f}" 2>/dev/null
  sed --in-place --regexp-extended  "s@"${OLD_SURVEY,,}"@"${NEW_SURVEY,,}"@g" "${f}" 2>/dev/null
done

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

rm -f "${PRJ:?}"/*.txt                  2>/dev/null
rm -f "${PRJ:?}"/*.sbatch               2>/dev/null
rm -f "${PRJ:?}"/*.ipynb                2>/dev/null
rm -f "${PRJ:?}"/interface/*.so         2>/dev/null
rm -f "${PRJ:?}"/interface/*.o          2>/dev/null
rm -f "${PRJ:?}"/chains/*.txt           2>/dev/null
rm -f "${PRJ:?}"/chains/*.1.txt         2>/dev/null
rm -f "${PRJ:?}"/chains/*.2.txt         2>/dev/null
rm -f "${PRJ:?}"/chains/*.3.txt         2>/dev/null
rm -f "${PRJ:?}"/chains/*.4.txt         2>/dev/null
rm -f "${PRJ:?}"/chains/*.5.txt         2>/dev/null
rm -f "${PRJ:?}"/chains/*.6.txt         2>/dev/null
rm -f "${PRJ:?}"/chains/*.7.txt         2>/dev/null
rm -f "${PRJ:?}"/chains/*.8.txt         2>/dev/null
rm -f "${PRJ:?}"/chains/*.progress      2>/dev/null
rm -f "${PRJ:?}"/chains/*.covmat        2>/dev/null
rm -f "${PRJ:?}"/chains/*.locked        2>/dev/null
rm -f "${PRJ:?}"/chains/*.checkpoint    2>/dev/null
rm -f "${PRJ:?}"/chains/*.py.           2>/dev/null
rm -f "${PRJ:?}"/chains/*.yaml.         2>/dev/null
rm -f "${PRJ:?}"/chains/*.input.yaml    2>/dev/null
rm -f "${PRJ:?}"/chains/*.updated.yaml  2>/dev/null
rm -f "${PRJ:?}"/chains/*.pyc           2>/dev/null
rm -rf "${PRJ:?}"/.git/                 2>/dev/null
rm -rf "${PRJ:?}"/scripts/random_scripts_used_by_dev 2>/dev/null

# ------------------------------------------------------------------------------------
