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

NEW_PROJECT="roman_real"
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

declare -a TMP=("scripts/compile_"${OLD_PROJECT}".sh"
                "scripts/start_"${OLD_PROJECT}".sh"
                "scripts/stop_"${OLD_PROJECT}".sh"
                "likelihood/${OLD_SURVEY,,}_3x2pt.py"
                "likelihood/${OLD_SURVEY,,}_3x2pt.yaml"
                "likelihood/${OLD_SURVEY,,}_2x2pt.py"
                "likelihood/${OLD_SURVEY,,}_2x2pt.yaml"
                "likelihood/${OLD_SURVEY,,}_ggl.py"
                "likelihood/${OLD_SURVEY,,}_ggl.yaml"
                "likelihood/${OLD_SURVEY,,}_cosmic_shear.py"
                "likelihood/${OLD_SURVEY,,}_cosmic_shear.yaml"
                "likelihood/${OLD_SURVEY,,}_xi_ggl.py"
                "likelihood/${OLD_SURVEY,,}_xi_ggl.yaml"
                "likelihood/${OLD_SURVEY,,}_clustering.py"
                "likelihood/${OLD_SURVEY,,}_clustering.yaml"
                "likelihood/${OLD_SURVEY,,}_xi_gg.py"
                "likelihood/${OLD_SURVEY,,}_xi_gg.yaml"
                "likelihood/params_${OLD_SURVEY,,}_lens.yaml"
                "likelihood/params_${OLD_SURVEY,,}_source.yaml"
                "data/${OLD_PROJECT}.dataset"
                "interface/cosmolike_${OLD_PROJECT}_interface.py"
                ) 

declare -a TMP2=("scripts/compile_"${NEW_PROJECT,,}".sh"
                 "scripts/start_"${NEW_PROJECT,,}".sh"
                 "scripts/stop_"${NEW_PROJECT,,}".sh"
                 "likelihood/${NEW_SURVEY,,}_3x2pt.py"
                 "likelihood/${NEW_SURVEY,,}_3x2pt.yaml"
                 "likelihood/${NEW_SURVEY,,}_2x2pt.py"
                 "likelihood/${NEW_SURVEY,,}_2x2pt.yaml"
                 "likelihood/${NEW_SURVEY,,}_ggl.py"
                 "likelihood/${NEW_SURVEY,,}_ggl.yaml"
                 "likelihood/${NEW_SURVEY,,}_cosmic_shear.py"
                 "likelihood/${NEW_SURVEY,,}_cosmic_shear.yaml"
                 "likelihood/${NEW_SURVEY,,}_xi_ggl.py"
                 "likelihood/${NEW_SURVEY,,}_xi_ggl.yaml"
                 "likelihood/${NEW_SURVEY,,}_clustering.py"
                 "likelihood/${NEW_SURVEY,,}_clustering.yaml"
                 "likelihood/${NEW_SURVEY,,}_xi_gg.py"
                 "likelihood/${NEW_SURVEY,,}_xi_gg.yaml"
                 "likelihood/params_${NEW_SURVEY,,}_lens.yaml"
                 "likelihood/params_${NEW_SURVEY,,}_source.yaml"
                 "data/${NEW_PROJECT,,}.dataset"
                 "interface/cosmolike_${NEW_PROJECT,,}_interface.py"
                )

for (( i=0; i<${#TMP[@]}; i++ ));
do
  mv "${PRJ:?}/${TMP[$i]:?}" "${PRJ:?}/${TMP2[$i]:?}"     2>/dev/null

  sed --in-place --regexp-extended "s@${OLD_PROJECT}@${NEW_PROJECT,,}@g" \
    "${PRJ:?}/${TMP2[$i]}" 2>/dev/null
  sed --in-place --regexp-extended "s@${OLD_PROJECT^^}@${NEW_PROJECT,,}@g" \
    "${PRJ:?}/${TMP2[$i]}" 2>/dev/null
  sed --in-place --regexp-extended "s@${OLD_PROJECT,,}@${NEW_PROJECT,,}@g" \
    "${PRJ:?}/${TMP2[$i]}" 2>/dev/null

  sed --in-place --regexp-extended  "s@"${OLD_SURVEY}"@"${NEW_SURVEY,,}"@g" \
    "${PRJ:?}/${TMP2[$i]}" 2>/dev/null
  sed --in-place --regexp-extended  "s@"${OLD_SURVEY^^}"@"${NEW_SURVEY,,}"@g" \
    "${PRJ:?}/${TMP2[$i]}" 2>/dev/null
  sed --in-place --regexp-extended  "s@"${OLD_SURVEY,,}"@"${NEW_SURVEY,,}"@g" \
    "${PRJ:?}/${TMP2[$i]}" 2>/dev/null
done

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
declare -a TMP=("EXAMPLE_EVALUATE1.yaml"
                "EXAMPLE_EVALUATE2.yaml"
                "EXAMPLE_EVALUATE3.yaml"
                "EXAMPLE_EVALUATE4.yaml"
                "EXAMPLE_EVALUATE5.yaml"
                "EXAMPLE_MCMC1.yaml"
                "EXAMPLE_MCMC2.yaml"
                "EXAMPLE_MCMC3.yaml"
                "EXAMPLE_MCMC4.yaml"
                "EXAMPLE_PROFILE1.yaml"
                "interface/interface.cpp"
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

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

sed --in-place --regexp-extended \
  "s@cosmolike_"${OLD_PROJECT:?}"_interface@cosmolike_"${NEW_PROJECT}"_interface@g" \
  "${PRJ:?}/likelihood/_cosmolike_prototype_base.py" 2>/dev/null

sed --in-place --regexp-extended "s@${OLD_SURVEY:?}@${NEW_SURVEY:?}@g" \
  "${PRJ:?}/likelihood/_cosmolike_prototype_base.py" 2>/dev/null

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
rm -f "${PRJ:?}"/chains/*.py.           2>/dev/null
rm -f "${PRJ:?}"/chains/*.yaml.         2>/dev/null
rm -f "${PRJ:?}"/chains/*.input.yaml    2>/dev/null
rm -f "${PRJ:?}"/chains/*.updated.yaml  2>/dev/null
rm -f "${PRJ:?}"/chains/*.pyc           2>/dev/null
rm -rf "${PRJ:?}"/.git/                 2>/dev/null
rm -rf "${PRJ:?}"/scripts/random_scripts_used_by_dev 2>/dev/null

# ------------------------------------------------------------------------------------
