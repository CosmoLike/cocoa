#!/bin/bash

OLD_PROJECT="des_y3"
OLD_SURVEY="DES"

NEW_PROJECT="XXX"
NEW_SURVEY="LSST" 

if [ -d $ROOTDIR/projects/$NEW_PROJECT ]; then
  >&2 echo "error: project $NEW_PROJECT already exists at $ROOTDIR/projects/$NEW_PROJECT"
  exit
fi

cp -r $ROOTDIR/projects/$OLD_PROJECT $ROOTDIR/projects/$NEW_PROJECT

# ------------------------------------------------------------------------------------

rm -rf $ROOTDIR/projects/$NEW_PROJECT/.git/ 2>/dev/null

# ------------------------------------------------------------------------------------

sed --in-place --regexp-extended  \
"s@cosmolike_"$OLD_PROJECT"_interface@cosmolike_"$NEW_PROJECT"_interface@g" \
$ROOTDIR/projects/$NEW_PROJECT/interface/MakefileCosmolike 2>/dev/null

# ------------------------------------------------------------------------------------

TMP="cosmolike_"$OLD_PROJECT"_interface.py"
TMP2="cosmolike_"$NEW_PROJECT"_interface.py"
mv $ROOTDIR/projects/$NEW_PROJECT/interface/$TMP \
$ROOTDIR/projects/$NEW_PROJECT/interface/$TMP2 2>/dev/null

sed --in-place --regexp-extended  \
"s@cosmolike_"$OLD_PROJECT"_interface@cosmolike_"$NEW_PROJECT"_interface@g" \
$ROOTDIR/projects/$NEW_PROJECT/interface/$TMP2 2>/dev/null

# ------------------------------------------------------------------------------------

TMP="cosmolike_"$OLD_PROJECT"_interface.so"
rm -rf $ROOTDIR/projects/$NEW_PROJECT/interface/$TMP 2>/dev/null

# ------------------------------------------------------------------------------------

sed --in-place --regexp-extended  \
"s@cosmolike_"$OLD_PROJECT"_interface@cosmolike_"$NEW_PROJECT"_interface@g" \
$ROOTDIR/projects/$NEW_PROJECT/interface/interface.cpp 2>/dev/null


# ------------------------------------------------------------------------------------

TMP="compile_"$OLD_PROJECT
TMP2="compile_"$NEW_PROJECT
mv $ROOTDIR/projects/$NEW_PROJECT/scripts/$TMP \
$ROOTDIR/projects/$NEW_PROJECT/scripts/$TMP2 2>/dev/null

sed --in-place --regexp-extended "s@$OLD_PROJECT@$NEW_PROJECT@g" \
$ROOTDIR/projects/$NEW_PROJECT/scripts/$TMP2 2>/dev/null

# ------------------------------------------------------------------------------------

TMP="start_"$OLD_PROJECT
TMP2="start_"$NEW_PROJECT
mv $ROOTDIR/projects/$NEW_PROJECT/scripts/$TMP \
$ROOTDIR/projects/$NEW_PROJECT/scripts/$TMP2 2>/dev/null

sed --in-place --regexp-extended "s@$OLD_PROJECT@$NEW_PROJECT@g" \
$ROOTDIR/projects/$NEW_PROJECT/scripts/$TMP2 2>/dev/null

# ------------------------------------------------------------------------------------

TMP="stop_"$OLD_PROJECT
TMP2="stop_"$NEW_PROJECT
mv $ROOTDIR/projects/$NEW_PROJECT/scripts/$TMP \
$ROOTDIR/projects/$NEW_PROJECT/scripts/$TMP2 2>/dev/null

# ------------------------------------------------------------------------------------

sed --in-place --regexp-extended \
"s@cosmolike_"$OLD_PROJECT"_interface@cosmolike_"$NEW_PROJECT"_interface@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/_cosmolike_prototype_base.py 2>/dev/null

sed --in-place --regexp-extended "s@$OLD_SURVEY@$NEW_SURVEY@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/_cosmolike_prototype_base.py 2>/dev/null

# ------------------------------------------------------------------------------------

TMP=$OLD_SURVEY"_3x2pt"
TMP2=$NEW_SURVEY"_3x2pt"

mv $ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP,,}".py" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".py" 2>/dev/null

sed --in-place --regexp-extended  "s@$OLD_PROJECT@"$NEW_PROJECT"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".py" 2>/dev/null

sed --in-place --regexp-extended  "s@"${OLD_SURVEY,,}"@"${NEW_SURVEY,,}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".py" 2>/dev/null

mv $ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP,,}".yaml" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null 

sed --in-place --regexp-extended  "s@$OLD_PROJECT@"$NEW_PROJECT"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null

sed --in-place --regexp-extended  "s@${OLD_SURVEY,,}@"${NEW_SURVEY,,}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null

sed --in-place --regexp-extended  "s@${OLD_PROJECT^^}@"${NEW_PROJECT^^}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null

sed --in-place --regexp-extended  "s@${OLD_SURVEY}@"${NEW_SURVEY}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null

mv $ROOTDIR/projects/$NEW_PROJECT/likelihood/"params_"${TMP,,}".yaml" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/"params_"${TMP2,,}".yaml" 2>/dev/null

sed --in-place --regexp-extended  "s@${OLD_SURVEY}@"${NEW_SURVEY}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/"params_"${TMP2,,}".yaml" 2>/dev/null
# ------------------------------------------------------------------------------------

TMP=$OLD_SURVEY"_2x2pt"
TMP2=$NEW_SURVEY"_2x2pt"

mv $ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP,,}".py" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".py" 2>/dev/null

sed --in-place --regexp-extended  "s@$OLD_PROJECT@"$NEW_PROJECT"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".py" 2>/dev/null

sed --in-place --regexp-extended  "s@"${OLD_SURVEY,,}"@"${NEW_SURVEY,,}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".py" 2>/dev/null

mv $ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP,,}".yaml" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null

sed --in-place --regexp-extended  "s@$OLD_PROJECT@"$NEW_PROJECT"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null 

sed --in-place --regexp-extended  "s@${OLD_SURVEY,,}@"${NEW_SURVEY,,}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null 

sed --in-place --regexp-extended  "s@${OLD_PROJECT^^}@"${NEW_PROJECT^^}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null

sed --in-place --regexp-extended  "s@${OLD_SURVEY}@"${NEW_SURVEY}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null
# ------------------------------------------------------------------------------------

TMP=$OLD_SURVEY"_cosmic_shear"
TMP2=$NEW_SURVEY"_cosmic_shear"

mv $ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP,,}".py" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".py" 2>/dev/null

sed --in-place --regexp-extended  "s@$OLD_PROJECT@"$NEW_PROJECT"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".py" 2>/dev/null

sed --in-place --regexp-extended  "s@"${OLD_SURVEY,,}"@"${NEW_SURVEY,,}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".py" 2>/dev/null

mv $ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP,,}".yaml" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null

sed --in-place --regexp-extended  "s@$OLD_PROJECT@"$NEW_PROJECT"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null

sed --in-place --regexp-extended  "s@${OLD_SURVEY,,}@"${NEW_SURVEY,,}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null

sed --in-place --regexp-extended  "s@${OLD_PROJECT^^}@"${NEW_PROJECT^^}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null

sed --in-place --regexp-extended  "s@${OLD_SURVEY}@"${NEW_SURVEY}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null

mv $ROOTDIR/projects/$NEW_PROJECT/likelihood/"params_"${TMP,,}".yaml" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/"params_"${TMP2,,}".yaml" 2>/dev/null

sed --in-place --regexp-extended  "s@${OLD_SURVEY}@"${NEW_SURVEY}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/"params_"${TMP2,,}".yaml" 2>/dev/null

# ------------------------------------------------------------------------------------

TMP=$OLD_SURVEY"_xi_ggl"
TMP2=$NEW_SURVEY"_xi_ggl"

mv $ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP,,}".py" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".py" 2>/dev/null

sed --in-place --regexp-extended  "s@$OLD_PROJECT@"$NEW_PROJECT"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".py" 2>/dev/null

sed --in-place --regexp-extended  "s@"${OLD_SURVEY,,}"@"${NEW_SURVEY,,}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".py" 2>/dev/null

mv $ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP,,}".yaml" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null

sed --in-place --regexp-extended  "s@$OLD_PROJECT@"$NEW_PROJECT"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null

sed --in-place --regexp-extended  "s@${OLD_SURVEY,,}@"${NEW_SURVEY,,}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null

sed --in-place --regexp-extended  "s@${OLD_PROJECT^^}@"${NEW_PROJECT^^}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null

sed --in-place --regexp-extended  "s@${OLD_SURVEY}@"${NEW_SURVEY}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null

# ------------------------------------------------------------------------------------

TMP=$OLD_SURVEY"_ggl"
TMP2=$NEW_SURVEY"_ggl"

mv $ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP,,}".py" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".py" 2>/dev/null

sed --in-place --regexp-extended  "s@$OLD_PROJECT@"$NEW_PROJECT"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".py" 2>/dev/null

sed --in-place --regexp-extended  "s@"${OLD_SURVEY,,}"@"${NEW_SURVEY,,}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".py" 2>/dev/null

mv $ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP,,}".yaml" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null

sed --in-place --regexp-extended  "s@$OLD_PROJECT@"$NEW_PROJECT"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null

sed --in-place --regexp-extended  "s@${OLD_SURVEY,,}@"${NEW_SURVEY,,}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null

sed --in-place --regexp-extended  "s@${OLD_PROJECT^^}@"${NEW_PROJECT^^}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null

sed --in-place --regexp-extended  "s@${OLD_SURVEY}@"${NEW_SURVEY}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null

# ------------------------------------------------------------------------------------

TMP=$OLD_SURVEY"_clustering"
TMP2=$NEW_SURVEY"_clustering"

mv $ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP,,}".py" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".py" 2>/dev/null

sed --in-place --regexp-extended  "s@$OLD_PROJECT@"$NEW_PROJECT"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".py" 2>/dev/null

sed --in-place --regexp-extended  "s@"${OLD_SURVEY,,}"@"${NEW_SURVEY,,}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".py" 2>/dev/null

mv $ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP,,}".yaml" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null 

sed --in-place --regexp-extended  "s@$OLD_PROJECT@"$NEW_PROJECT"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null

sed --in-place --regexp-extended  "s@${OLD_SURVEY,,}@"${NEW_SURVEY,,}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null

sed --in-place --regexp-extended  "s@${OLD_PROJECT^^}@"${NEW_PROJECT^^}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null

sed --in-place --regexp-extended  "s@${OLD_SURVEY}@"${NEW_SURVEY}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/likelihood/${TMP2,,}".yaml" 2>/dev/null

# ------------------------------------------------------------------------------------

TMP=${OLD_PROJECT^^}
TMP2=${NEW_PROJECT^^}

mv $ROOTDIR/projects/$NEW_PROJECT/data/$TMP".dataset" \
$ROOTDIR/projects/$NEW_PROJECT/data/$TMP2".dataset" 2>/dev/null

rm -r $ROOTDIR/projects/$NEW_PROJECT/data/CosmoSiS 2>/dev/null
rm $ROOTDIR/projects/$NEW_PROJECT/data/*.txt 2>/dev/null
rm $ROOTDIR/projects/$NEW_PROJECT/data/*.mask 2>/dev/null

# ------------------------------------------------------------------------------------

TMP='EXAMPLE_EVALUATE1'

sed --in-place --regexp-extended  "s@${OLD_PROJECT}@"${NEW_PROJECT}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/$TMP".sbatch" 2>/dev/null

sed --in-place --regexp-extended  "s@${OLD_PROJECT}@"${NEW_PROJECT}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/$TMP".yaml" 2>/dev/null

sed --in-place --regexp-extended  "s@${OLD_PROJECT^^}@"${NEW_PROJECT^^}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/$TMP".yaml" 2>/dev/null

sed --in-place --regexp-extended  "s@${OLD_SURVEY}@"${NEW_SURVEY}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/$TMP".yaml" 2>/dev/null

# ------------------------------------------------------------------------------------

TMP='example_mcmc'
sed --in-place --regexp-extended  "s@${OLD_PROJECT}@"${NEW_PROJECT}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/examples/SLURM/$TMP".sbatch" 2>/dev/null

TMP='example_poly_mcmc'
sed --in-place --regexp-extended  "s@${OLD_PROJECT}@"${NEW_PROJECT}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/examples/SLURM/$TMP".sbatch" 2>/dev/null

TMP='example_poly'
sed --in-place --regexp-extended  "s@${OLD_PROJECT}@"${NEW_PROJECT}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/examples/SLURM/$TMP".sbatch" 2>/dev/null

TMP='example_post'
sed --in-place --regexp-extended  "s@${OLD_PROJECT}@"${NEW_PROJECT}"@g" \
$ROOTDIR/projects/$NEW_PROJECT/examples/SLURM/$TMP".sbatch" 2>/dev/null

# ------------------------------------------------------------------------------------
