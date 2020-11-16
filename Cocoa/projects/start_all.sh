cd $ROOTDIR/projects

# create tmp .gitignore and include all projects in it
touch $ROOTDIR/projects/.gitignore

for NAME in $(find . -mindepth 1 -maxdepth 1 -type d ! -name 'example'); do
  cd $ROOTDIR/projects/$NAME/scripts
  # https://stackoverflow.com/a/11231970/2472169
  # we don't want the sh to crash if source
  declare FNAME=start_$(echo "${NAME}" | sed -E 's@./@@')
  if [ -f "$FNAME" ]; then
    source $FNAME
  fi
  cd $ROOTDIR/projects
  declare NAME2=$(echo "${NAME}" | sed -E 's@./@@')
  echo "${NAME2}" >> $ROOTDIR/projects/.gitignore
done

source $ROOTDIR/projects/create_link.sh

cp $ROOTDIR/projects/.gitignore $ROOTDIR/external_modules/data/
cp $ROOTDIR/projects/.gitignore $ROOTDIR/external_modules/code/
cp $ROOTDIR/projects/.gitignore $ROOTDIR/cobaya/cobaya/likelihoods

#cd $ROOTDIR