cd $ROOTDIR/projects

# create tmp .gitignore and include all projects in it
touch $ROOTDIR/projects/.gitignore

for NAME in $(find . -mindepth 1 -maxdepth 1 -type d ! -name 'example'); do
  cd $ROOTDIR/projects/$NAME/scripts
  # https://stackoverflow.com/a/11231970/2472169
  # we don't want the sh to crash if source
  echo $NAME
  declare FNAME=start_$NAME
  if [ -f "$FNAME" ]; then
    source $FNAME
  fi
  cd $ROOTDIR/projects
  echo "${NAME}/" >> $ROOTDIR/projects/.gitignore
done

source $ROOTDIR/projects/create_link.sh

cd $ROOTDIR