cd $ROOTDIR/projects

for NAME in $(find . -mindepth 1 -maxdepth 1 -type d ! -name 'example'); do
  cd $ROOTDIR/projects/$NAME/scripts
  # https://stackoverflow.com/a/11231970/2472169
  # we don't want the sh to crash if source
  FNAME = stop_$NAME
  if [ -f "$FNAME" ]; then
    source $FNAME
  fi
done

rm -rf $ROOTDIR/projects/.gitignore

source delete_link.sh

cd $ROOTDIR