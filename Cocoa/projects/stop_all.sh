cd $ROOTDIR/projects

for NAME in $(find . -mindepth 1 -maxdepth 1 -type d ! -name 'example'); do
  cd $ROOTDIR/projects/$NAME/scripts
  # https://stackoverflow.com/a/11231970/2472169
  # we don't want the sh to crash if source
  declare FNAME=stop_$NAME
  if [ -f "$FNAME" ]; then
    source $FNAME
  fi
  cd $ROOTDIR/projects
done

source delete_link.sh

rm -rf $ROOTDIR/external_modules/data/.gitignore
rm -rf $ROOTDIR/external_modules/code/.gitignore
rm -rf $ROOTDIR/cobaya/cobaya/likelihoods/.gitignore

cd $ROOTDIR