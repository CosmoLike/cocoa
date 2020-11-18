if [ -z "${ROOTDIR}" ]; then
    echo 'ERROR ROOTDIR not define'
    return
fi

cd $ROOTDIR/projects

for NAME in $(find . -mindepth 1 -maxdepth 1 -type d ! -name 'example'); do
  cd $ROOTDIR/projects/$NAME/scripts
  # https://stackoverflow.com/a/11231970/2472169
  # we don't want the sh to crash if source
  declare FNAME=stop_$(echo "${NAME}" | sed -E 's@./@@')
  if [ -f "$FNAME" ]; then
    source $FNAME
  fi
  cd $ROOTDIR/projects
done

unset FNAME

source $ROOTDIR/projects/.delete_link.sh

rm -f $ROOTDIR/projects/.gitignore
rm -f $ROOTDIR/external_modules/data/.gitignore
rm -f $ROOTDIR/external_modules/code/.gitignore
rm -f $ROOTDIR/cobaya/cobaya/likelihoods/.gitignore

cd $ROOTDIR