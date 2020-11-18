if [ -z "${ROOTDIR}" ]; then
    echo 'ERROR ROOTDIR not define'
    return
fi

cd $ROOTDIR/projects

for NAME in $(find . -mindepth 1 -maxdepth 1 -type d ! -name 'example'); do
  declare NAME2=$(echo "${NAME}" | sed -E 's@./@@')
  if [ ! -d "$ROOTDIR/projects/$NAME2/scripts" ]; then
    echo "Error: directory $ROOTDIR/projects/$NAME2/scripts DOES NOT exists."
  else
    cd $ROOTDIR/projects/$NAME2/scripts
    # https://stackoverflow.com/a/11231970/2472169
    # we don't want the sh to crash if source
    declare FNAME=stop_$NAME2
    if [ -f "$FNAME" ]; then
      source $FNAME
    fi
    cd $ROOTDIR/projects
  fi
done

unset FNAME
unset NAME2

source $ROOTDIR/projects/.delete_link.sh

rm -f $ROOTDIR/projects/.gitignore
rm -f $ROOTDIR/external_modules/data/.gitignore
rm -f $ROOTDIR/external_modules/code/.gitignore
rm -f $ROOTDIR/cobaya/cobaya/likelihoods/.gitignore

cd $ROOTDIR