if [ -z "${ROOTDIR}" ]; then
    echo 'ERROR ROOTDIR not define'
    exit 1
fi

cd $ROOTDIR/projects

for NAME in $(find . -mindepth 1 -maxdepth 1 -type d ! -name 'example'); do
  declare NAME2=$(echo "${NAME}" | sed -E 's@./@@')
  if [ ! -d "$ROOTDIR/projects/$NAME2/scripts" ]; then
    echo "Error: directory $ROOTDIR/projects/$NAME2/scripts DOES NOT exists."
  else
    cd $ROOTDIR/projects/$NAME2/scripts
    echo $ROOTDIR/projects/$NAME2/scripts
    # https://stackoverflow.com/a/11231970/2472169
    # we don't want the sh to crash if source
    declare FNAME=compile_$NAME2
    if [ -f "$FNAME" ]; then
      sh $FNAME
    fi
  fi
done

unset NAME2
unset FNAME

cd $ROOTDIR