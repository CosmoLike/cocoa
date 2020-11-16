cd $ROOTDIR/projects

for NAME in $(find . -mindepth 1 -maxdepth 1 -type d ! -name 'example'); do
  cd $ROOTDIR/projects/$NAME/scripts
  # https://stackoverflow.com/a/11231970/2472169
  # we don't want the sh to crash if source
  declare FNAME=setup_$(echo "${NAME}" | sed -E 's@./@@')
  echo $FNAME
  if [ -f "$FNAME" ]; then
    source $FNAME
  fi
done

unset $FNAME

cd $ROOTDIR