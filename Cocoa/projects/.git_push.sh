if [ -z "${ROOTDIR}" ]; then
    echo 'ERROR ROOTDIR not define'
    return
fi

cd $ROOTDIR/projects

for NAME in $(find . -mindepth 1 -maxdepth 1 -type d ! -name 'example'); do
  declare NAME2=$(echo "${NAME}" | sed -E 's@./@@')
  if [ ! -d "$ROOTDIR/projects/$NAME2/" ]; then
    echo "Error: directory $ROOTDIR/projects/$NAME2/ DOES NOT exists."
  else
    cd $ROOTDIR/projects/$NAME2/
    git push origin main
  fi
done

unset NAME2
unset FNAME

cd $ROOTDIR