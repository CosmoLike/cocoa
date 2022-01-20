if [ -z "${ROOTDIR}" ]; then
    echo 'ERROR ROOTDIR not defined'
    exit 1
fi

cd $ROOTDIR/projects

for NAME in $(find . -mindepth 1 -maxdepth 1 -type d ! -name 'example'); do
  declare NAME2=$(echo "${NAME}" | sed -E 's@./@@')
  rm -rf $ROOTDIR/projects/$NAME2
done

unset NAME2

cd $ROOTDIR