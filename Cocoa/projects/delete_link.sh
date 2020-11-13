cd $ROOTDIR/projects

for NAME in $(find . -mindepth 1 -maxdepth 1 -type d); do
  rm $ROOTDIR/external_modules/code/$NAME

  rm $ROOTDIR/cobaya/cobaya/likelihoods/$NAME

  rm $ROOTDIR/external_modules/data/$NAME
done

cd $ROOTDIR