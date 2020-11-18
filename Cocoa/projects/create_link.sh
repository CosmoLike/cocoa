cd $ROOTDIR/projects

for NAME in $(find . -mindepth 1 -maxdepth 1 -type d ! -name 'example'); do
  declare NAME2=$(echo "${NAME}" | sed -E 's@./@@')

  ln -s $ROOTDIR/projects/$NAME/interface $ROOTDIR/external_modules/code/
  mv $ROOTDIR/external_modules/code/interface $ROOTDIR/external_modules/code/$NAME2

  ln -s $ROOTDIR/projects/$NAME/likelihood $ROOTDIR/cobaya/cobaya/likelihoods/
  mv $ROOTDIR/cobaya/cobaya/likelihoods/likelihood $ROOTDIR/cobaya/cobaya/likelihoods/$NAME2

  ln -s $ROOTDIR/projects/$NAME/data $ROOTDIR/external_modules/data/
  mv $ROOTDIR/external_modules/data/data $ROOTDIR/external_modules/data/$NAME2
done

cd $ROOTDIR