for NAME in $(find . -mindepth 1 -maxdepth 1 -type d ! -name 'example'); do
  declare NAME2=$(echo "${NAME}" | sed -E 's@./@@')

  rm -f $ROOTDIR/external_modules/code/interface
  rm -f $ROOTDIR/external_modules/code/$NAME2
  if [ ! -d "$ROOTDIR/projects/$NAME2/interface" ]; then
    echo "Error: directory $ROOTDIR/projects/$NAME2/interface DOES NOT exists."
  else
    ln -sf $ROOTDIR/projects/$NAME2/interface $ROOTDIR/external_modules/code/
    mv $ROOTDIR/external_modules/code/interface $ROOTDIR/external_modules/code/$NAME2
  fi

  rm -f $ROOTDIR/cobaya/cobaya/likelihoods/likelihood
  rm -f $ROOTDIR/cobaya/cobaya/likelihoods/$NAME2
  if [ ! -d "$ROOTDIR/projects/$NAME2/likelihood" ]; then
    echo "Error: directory $ROOTDIR/projects/$NAME2/likelihood DOES NOT exists."
  else
    ln -sf $ROOTDIR/projects/$NAME2/likelihood $ROOTDIR/cobaya/cobaya/likelihoods/
    mv $ROOTDIR/cobaya/cobaya/likelihoods/likelihood $ROOTDIR/cobaya/cobaya/likelihoods/$NAME2
  fi

  rm -f $ROOTDIR/external_modules/data/data
  rm -f $ROOTDIR/external_modules/data/$NAME2
  if [ ! -d "$ROOTDIR/projects/$NAME2/data" ]; then
    echo "Error: directory $ROOTDIR/projects/$NAME2/data DOES NOT exists."
  else
    ln -sf $ROOTDIR/projects/$NAME2/data $ROOTDIR/external_modules/data/
    mv $ROOTDIR/external_modules/data/data $ROOTDIR/external_modules/data/$NAME2
  fi
done

unset NAME2