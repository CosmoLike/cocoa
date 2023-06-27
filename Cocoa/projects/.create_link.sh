for NAME in $(find . -mindepth 1 -maxdepth 1 -type d ! -name 'example'); do
  declare NAME2=$(echo "${NAME}" | sed -E 's@./@@')

  if [ ! -d "$ROOTDIR/projects/$NAME2/interface" ]; then
    echo "Warning: directory $ROOTDIR/projects/$NAME2/interface DOES NOT exists."
  else
    if [ ! -d "$ROOTDIR/external_modules/code/$NAME2" ]; then
      ln -sf $ROOTDIR/projects/$NAME2/interface $ROOTDIR/external_modules/code/
      mv -T $ROOTDIR/external_modules/code/interface $ROOTDIR/external_modules/code/$NAME2
    fi
  fi

  if [ ! -d "$ROOTDIR/projects/$NAME2/likelihood" ]; then
    echo "Warning: directory $ROOTDIR/projects/$NAME2/likelihood DOES NOT exists."
  else
    if [ ! -d "$ROOTDIR/cobaya/cobaya/likelihoods/$NAME2" ]; then
      ln -sf $ROOTDIR/projects/$NAME2/likelihood $ROOTDIR/cobaya/cobaya/likelihoods/
      mv -T $ROOTDIR/cobaya/cobaya/likelihoods/likelihood $ROOTDIR/cobaya/cobaya/likelihoods/$NAME2
    fi
  fi

  if [ ! -d "$ROOTDIR/projects/$NAME2/data" ]; then
    echo "Warning: directory $ROOTDIR/projects/$NAME2/data DOES NOT exists."
  else
    if [ ! -d "$ROOTDIR/external_modules/data/$NAME2" ]; then
      ln -sf $ROOTDIR/projects/$NAME2/data $ROOTDIR/external_modules/data/
      mv -T $ROOTDIR/external_modules/data/data $ROOTDIR/external_modules/data/$NAME2
    fi
  fi
done

unset NAME2