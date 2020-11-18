for NAME in $(find . -mindepth 1 -maxdepth 1 -type d ! -name 'example'); do
  declare NAME2=$(echo "${NAME}" | sed -E 's@./@@')

  rm -f $ROOTDIR/external_modules/code/$NAME2

  rm -f $ROOTDIR/cobaya/cobaya/likelihoods/$NAME2

  rm -f $ROOTDIR/external_modules/data/$NAME2
done

unset NAME2