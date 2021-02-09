if [ -z "${ROOTDIR}" ]; then
    echo 'ERROR ROOTDIR not define'
    exit 1
fi

cd $ROOTDIR/projects

URLS="git@github.com:CosmoLike/cocoa_des_y3.git
      git@github.com:CosmoLike/cocoa_desxplanck.git
      git@github.com:CosmoLike/cocoa_lsst_fourier.git"

for NAME in $URLS; do
  # https://stackoverflow.com/a/11231970/2472169
  # we don't want the sh to crash if git fails
  rm -rf $NAME
  git clone $NAME || true
done

# take the cocoa_ out of the dir names
for DIR in $(find . -mindepth 1 -maxdepth 1 -type d ! -name 'example'); do
    # https://unix.stackexchange.com/a/61402
    mv -T "${DIR}" $(echo "${DIR}" | sed -E 's@cocoa_@@') 2> /dev/null
done

cd $ROOTDIR