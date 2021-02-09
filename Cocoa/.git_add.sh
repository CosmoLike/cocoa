if [ -z "${ROOTDIR}" ]; then
    echo 'ERROR ROOTDIR not defined'
    exit 1
fi

git add --all
source $ROOTDIR/projects/.git_add.sh