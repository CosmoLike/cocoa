if [ -z "${ROOTDIR}" ]; then
    echo 'ERROR ROOTDIR not defined'
    return
fi

git add --all
source $ROOTDIR/projects/.git_add.sh