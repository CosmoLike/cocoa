if [ -z "${ROOTDIR}" ]; then
    echo 'ERROR ROOTDIR not defined'
    return
fi

git commit -am "update"
source $ROOTDIR/projects/.git_commit.sh