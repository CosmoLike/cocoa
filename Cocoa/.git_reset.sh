if [ -z "${ROOTDIR}" ]; then
    echo 'ERROR ROOTDIR not defined'
    exit 1
fi

git reset --soft HEAD~1
sh $ROOTDIR/projects/.git_reset.sh