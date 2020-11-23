if [ -z "${ROOTDIR}" ]; then
    echo 'ERROR ROOTDIR not defined'
    return
fi

git push origin main
source $ROOTDIR/projects/.git_push.sh