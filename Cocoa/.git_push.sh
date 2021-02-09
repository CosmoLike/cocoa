if [ -z "${ROOTDIR}" ]; then
    echo 'ERROR ROOTDIR not defined'
    exit 1
fi

git push origin main
source $ROOTDIR/projects/.git_push.sh