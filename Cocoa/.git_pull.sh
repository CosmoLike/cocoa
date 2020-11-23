if [ -z "${ROOTDIR}" ]; then
    echo 'ERROR ROOTDIR not defined'
    return
fi

git pull origin main
source $ROOTDIR/projects/.git_pull.sh