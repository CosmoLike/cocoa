if [ -z "${ROOTDIR}" ]; then
    echo 'ERROR ROOTDIR not defined'
    exit 1
fi

git update-index --assume-unchanged set_installation_options
git update-index --assume-unchanged setup_cocoa_installation_packages
git update-index --assume-unchanged compile_external_modules
git update-index --assume-unchanged start_cocoa
git update-index --assume-unchanged stop_cocoa
git update-index --assume-unchanged clean_all

git commit -am "update"
sh $ROOTDIR/projects/.git_commit.sh