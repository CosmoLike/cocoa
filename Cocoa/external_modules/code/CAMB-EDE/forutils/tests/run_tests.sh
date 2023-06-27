
gfortran --version


cd forutils/tests

set -e
set -o pipefail

make

echo ''
echo 'Testing Debug build..'

./tester_Debug

echo ''
echo 'Testing Release build..'

./tester_Release

