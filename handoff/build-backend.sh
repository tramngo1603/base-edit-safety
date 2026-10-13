#!/bin/sh
set -eu
export PATH=/usr/bin:/bin:/usr/sbin:/sbin
unset PYTHONPATH PYTHONHOME VIRTUAL_ENV CPATH LIBRARY_PATH CMAKE_PREFIX_PATH
BES_BUNDLE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)
cd "$BES_BUNDLE"
[ "$#" -eq 1 ] || { echo 'Usage: ./build-backend.sh NEW_BUILD_FOLDER'; exit 2; }
BES_BUILD="$1"
[ ! -e "$BES_BUILD" ] || { echo 'Choose a new build folder.'; exit 2; }
/usr/bin/xcode-select -p >/dev/null || { echo 'Install Apple Command Line Tools with xcode-select --install, then retry.'; exit 2; }
shasum -a 256 -c SHA256SUMS >/dev/null
mkdir -p "$BES_BUILD"
BES_BUILD=$(CDPATH= cd -- "$BES_BUILD" && pwd -P)
.venv/bin/python -I -m pip --isolated install --no-index --require-hashes --find-links vendor/wheels -r build-requirements.lock
tar -xzf vendor/cas-offinder-source.tar.gz -C "$BES_BUILD"
BES_SOURCE="$BES_BUILD/cas-offinder-9816b94c20c4cba2e79b039e1e2a6dee684b7b66"
(cd "$BES_SOURCE" && /usr/bin/patch -p1 -i "$BES_BUNDLE/cas-offinder-2.4.1-compat.patch")
.venv/bin/python -I -m cmake -S "$BES_SOURCE" -B "$BES_BUILD/build" \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
    -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ \
    -DCMAKE_OSX_DEPLOYMENT_TARGET=14.4 -DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=TRUE
.venv/bin/python -I -m cmake --build "$BES_BUILD/build" --parallel 2
shasum -a 256 "$BES_BUILD/build/cas-offinder"
/usr/bin/otool -L "$BES_BUILD/build/cas-offinder"
echo 'Build complete. Compare its hash with backend.sha256 in RELEASE.json.'
