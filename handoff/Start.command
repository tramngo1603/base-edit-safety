#!/bin/sh
set -eu
export PATH=/usr/bin:/bin:/usr/sbin:/sbin
export LANG=en_US.UTF-8
unset PYTHONPATH PYTHONHOME VIRTUAL_ENV
umask 077
BES_BUNDLE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)
cd "$BES_BUNDLE"
fail() {
    echo "Setup needs attention: $1"
    echo "Read QUICKSTART.html in this folder for the next step."
    if [ -t 0 ]; then printf 'Press Return to close. '; read -r BES_REPLY; fi
    exit 1
}
[ "$(uname -m)" = arm64 ] || fail 'Use an Apple silicon Mac (arm64), without Rosetta.'
[ "$(uname -s)" = Darwin ] || fail 'This handoff is for macOS.'
echo 'Checking the handoff files…'
shasum -a 256 -c SHA256SUMS > /dev/null || fail 'A handoff file is missing or changed. Extract a fresh copy of the bundle.'
if [ ! -x .runtime/python/bin/python3 ]; then
    BES_FREE_KB=$(df -k . | awk 'END {print $4}')
    [ "$BES_FREE_KB" -ge 2097152 ] || fail 'Free at least 2 GiB in this folder for setup.'
    echo 'Installing the included Python 3.12.12 runtime…'
    mkdir -p .runtime
    tar -xzf vendor/python-runtime.tar.gz -C .runtime
fi
if [ -f .installed-at ] && [ "$(cat .installed-at)" != "$BES_BUNDLE" ]; then
    mv .venv ".venv-before-move-$(date +%s)"
    rm .installed-at
fi
if [ ! -f .installed-at ]; then
    echo 'Installing pinned dependencies from the included wheels (offline)…'
    .runtime/python/bin/python3 -I -m venv .venv
    .venv/bin/python -I -m pip --isolated install --no-index --require-hashes --find-links vendor/wheels -r requirements.lock
    .venv/bin/python -I -m pip --isolated check
    .venv/bin/python -I -c 'from base_edit_safety.launcher import configure_tk; configure_tk(); import tkinter; root = tkinter.Tk(); root.withdraw(); root.update(); root.destroy()'
    echo 'Testing the native Cas-OFFinder backend and the Apple OpenCL GPU…'
    BES_CHECK="setup-check-$(date +%Y%m%d-%H%M%S)"
    .venv/bin/python -I -m base_edit_safety.launcher --bundle "$BES_BUNDLE" --run small --output "$BES_BUNDLE/$BES_CHECK" \
        > "$BES_BUNDLE/$BES_CHECK.log" 2>&1 || fail "Backend check failed. Read $BES_CHECK.log; the Python wheel alone is not a working backend."
    printf '%s\n' "$BES_BUNDLE" > .installed-at
fi
if [ "${1:-}" = --setup-only ]; then
    echo 'Setup complete; the real native search and sample example passed.'
    exit 0
fi
echo 'Opening Base Edit Safety. Keep this terminal open while the app runs.'
exec .venv/bin/python -I -m base_edit_safety.launcher --bundle "$BES_BUNDLE"
