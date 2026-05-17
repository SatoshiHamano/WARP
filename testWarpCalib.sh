#!/bin/sh
set -eu

PYTHON=${PYTHON:-python3}
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
TMPBASE=${TMPDIR:-/tmp}
WORKDIR=$(mktemp -d "$TMPBASE/warp-calib-wide.XXXXXX")

cleanup() {
    if [ "${KEEP_WARP_TEST_OUTPUT:-0}" = "1" ]; then
        echo "Keeping test output in $WORKDIR"
    else
        rm -rf "$WORKDIR"
    fi
}
trap cleanup EXIT

cp -R "$ROOT/TEST/WIDE-calib/." "$WORKDIR/"

if [ -f "$HOME/login.cl" ]; then
    cp "$HOME/login.cl" "$WORKDIR/login.cl"
fi
mkdir -p "$WORKDIR/uparm"

export PYRAF_NO_DISPLAY=${PYRAF_NO_DISPLAY:-1}
if [ -n "${WARP_IRAFARCH:-}" ]; then
    IRAFARCH=$WARP_IRAFARCH
elif [ -z "${IRAFARCH:-}" ] || [ "${IRAFARCH:-}" = "macintel" ]; then
    IRAFARCH=macos64
fi
export IRAFARCH

if [ -z "${iraf:-}" ]; then
    echo "Warning: iraf is not set. PyRAF/IRAF may fail unless it can find a default IRAF installation."
fi

(
    cd "$WORKDIR"
    "$PYTHON" "$ROOT/Warp_calib.py" input.list
)

test -f "$WORKDIR/calibration_LCO22b_setting3_WIDE100/input_files.txt"
test -f "$WORKDIR/calibration_LCO22b_setting3_WIDE100/flat_WIDE100_20220915_mscmn.fits"
test -f "$WORKDIR/calibration_LCO22b_setting3_WIDE100/mask_flat_WIDE100_20220915.fits"
test -f "$WORKDIR/calibration_LCO22b_setting3_WIDE100/comp_WIDE100_20220915_fm_ecall.fits"

echo "Warp_calib WIDE smoke test finished successfully."
