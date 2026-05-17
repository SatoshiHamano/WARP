#!/bin/sh
set -eu

PYTHON=${PYTHON:-python3}
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
TMPBASE=${TMPDIR:-/tmp}
WORKDIR=$(mktemp -d "$TMPBASE/warp-sci-full.XXXXXX")

export PYRAF_NO_DISPLAY=${PYRAF_NO_DISPLAY:-1}
if [ -n "${WARP_IRAFARCH:-}" ]; then
    IRAFARCH=$WARP_IRAFARCH
elif [ -z "${IRAFARCH:-}" ] || [ "${IRAFARCH:-}" = "macintel" ]; then
    IRAFARCH=macos64
fi
export IRAFARCH

cleanup() {
    if [ "${KEEP_WARP_TEST_OUTPUT:-0}" = "1" ]; then
        echo "Keeping test output in $WORKDIR"
    else
        rm -rf "$WORKDIR"
    fi
}
trap cleanup EXIT

"$PYTHON" "$ROOT/Warp_sci.py" \
    "$ROOT/TEST/WIDE/4_Ari_list.txt" \
    -r "$ROOT/TEST/WIDE/" \
    -c "$ROOT/TEST/WIDE/WINERED_calibration_LCO22b_WIDE100_20220914_v2/" \
    -v "$ROOT/TEST/WIDE/" \
    -d "$WORKDIR/4_Ari_WIDE_test_def"

"$PYTHON" "$ROOT/Warp_sci.py" \
    "$ROOT/TEST/WIDE/4_Ari_list.txt" \
    -r "$ROOT/TEST/WIDE/" \
    -c "$ROOT/TEST/WIDE/WINERED_calibration_LCO22b_WIDE100_20220914_v2/" \
    -v "$ROOT/TEST/WIDE/" \
    -d "$WORKDIR/4_Ari_WIDE_test_skysub" \
    -p "$ROOT/TEST/WIDE/paramSample.txt"

"$PYTHON" "$ROOT/Warp_sci.py" \
    "$ROOT/TEST/HIRES-J/21_Peg_list.txt" \
    -r "$ROOT/TEST/HIRES-J/" \
    -c "$ROOT/TEST/HIRES-J/WINERED_calibration_NTTrun_HIRESJ_20170728/" \
    -v "$ROOT/TEST/HIRES-J/" \
    -d "$WORKDIR/21_Peg_HIRESJ_test"

"$PYTHON" "$ROOT/Warp_sci.py" \
    "$ROOT/TEST/HIRES-Y/HD_163336_list.txt" \
    -r "$ROOT/TEST/HIRES-Y/" \
    -c "$ROOT/TEST/HIRES-Y/WINERED_calibration_NTTrun_HIRESY_20170727/" \
    -v "$ROOT/TEST/HIRES-Y/" \
    -d "$WORKDIR/HD_163336_HIRESY_test"

echo "Full Warp_sci smoke tests finished successfully."
