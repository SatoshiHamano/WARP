#!/bin/sh
set -eu

PYTHON=${PYTHON:-python3}
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
TMPBASE=${TMPDIR:-/tmp}
WORKDIR=$(mktemp -d "$TMPBASE/warp-sci-wide.XXXXXX")

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
    -d "$WORKDIR/4_Ari_WIDE_test" \
    -f

test -f "$WORKDIR/4_Ari_WIDE_test/calibration_data/input_files.txt"
test -f "$WORKDIR/4_Ari_WIDE_test/reduction_log/status.txt"

echo "Warp_sci WIDE fast smoke test finished successfully."
