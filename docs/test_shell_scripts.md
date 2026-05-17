# Test Shell Scripts

This document summarizes the repository-level shell smoke tests.

## Common Behavior

The shell scripts can be run from the repository root:

```sh
./testWarpSci.sh
./testWarpSciFull.sh
./testWarpCalib.sh
```

They use paths relative to the script location, so they do not depend on the
caller manually typing `./TEST/...` paths.

The scripts set `PYRAF_NO_DISPLAY=1` by default for non-interactive PyRAF runs.
On local macOS runs, they also map an inherited `IRAFARCH=macintel` to
`IRAFARCH=macos64`.  Override this with `WARP_IRAFARCH` if needed.

By default, generated outputs are written under a temporary directory and
removed at the end of the run.  To keep the outputs for inspection:

```sh
KEEP_WARP_TEST_OUTPUT=1 ./testWarpSci.sh
```

To select a Python interpreter:

```sh
PYTHON=/path/to/python ./testWarpSci.sh
```

For PyRAF/IRAF-backed runs, also set `iraf` when it is not already available in
the shell:

```sh
PYTHON=/path/to/python iraf=/path/to/iraf/ ./testWarpCalib.sh
```

## Scripts

- `testWarpSci.sh`
  - WIDE science fast smoke test.
  - Uses `TEST/WIDE/4_Ari_list.txt`.
  - Runs with `-f`.
- `testWarpSciFull.sh`
  - Longer science smoke suite.
  - Runs WIDE, WIDE with sky-subtraction parameters, HIRES-J, and HIRES-Y.
  - Intended for manual validation, not quick review.
- `testWarpCalib.sh`
  - WIDE calibration smoke test.
  - Uses `TEST/WIDE-calib/input.list`.
  - Copies the fixture into a temporary working directory because
    `Warp_calib.py` creates many intermediate products in the current
    directory.

## Current Limits

These scripts are smoke tests.  They check that representative pipeline runs
complete and that selected output files exist.  They do not yet compare
numerical output values against reference summaries.
