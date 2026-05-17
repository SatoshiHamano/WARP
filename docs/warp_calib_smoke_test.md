# Warp_calib Smoke Test

This document records the first lightweight execution coverage for
`Warp_calib.py`.

## Scope

The smoke test runs the WIDE calibration fixture:

- input directory: `TEST/WIDE-calib/`
- input list: `input.list`
- command target: `Warp_calib.py`

The test checks that calibration output is created, but it does not yet compare
numerical output values against a reference baseline.

## Script

Use:

```sh
./testWarpCalib.sh
```

The script:

- copies `TEST/WIDE-calib/` into a temporary working directory
- copies `$HOME/login.cl` into that working directory when available
- creates a local `uparm/` directory
- sets `PYRAF_NO_DISPLAY=1` by default
- maps an inherited `IRAFARCH=macintel` to `IRAFARCH=macos64`
- runs `Warp_calib.py input.list`
- verifies that the expected calibration products exist

To choose a Python interpreter explicitly:

```sh
PYTHON=/path/to/python iraf=/path/to/iraf/ ./testWarpCalib.sh
```

To keep the temporary output directory for inspection:

```sh
KEEP_WARP_TEST_OUTPUT=1 PYTHON=/path/to/python iraf=/path/to/iraf/ ./testWarpCalib.sh
```

## Verified Locally

The WIDE calibration smoke test completed successfully in:

- Python 3.7.16 Astroconda/PyRAF environment
- Python 3.13.5 test environment with PyRAF 2.2.4

Both runs used:

```sh
iraf=/Users/hamano/iraf/iraf-2.17.1/
IRAFARCH=macos64
```

## Environment Notes

`Warp_calib.py` is more sensitive to IRAF initialization than the current
`Warp_sci.py` smoke test.

Without a local `login.cl` and `uparm/`, the Astroconda/PyRAF run reached
`apscatter` and failed because IRAF variable `uparm` was undefined.

With a stale inherited `IRAFARCH=macintel`, PyRAF looked for:

```text
bin.macintel/x_system.e
```

The working Community IRAF installation provides the executable under:

```text
bin.macos64/x_system.e
```

The smoke script handles this local macOS case by defaulting to `macos64`.

## Remaining Work

- Add numerical regression summaries for the generated calibration products.
- Add HIRES-Y and HIRES-J calibration smoke coverage.
- Refactor `Warp_calib.py` only after the calibration baseline is stronger.
