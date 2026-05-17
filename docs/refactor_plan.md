# WARP Refactor Plan

This document describes a low-risk modernization plan for the WINERED Automatic
Reduction Pipeline.

## Guiding Principles

WARP is already used by multiple people for real data reduction, so changes
should be small, reviewable, and separated by purpose.

- Behavior-preserving maintenance
- Test, packaging, and environment improvements
- Changes that may affect scientific results

Changes that may affect scientific results should be protected by reference
outputs or explicit before/after comparisons.

## Completed

- Fixed the local branch pull configuration.
  - `main` now pulls only from `origin/main`.
  - The earlier `config.py` conflict was caused by an untracked local
    `config.py` plus stale/duplicated upstream merge settings.
- Added `.gitignore`.
  - macOS/editor files
  - Python caches
  - IRAF/PyRAF local working files
  - generated documentation/output
  - root-level runtime `database/`
  - test-run output directories
- Removed tracked local artifacts.
  - `.DS_Store`
  - `._*`
  - `*~`
  - old `.pyc` files
- Added lightweight pytest coverage.
  - WIDE input-list parsing
  - WIDE 100-line long-input regression fixture
  - HIRES-J / HIRES-Y input-list parsing
  - required aperture fields when manual aperture mode is enabled
  - `TEST/WIDE/paramSample.txt` parameter parsing
  - `warp.config` import does not import PyRAF
  - `make_rawdata_list` import does not import PyRAF
- Updated science smoke-test shell scripts.
  - default to `python3`
  - allow selecting a Python executable, for example
    `PYTHON=python3.11 ./testWarpSci.sh`
  - stop on shell errors with `set -eu`
- Documented test commands in README.
- Isolated a small FITS-header helper from PyRAF-dependent code.
  - Moved `header_key_read` to `warp/fits_utils.py`.
  - Kept `warp.Spec2Dtools.header_key_read` available for compatibility.
  - `warp.config` can now be imported without importing PyRAF.
  - Updated direct `header_key_read` imports to use `warp.fits_utils`.
- Verified lightweight tests on newer Python versions.
  - Python 3.7.16 -> 8 passed
  - Python 3.12.11 -> 8 passed
  - Python 3.13.5 -> 8 passed
- Verified one full WIDE fast smoke run with the current Astroconda/PyRAF
  environment.
  - Explicit environment:
    `IRAFARCH=macos64 iraf=/Users/hamano/iraf/iraf-2.17.1/`
  - Command shape:
    `python Warp_sci.py ./TEST/WIDE/4_Ari_list.txt ... -d /private/tmp/warp_codex_4_Ari_WIDE_test -f`
  - Result: exited successfully and printed `=== Finished. ===`
  - Re-verified after direct `header_key_read` import cleanup.
- Created an isolated Python 3.13 PyRAF test environment without modifying the
  existing Astroconda environment.
  - Environment name: `warp-py313-pyraf`
  - Base environment: clone of `warp-py313`
  - Added package: `pyraf`
  - Verified PyRAF 2.2.4 import on Python 3.13.5.
  - Verified `iraf.noao()` and `iraf.onedspec()` package loading.
  - Verified lightweight pytest suite: 8 passed.
  - Verified the same WIDE fast smoke case with Python 3.13.5, PyRAF 2.2.4,
    and the existing IRAF 2.17.1 installation.
  - Required explicit non-interactive environment:
    `PYRAF_NO_DISPLAY=1 IRAFARCH=macos64 iraf=/Users/hamano/iraf/iraf-2.17.1/`

## Related Commits

- `83d9749 Add gitignore for local artifacts`
- `6623e6f Remove tracked local artifacts`
- `60ce277 Add lightweight input parsing tests`
- `a52abd8 Document refactor plan`
- `1f9e977 Document PyRAF environment caveat`

## Next Steps

1. Open a pull request for the current low-risk maintenance and
   PyRAF import-isolation changes.
2. Add a small environment definition.
   - Start by recording the currently working environment.
   - Treat newer Python targets as a separate follow-up.
3. Split tests into two clear layers.
   - Lightweight tests that do not require IRAF/PyRAF
   - Full smoke tests that require PyRAF/IRAF and the TEST data
4. Add a reference-output strategy for full pipeline runs.
   - Record expected output file lists.
   - Compare important FITS headers.
   - Compare array shapes and simple statistics.
   - Later, compare selected numerical outputs with tolerances.
5. Continue modernizing Python compatibility in small patches.
   - Identify remaining PyRAF side effects that happen at import time.
   - Separate IRAF-dependent processing from pure parsing/utilities.
   - Run lightweight tests on newer Python versions after each change.
6. Once test coverage is stronger, start behavior-preserving refactors.
   - Clean up import paths.
   - Clean up configuration parsing.
   - Clean up logging and error handling.
   - Split large functions.

## Planned Output Layout Refactor

The WARP output directory tree is scientifically important because many users
may already have scripts that consume the generated products.  For that reason,
changing the actual output paths should not be an early refactor.

The current low-risk step is to keep the existing output paths unchanged while
collecting the output directory names in code.  A better later step is to
replace ad-hoc path string construction with a small layout object, for example
`ReductionLayout`, that returns paths for each output class.

Candidate first scope:

- 1D science spectra only.
  - `*_NO*/onedspec/{AIR_flux,AIR_norm,AIR_cont,VAC_flux,VAC_norm,VAC_cont}/fsr*`
  - `*_sum/{AIR_flux,AIR_norm,AIR_cont,VAC_flux,VAC_norm,VAC_cont}/fsr*`
- Keep generated paths byte-for-byte compatible with the current output.
- Add unit tests that assert the layout object returns the legacy paths.
- Use the 1D spectrum regression summary to confirm that a real run is
  unchanged.

Possible design:

- `warp.output_layout.ReductionLayout`
  - owns `Path`-based output path construction.
  - formats `fsr` directories in one place.
  - separates frame-level and combined-output paths.
  - avoids scattering literal directory names and `fsr%.2f` formatting across
    `Warp_sci.py`.

Priority:

- Useful, but not urgent.
- Do this after the current stacked test/import/package-boundary PRs are
  reviewed.
- Prefer a dedicated PR because output paths are part of the user-facing
  contract, even if the intended change is behavior-preserving.

## Planned Numerical Regression Tests

The first numerical regression layer should avoid committing large generated
FITS products. Instead, it should commit compact summaries derived from a known
good pipeline run.

Initial target:

- WIDE fast smoke case:
  - input list: `TEST/WIDE/4_Ari_list.txt`
  - calibration directory:
    `TEST/WIDE/WINERED_calibration_LCO22b_WIDE100_20220914_v2/`
  - command style:
    `python Warp_sci.py ./TEST/WIDE/4_Ari_list.txt ... -f`
- Run this only in a PyRAF/IRAF-capable environment.
- Keep it separate from the default lightweight pytest suite.

Proposed files:

- `tools/summarize_warp_output.py`
  - Reads a WARP output directory.
  - Produces a deterministic JSON summary.
- `tests/reference/wide_4_ari_fast_summary.json`
  - Compact reference summary generated from the current known-good output.
- `tests/test_reference_outputs.py`
  - Compares a newly generated summary with the committed reference.
  - Mark this test as `regression` or `smoke` so it is opt-in.

Suggested summary content for selected FITS outputs:

- output file list
- FITS HDU count
- selected stable header values
  - examples: `NAXIS*`, `CRVAL*`, `CRPIX*`, `CDELT*`, `AIRORVAC`
- array shape and dtype
- numerical fingerprints
  - `nanmin`
  - `nanmax`
  - `nanmean`
  - `nanmedian`
  - `nanstd`
  - sampled `nansum`, for example every 50th pixel
  - a few fixed pixel or small-window values

Comparison policy:

- file lists and array shapes should match exactly.
- stable string headers should match exactly.
- floating point values should use explicit `rtol` and `atol`.
- volatile headers and run-specific paths should be excluded.
- the first implementation should compare summaries only, not the full FITS
  files.

Execution policy:

- Default command remains lightweight:
  - `python -m pytest`
- Numerical regression command is opt-in:
  - example: `python -m pytest -m regression`
- If the full pipeline execution is too slow for regular review, split it into
  two manual steps:
  - run WARP and generate a summary JSON
  - run pytest to compare that summary with the committed reference

## Known Risks / Notes

- PyRAF/IRAF behavior depends on shell environment variables, not only on the
  current working directory or `login.cl`.
  In the normal interactive terminal, PyRAF worked with:
  - `iraf=/Users/hamano/iraf/iraf-2.17.1/`
  - `IRAFARCH=macos64`
  In one Codex shell session, `IRAFARCH` was `macintel`, so PyRAF looked for
  `bin.macintel/x_system.e` and failed with:
  - `Cannot find executable for task pathnames`
  The correct executable exists under `bin.macos64/x_system.e`.
  When running PyRAF/WARP tests from non-interactive tools, explicitly set:
  - `IRAFARCH=macos64`
  - `iraf=/Users/hamano/iraf/iraf-2.17.1/`
- On Python 3.13 with PyRAF 2.2.4 on macOS, PyRAF import failed in one
  non-interactive shell unless `PYRAF_NO_DISPLAY=1` was set.
  This does not indicate that WARP itself has a GUI dependency.  It came from
  PyRAF's macOS display/focus helper during import.  WARP is still operated as a
  command-line pipeline, but PyRAF may need this variable in headless or
  tool-driven test runs.
- `warp.config` previously imported `warp.Spec2Dtools`, and
  `warp.Spec2Dtools` imports PyRAF at module import time.
  This made simple config tests depend on a working IRAF installation.
  `header_key_read` has now been moved to `warp/fits_utils.py`, so
  `warp.config` no longer imports PyRAF.
  Direct users of `header_key_read` have been moved to `warp.fits_utils`,
  but there are still other modules with PyRAF import-time side effects.
- `testWarpSci*.sh` still runs the real pipeline.
  It may fail if output directories already exist.
  That behavior is unchanged.
- The first lightweight tests were verified with Python 3.7.16, 3.12.11, and
  3.13.5.
  - `python3 -m pytest -q` -> 8 passed
- Current smoke-test coverage is weighted toward `Warp_sci.py`.
  `Warp_calib.py` has not yet received the same execution coverage or
  refactoring attention in this PR.
