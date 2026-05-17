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
- Verified lightweight tests on newer Python versions.
  - Python 3.7.16 -> 7 passed
  - Python 3.12.11 -> 7 passed
  - Python 3.13.5 -> 7 passed
- Verified one full WIDE fast smoke run with the current Astroconda/PyRAF
  environment.
  - Explicit environment:
    `IRAFARCH=macos64 iraf=/Users/hamano/iraf/iraf-2.17.1/`
  - Command shape:
    `python Warp_sci.py ./TEST/WIDE/4_Ari_list.txt ... -d /private/tmp/warp_codex_4_Ari_WIDE_test -f`
  - Result: exited successfully and printed `=== Finished. ===`

## Related Commits

- `83d9749 Add gitignore for local artifacts`
- `6623e6f Remove tracked local artifacts`
- `60ce277 Add lightweight input parsing tests`
- `a52abd8 Document refactor plan`
- `1f9e977 Document PyRAF environment caveat`

## Next Steps

1. Commit and push the current PyRAF import-isolation change.
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
- `warp.config` previously imported `warp.Spec2Dtools`, and
  `warp.Spec2Dtools` imports PyRAF at module import time.
  This made simple config tests depend on a working IRAF installation.
  `header_key_read` has now been moved to `warp/fits_utils.py`, so
  `warp.config` no longer imports PyRAF.
  There may still be other modules with PyRAF import-time side effects.
- `testWarpSci*.sh` still runs the real pipeline.
  It may fail if output directories already exist.
  That behavior is unchanged.
- The first lightweight tests were verified with Python 3.7.16, 3.12.11, and
  3.13.5.
  - `python3 -m pytest -q` -> 7 passed
