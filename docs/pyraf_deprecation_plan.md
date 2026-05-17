# PyRAF Deprecation Plan

This document lists the current computational PyRAF/IRAF dependencies in WARP
and proposes a staged plan to remove them without changing scientific results
accidentally.

## Goal

The long-term goal is to make WARP runnable on modern Python without PyRAF.
This requires replacing IRAF task calls with Python implementations or small
compatibility wrappers that can be tested numerically.

This is larger and riskier than import cleanup. Any replacement of a numerical
IRAF task should be protected by reference-output comparisons.

## Current Ecosystem Check

Astropy should be the first-choice replacement foundation for most low-level
WARP I/O and metadata operations.

Current evidence:

- `astropy.io.fits` is the actively documented core Python interface for FITS
  files, headers, and array data:
  <https://docs.astropy.org/en/stable/io/fits/>
- `astropy.wcs` is the core Astropy interface for FITS WCS transformations:
  <https://docs.astropy.org/en/stable/wcs/>
- `ccdproc` is an Astropy-affiliated package for CCD image reduction:
  <https://ccdproc.readthedocs.io/en/stable/>
- `specutils` is an Astropy-affiliated spectroscopy package, but its own
  project description says it is not intended to cover all spectroscopic
  analysis or reduction needs:
  <https://github.com/astropy/specutils>
- `gwcs` is an Astropy-affiliated generalized WCS package and may be useful
  later for non-trivial detector-to-world transformations:
  <https://gwcs.readthedocs.io/en/stable/>
- STScI documentation recommends Python 3 tools and Astropy for modern
  analysis, and notes that STScI stopped supporting IRAF/PyRAF in 2019:
  <https://hst-docs.stsci.edu/hstdhb/4-hst-data-analysis/4-1-analysis-options-for-hst-data>
- PyRAF is not completely abandoned: the IRAF Community now maintains PyRAF,
  and recent community releases mention Python 3.13 compatibility. However,
  PyRAF remains a command-language bridge to IRAF rather than a native Python
  scientific stack:
  <https://iraf-community.github.io/pyraf.html>
  <https://zenodo.org/records/17341450>

Conclusion:

- Astropy remains the right default migration target for FITS I/O, headers,
  WCS metadata, and simple numerical array operations.
- Astropy-affiliated packages such as `ccdproc`, `specutils`, and `gwcs`
  should be evaluated case by case, not adopted wholesale.
- PyRAF can remain a compatibility backend while WARP builds numerical
  baselines and replaces task families incrementally.
- The motivation should be long-term maintainability, reproducibility, and
  modern Python compatibility, not the claim that PyRAF is currently unusable.

## Current Dependency Inventory

Files that import PyRAF or IRAF directly:

- `ArrayMonitor.py`
- `Warp_calib.py`
- `Warp_sci.py`
- `warp/ECtoID.py`
- `warp/PyContinuum.py`
- `warp/SNratio_estimate.py`
- `warp/Spec1Dtools.py`
- `warp/Spec2Dtools.py`
- `warp/angle_measure.py`
- `warp/apscatter.py`
- `warp/auto_ecidentify.py`
- `warp/badpixmask.py`
- `warp/ccwaveshift.py`
- `warp/centersearch_fortrans.py`
- `warp/cutransform.py`
- `warp/vac2air_spec.py`

Direct IRAF task usage found in the codebase:

- FITS header editing:
  - `iraf.hedit`
  - `iraf.images.imutil.hedit`
- FITS image arithmetic and copying:
  - `iraf.imarith`
  - `iraf.imcopy`
  - `iraf.imreplace`
- 1D spectral operations:
  - `iraf.scopy`
  - `iraf.scombine`
  - `iraf.sarith`
  - `iraf.specshift`
  - `iraf.continuum`
- 2D spectral transformations:
  - `iraf.transform`
- Aperture and mask operations:
  - `iraf.apmask`
  - `iraf.apnormalize`
  - `iraf.fixpix`
- Wavelength calibration:
  - `iraf.ecidentify`
- IRAF package initialization:
  - `iraf.noao`
  - `iraf.onedspec`
  - `iraf.imred`
  - `iraf.echelle`
  - `iraf.twodspec`
  - `iraf.longslit`
  - `iraf.images`
  - `iraf.imutil`
  - `iraf.imgeom`

Some files import PyRAF but do not currently show direct `iraf.*` task calls:

- `warp/angle_measure.py`
- `warp/auto_ecidentify.py`

These should be checked first. If the imports are unused, remove them in a
small import-only patch.

## Replacement Strategy By Task Family

### 1. FITS Header Editing

Current use:

- `iraf.hedit`
- `iraf.images.imutil.hedit`

Likely replacement:

- `astropy.io.fits` opened in update mode
- a helper such as `set_header_value(fits_path, key, value, add=True)`

Risk level: low to medium.

Notes:

- This is one of the best first targets.
- Need to preserve the old behavior for extension names, missing keys, string
  formatting, and files with or without `.fits` suffixes.

### 2. FITS Image Copying, Slicing, Arithmetic, And Pixel Replacement

Current use:

- `iraf.imcopy`
- `iraf.imarith`
- `iraf.imreplace`

Likely replacement:

- `astropy.io.fits` for reading/writing FITS
- `numpy` for arithmetic and masks
- evaluate `ccdproc` only where WARP needs CCD-style reduction primitives
  with uncertainty or mask propagation
- explicit IRAF-section parser for expressions such as `[1:406,*]`

Risk level: medium.

Notes:

- Simple full-image arithmetic can be replaced early.
- IRAF image-section syntax is a hidden compatibility risk.
- Replacements need tests for header propagation and data type behavior.

### 3. 1D Spectral Copying, Combining, And Arithmetic

Current use:

- `iraf.scopy`
- `iraf.scombine`
- `iraf.sarith`
- `iraf.specshift`

Likely replacement:

- `astropy.io.fits`
- `astropy.wcs` for FITS WCS interpretation where applicable
- `numpy`
- possibly `specutils` for spectrum-aware operations, if its behavior can be
  pinned and tested

Risk level: medium to high.

Notes:

- `scombine` and `scopy` may encode IRAF-specific WCS and multispec behavior.
- `specutils` is useful as a representation and analysis layer, but should not
  be assumed to replace IRAF reduction tasks directly.
- These should be migrated only after reference summaries exist.
- Start with narrow wrappers and compare output headers and sampled numerical
  values before changing call sites.

### 4. Continuum Fitting

Current use:

- `iraf.continuum`

Likely replacement:

- `numpy.polynomial`
- `scipy.interpolate`
- `astropy.modeling`
- custom rejection loop matching the current IRAF parameters

Risk level: high.

Notes:

- This directly affects normalized spectra.
- Treat as a scientific-behavior change unless numerical equivalence is shown.
- Implement only after the regression summary framework is in place.

### 5. 2D Transformations

Current use:

- `iraf.transform`

Likely replacement:

- `scipy.ndimage.map_coordinates`
- `astropy.modeling` or existing transformation database parsing
- evaluate `gwcs` if the transform can be represented cleanly as a reusable
  detector-to-world model
- custom reader for IRAF database transform files

Risk level: high.

Notes:

- This is a core calibration/reduction step.
- It probably depends on IRAF database files produced by calibration tasks.
- Do not replace until calibration reference outputs and transformation
  metadata are well understood.

### 6. Aperture Masking, Normalization, And Bad Pixel Fixing

Current use:

- `iraf.apmask`
- `iraf.apnormalize`
- `iraf.fixpix`

Likely replacement:

- custom mask generation based on existing aperture database parsing
- `numpy` / `scipy.ndimage`
- interpolation along selected axes for `fixpix`

Risk level: high.

Notes:

- These tasks are tied to aperture database semantics.
- Replace only after simpler FITS utility replacements are stable.

### 7. Wavelength Calibration

Current use:

- `iraf.ecidentify`

Likely replacement:

- keep current database files as inputs where possible
- later replace line identification with explicit Python fitting code

Risk level: very high.

Notes:

- This is likely the hardest PyRAF dependency to remove.
- It may require an interactive or semi-interactive alternative.
- For the first deprecation stage, keep this behind an adapter rather than
  replacing it.

## Proposed Staged Plan

### Stage 0: Stabilize Tests And Baselines

- Add numerical regression summaries for the WIDE fast smoke run.
- Record selected output FITS headers, shapes, and numerical fingerprints.
- Keep regression tests opt-in because they require PyRAF/IRAF and test data.

Exit criteria:

- A known-good PyRAF run can generate a compact reference summary.
- A new run can be compared against that summary.

### Stage 1: Remove Unused PyRAF Imports

- Check `warp/angle_measure.py` and `warp/auto_ecidentify.py`.
- Remove PyRAF imports if unused.
- Add import-boundary tests where the module can reasonably import without
  PyRAF.

Exit criteria:

- No behavior changes.
- Lightweight tests still pass on Python 3.7, 3.12, and 3.13.

### Stage 2: Introduce An IRAF Adapter Layer

- Add a small module such as `warp/iraf_adapter.py`.
- Move direct `iraf.*` calls behind named Python functions.
- Keep the implementation PyRAF-backed at first.
- Do not change numerical behavior in this stage.

Example wrappers:

- `hedit(path, key, value, add=True)`
- `imarith(left, operator, right, output)`
- `imcopy(input_spec, output)`
- `scopy(input_spec, output, **kwargs)`
- `scombine(inputs, output, **kwargs)`
- `continuum(input, output, **kwargs)`
- `transform(input, output, reference, **kwargs)`

Exit criteria:

- Direct PyRAF usage is concentrated in one adapter module.
- High-level pipeline files no longer import PyRAF directly.
- Existing smoke test output remains numerically equivalent.

### Stage 3: Replace Low-Risk FITS Utilities

- Replace adapter implementations for `hedit`, simple `imcopy`, simple
  `imarith`, and `imreplace` with Astropy/Numpy code.
- Keep the PyRAF implementation available as a fallback during migration.

Exit criteria:

- Numerical summaries match the PyRAF baseline within agreed tolerances.
- Header changes are explicitly reviewed.

### Stage 4: Replace Medium-Risk Spectral Operations

- Replace limited, well-understood usages of `scopy`, `sarith`, `scombine`,
  and `specshift`.
- Avoid broad generality; implement only the modes WARP actually uses.

Exit criteria:

- WIDE fast regression summary matches.
- Representative HIRES cases have at least smoke-test coverage.

### Stage 5: Replace High-Risk Scientific Tasks

- Replace `continuum`, `transform`, `apmask`, `apnormalize`, `fixpix`, and
  eventually `ecidentify`.
- Treat each task family as a separate PR or set of PRs.

Exit criteria:

- Reference-output comparison exists before replacement.
- Scientific owner reviews numerical differences.
- Differences are documented, either as equivalent within tolerance or as an
  intentional algorithmic change.

## Suggested PR Boundaries

Good PR boundaries:

- remove unused PyRAF imports
- add `iraf_adapter.py` without changing behavior
- move one task family behind the adapter
- replace one low-risk adapter implementation
- add or expand numerical regression summaries

Avoid in one PR:

- replacing `transform` and `continuum` together
- changing calibration logic while also changing file layout
- changing PyRAF behavior without a reference-output comparison

## Open Questions

- Which output products are scientifically essential for equivalence checks?
- How much numerical drift is acceptable for normalized spectra?
- Should the no-PyRAF backend aim for exact IRAF compatibility or documented
  scientific equivalence?
- Are there users who depend on intermediate files, not only final products?
- Should interactive tasks such as `ecidentify` remain supported, or should
  WARP move toward fully scripted calibration?
