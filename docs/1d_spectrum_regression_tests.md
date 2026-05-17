# 1D spectrum regression tests

The most important scientific product of `Warp_sci.py` is the extracted 1D
spectrum.  The lightweight regression helper added here records numerical
summaries of those FITS files instead of committing the generated FITS products
themselves.

The summary currently covers the final per-order and combined files below each
`*_sum` directory:

- `AIR_norm/fsr*/*.fits`
- `VAC_norm/fsr*/*.fits`
- `AIR_flux/fsr*/*.fits`
- `VAC_flux/fsr*/*.fits`

For each 1D FITS file, the summary records the wavelength calibration header,
derived wavelength range and step, flux statistics, percentiles, and a few
sample/window values across the spectrum.  This keeps the reference file small
while still catching changes in wavelength coverage, sampling, flux scale, and
local spectral shape.

## Generate a summary

```sh
python tools/summarize_1d_spectra.py TEST/4_Ari_WIDE_test \
  --case-name wide_4_ari_fast \
  --input-list TEST/WIDE/4_Ari_list.txt \
  --calibration-path TEST/WIDE/WINERED_calibration_LCO22b_WIDE100_20220914_v2 \
  --rawdata-path TEST/WIDE \
  --viewer-path TEST/WIDE \
  --warp-version 3.8.15 \
  --metadata reduction_mode=fast \
  --output tests/reference/wide_4_ari_1d_summary.json
```

For a freshly generated temporary output from `testWarpSci.sh`, keep the output
directory first:

```sh
KEEP_WARP_TEST_OUTPUT=1 ./testWarpSci.sh
python tools/summarize_1d_spectra.py /path/to/4_Ari_WIDE_test \
  --case-name wide_4_ari_fast \
  --input-list TEST/WIDE/4_Ari_list.txt \
  --metadata reduction_mode=fast \
  --output /tmp/wide_4_ari_1d_summary.json
```

The summary JSON uses schema version 2.  It includes a `metadata` block so the
numerical reference is tied to the run conditions that produced it.  At minimum,
record the case name, input list, calibration directory, parameter file when one
is used, and relevant command-line options.  Use repeated `--metadata KEY=VALUE`
arguments for additional case-specific notes.

## Compare output against the committed reference

The regression comparison is opt-in because it needs a generated WARP output
tree:

```sh
WARP_1D_OUTPUT_ROOT=TEST/4_Ari_WIDE_test python -m pytest -q -m regression
```

Normal `pytest` runs still execute without generated WARP outputs; the
reference comparison is skipped unless `WARP_1D_OUTPUT_ROOT` is set.
