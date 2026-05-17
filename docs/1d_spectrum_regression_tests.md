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
  --output tests/reference/wide_4_ari_1d_summary.json
```

For a freshly generated temporary output from `testWarpSci.sh`, keep the output
directory first:

```sh
KEEP_WARP_TEST_OUTPUT=1 ./testWarpSci.sh
python tools/summarize_1d_spectra.py /path/to/4_Ari_WIDE_test \
  --output /tmp/wide_4_ari_1d_summary.json
```

## Compare output against the committed reference

The regression comparison is opt-in because it needs a generated WARP output
tree:

```sh
WARP_1D_OUTPUT_ROOT=TEST/4_Ari_WIDE_test python -m pytest -q -m regression
```

Normal `pytest` runs still execute without generated WARP outputs; the
reference comparison is skipped unless `WARP_1D_OUTPUT_ROOT` is set.
