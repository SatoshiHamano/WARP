import json
import os
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from tools.summarize_1d_spectra import summarize_fits, summarize_tree


REPO_ROOT = Path(__file__).resolve().parents[1]
REFERENCE_SUMMARY = REPO_ROOT / "tests/reference/wide_4_ari_1d_summary.json"


def write_spectrum(path, data, crval=1000.0, cdelt=0.5, air_or_vac="air"):
    path.parent.mkdir(parents=True, exist_ok=True)
    hdu = fits.PrimaryHDU(np.asarray(data, dtype=np.float32))
    hdu.header["CRVAL1"] = crval
    hdu.header["CRPIX1"] = 1.0
    hdu.header["CDELT1"] = cdelt
    hdu.header["CTYPE1"] = "LINEAR"
    hdu.header["AIRORVAC"] = air_or_vac
    hdu.writeto(path)


def test_summarize_fits_records_wavelength_and_flux_statistics(tmp_path):
    spectrum = tmp_path / "star_sum/AIR_norm/fsr1.05/star_m52_fsr1.05_AIR_norm.fits"
    write_spectrum(spectrum, [1.0, 2.0, np.nan, 4.0, 5.0])

    summary = summarize_fits(spectrum, tmp_path)

    assert summary["path"] == "star_sum/AIR_norm/fsr1.05/star_m52_fsr1.05_AIR_norm.fits"
    assert summary["shape"] == [5]
    assert summary["header"]["AIRORVAC"] == "air"
    assert summary["wavelength"]["start"] == 1000.0
    assert summary["wavelength"]["end"] == 1002.0
    assert summary["wavelength"]["step"] == 0.5
    assert summary["wavelength"]["monotonic_increasing"] is True
    assert summary["flux"]["finite_count"] == 4
    assert summary["flux"]["nan_count"] == 1
    assert summary["flux"]["median"] == 3.0
    assert summary["flux"]["samples"][0] == {"index": 0, "value": 1.0}


def assert_close(actual, expected, path="summary"):
    if isinstance(expected, dict):
        assert set(actual) == set(expected), path
        for key in expected:
            assert_close(actual[key], expected[key], f"{path}.{key}")
    elif isinstance(expected, list):
        assert len(actual) == len(expected), path
        for index, expected_item in enumerate(expected):
            assert_close(actual[index], expected_item, f"{path}[{index}]")
    elif isinstance(expected, float):
        assert actual == pytest.approx(expected, rel=1e-6, abs=1e-6), path
    else:
        assert actual == expected, path


@pytest.mark.regression
def test_wide_4_ari_1d_summary_matches_reference():
    output_root = os.environ.get("WARP_1D_OUTPUT_ROOT")
    if not output_root:
        pytest.skip("Set WARP_1D_OUTPUT_ROOT to compare a WARP output tree")

    actual = summarize_tree(Path(output_root))
    expected = json.loads(REFERENCE_SUMMARY.read_text())

    assert_close(actual, expected)
