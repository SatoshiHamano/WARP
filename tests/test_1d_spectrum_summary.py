import json
import os
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from tools.summarize_1d_spectra import build_metadata, summarize_fits, summarize_tree


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


def test_summarize_tree_records_run_metadata(tmp_path):
    input_list = tmp_path / "input.list"
    input_list.write_text("OBJ SKY ap=-1:1 bg=-5:5 ws=0.0\n")
    spectrum = tmp_path / "star_sum/AIR_flux/fsr1.05/star_m52_fsr1.05_AIR.fits"
    write_spectrum(spectrum, [10.0, 11.0, 12.0])

    metadata = build_metadata(
        case_name="unit_case",
        command=["Warp_sci.py input.list -f"],
        input_list=input_list,
        parameter_file=None,
        calibration_path="calib",
        rawdata_path="raw",
        viewer_path="viewer",
        warp_version="test-version",
        pyraf_version="test-pyraf",
        iraf_path="/iraf",
        irafarch="macos64",
        extra_metadata={"mode": "fast"},
    )
    summary = summarize_tree(tmp_path, metadata=metadata)

    assert summary["schema_version"] == 2
    assert summary["metadata"]["case_name"] == "unit_case"
    assert summary["metadata"]["command"] == ["Warp_sci.py input.list -f"]
    assert summary["metadata"]["input_list"]["text"] == "OBJ SKY ap=-1:1 bg=-5:5 ws=0.0\n"
    assert summary["metadata"]["parameter_file"] is None
    assert summary["metadata"]["calibration_path"] == "calib"
    assert summary["metadata"]["warp_version"] == "test-version"
    assert summary["metadata"]["pyraf_version"] == "test-pyraf"
    assert summary["metadata"]["iraf"] == "/iraf"
    assert summary["metadata"]["irafarch"] == "macos64"
    assert summary["metadata"]["extra"] == {"mode": "fast"}


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
def test_1d_summary_matches_reference():
    output_root = os.environ.get("WARP_1D_OUTPUT_ROOT")
    if not output_root:
        pytest.skip("Set WARP_1D_OUTPUT_ROOT to compare a WARP output tree")

    reference_summary = Path(os.environ.get("WARP_1D_REFERENCE_SUMMARY", REFERENCE_SUMMARY))
    expected = json.loads(reference_summary.read_text())
    actual = summarize_tree(Path(output_root), metadata=expected["metadata"])
    actual["root_name"] = expected["root_name"]

    assert_close(actual, expected)
