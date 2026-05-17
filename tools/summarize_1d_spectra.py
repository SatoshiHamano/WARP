#!/usr/bin/env python3
"""Summarize WARP 1D spectrum FITS outputs for regression tests."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from astropy.io import fits


HEADER_KEYS = (
    "NAXIS",
    "NAXIS1",
    "CRVAL1",
    "CRPIX1",
    "CDELT1",
    "CTYPE1",
    "CUNIT1",
    "BUNIT",
    "AIRORVAC",
)


def json_safe_float(value):
    if value is None:
        return None
    value = float(value)
    if not math.isfinite(value):
        return None
    return value


def summarize_array(values):
    flat = np.asarray(values, dtype=np.float64).reshape(-1)
    finite = flat[np.isfinite(flat)]
    summary = {
        "size": int(flat.size),
        "finite_count": int(finite.size),
        "nan_count": int(np.isnan(flat).sum()),
        "inf_count": int(np.isinf(flat).sum()),
    }
    if finite.size == 0:
        summary.update(
            {
                "min": None,
                "max": None,
                "mean": None,
                "median": None,
                "std": None,
                "sum": None,
                "percentiles": {},
                "samples": [],
                "window_medians": [],
            }
        )
        return summary

    summary.update(
        {
            "min": json_safe_float(np.min(finite)),
            "max": json_safe_float(np.max(finite)),
            "mean": json_safe_float(np.mean(finite)),
            "median": json_safe_float(np.median(finite)),
            "std": json_safe_float(np.std(finite)),
            "sum": json_safe_float(np.sum(finite)),
            "percentiles": {
                str(percentile): json_safe_float(np.percentile(finite, percentile))
                for percentile in (1, 5, 25, 50, 75, 95, 99)
            },
            "samples": sample_points(flat),
            "window_medians": window_medians(flat),
        }
    )
    return summary


def sample_indices(size):
    if size <= 0:
        return []
    indices = {0, size // 4, size // 2, (3 * size) // 4, size - 1}
    return sorted(indices)


def sample_points(values):
    samples = []
    for index in sample_indices(values.size):
        samples.append({"index": int(index), "value": json_safe_float(values[index])})
    return samples


def window_medians(values, half_width=2):
    windows = []
    for index in sample_indices(values.size):
        start = max(0, index - half_width)
        stop = min(values.size, index + half_width + 1)
        window = values[start:stop]
        finite = window[np.isfinite(window)]
        median = None if finite.size == 0 else json_safe_float(np.median(finite))
        windows.append(
            {
                "center_index": int(index),
                "start_index": int(start),
                "stop_index": int(stop),
                "median": median,
            }
        )
    return windows


def wavelength_summary(header, size):
    crval = header.get("CRVAL1")
    crpix = header.get("CRPIX1")
    cdelt = header.get("CDELT1")
    if size <= 0 or crval is None or crpix is None or cdelt is None:
        return {
            "available": False,
            "start": None,
            "end": None,
            "step": None,
            "monotonic_increasing": None,
        }

    indices = np.arange(size, dtype=np.float64) + 1.0
    wave = float(crval) + (indices - float(crpix)) * float(cdelt)
    diffs = np.diff(wave)
    return {
        "available": True,
        "start": json_safe_float(wave[0]),
        "end": json_safe_float(wave[-1]),
        "min": json_safe_float(np.min(wave)),
        "max": json_safe_float(np.max(wave)),
        "step": json_safe_float(float(cdelt)),
        "step_median": json_safe_float(np.median(diffs)) if diffs.size else None,
        "step_std": json_safe_float(np.std(diffs)) if diffs.size else None,
        "monotonic_increasing": bool(np.all(diffs > 0)) if diffs.size else True,
        "samples": sample_points(wave),
    }


def summarize_fits(path, root):
    with fits.open(path, memmap=False) as hdul:
        hdu = hdul[0]
        data = np.asarray(hdu.data)
        if data.ndim != 1:
            raise ValueError(f"{path} is not a 1D spectrum: shape={data.shape}")
        header = {key: hdu.header[key] for key in HEADER_KEYS if key in hdu.header}
        return {
            "path": path.relative_to(root).as_posix(),
            "shape": list(data.shape),
            "dtype": str(data.dtype),
            "header": header,
            "wavelength": wavelength_summary(hdu.header, data.size),
            "flux": summarize_array(data),
        }


def find_1d_spectra(root):
    patterns = (
        "*_sum/AIR_norm/fsr*/*.fits",
        "*_sum/VAC_norm/fsr*/*.fits",
        "*_sum/AIR_flux/fsr*/*.fits",
        "*_sum/VAC_flux/fsr*/*.fits",
    )
    paths = []
    for pattern in patterns:
        paths.extend(root.glob(pattern))
    return sorted(set(path for path in paths if path.is_file()))


def summarize_tree(root):
    root = Path(root).resolve()
    files = find_1d_spectra(root)
    if not files:
        raise FileNotFoundError(f"No 1D spectrum FITS files found below {root}")
    return {
        "schema_version": 1,
        "root_name": root.name,
        "file_count": len(files),
        "files": [summarize_fits(path, root) for path in files],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path, help="WARP output tree containing *_sum directories")
    parser.add_argument("--output", "-o", type=Path, help="Write summary JSON to this path")
    args = parser.parse_args(argv)

    summary = summarize_tree(args.root)
    text = json.dumps(summary, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n")
    else:
        print(text)


if __name__ == "__main__":
    main()
