# WARP
WINERED Automatic Reduction Pipeline

## What is WARP?
The pipeline software to reduce the astronomical spectroscopic data obtained with NIR high-resolution echelle spectrograph, WINERED. WARP is written with Python.

## How to install?

Using git:

`git clone https://github.com/SatoshiHamano/WARP`

WARP can also be installed by downloading zip from GitHub page. Just expand the zip to use it.

## Necessary environment

```
Python 3 (ver 3.6 or later)
Python libraries — numpy, matplotlib, PIL, astropy
PyRAF
```

## How to use?

See WARP_Manual_v?.?.pdf for detail.

## Command line usage

### Science reduction

```sh
python Warp_sci.py LISTFILE [options]
```

`LISTFILE` is a text file listing object/sky frame pairs.  In the current
format, each line is:

```text
OBJECT_FRAME SKY_FRAME [ap=LOW:HIGH] [bg=REGION] [ws=SHIFT]
```

- `ap=LOW:HIGH`: manual aperture limits, used when manual aperture mode is
  enabled.
- `bg=REGION`: background subtraction region, used when background subtraction
  is enabled.
- `ws=SHIFT`: manual wavelength-shift value, used when manual shift mode is
  enabled.

Example:

```text
WINA00036701 WINA00036702 ap=-7:3 bg=-22:-12,8:18 ws=0.
```

Main `Warp_sci.py` options:

| Option | Default | Meaning |
| --- | --- | --- |
| `listfile` | required | Input object/sky list. |
| `-r`, `--rawdatapath` | `../` | Directory containing raw science FITS files. |
| `-v`, `--viewerpath` | `../` | Directory containing slit-viewer FITS files. |
| `-c`, `--calibpath` | `./` | Directory containing WARP calibration products. |
| `-d`, `--destpath` | `./` | Output directory. WARP exits if this directory already exists. |
| `-q`, `--query` | off | Ask selected pipeline settings interactively. Cannot be used with `-p`. |
| `-s`, `--save` | off | Keep more intermediate products instead of moving them to `Trash`. |
| `-p`, `--parameterfile` | none | Read pipeline settings from a parameter file. |
| `-o`, `--oldformat` | off | Read old input-list format used by WARP <= 3.6. |
| `-f`, `--fastMode` | off | Use a fast smoke-test style setting: no cosmic-ray detection, no wavelength-shift measurement/correction, only `fsr1.05`. |
| `-a`, `--autoCalib` | off | Automatically choose a matching calibration directory below `--calibpath`. |
| `--noreport` | off | Skip PDF report generation. |

`-q` and `-p` are mutually exclusive.

### Calibration reduction

```sh
python Warp_calib.py LISTFILE [options]
```

`Warp_calib.py` expects a calibration input file listing the pinhole, flat-on,
flat-off, and comparison frames.  See `TEST/WIDE-calib/input.list` for an
example.

Main `Warp_calib.py` options:

| Option | Default | Meaning |
| --- | --- | --- |
| `listfile` | required | Calibration input list. |
| `-a`, `--aperturereplace` | off | Rebuild/replace aperture-related products. |
| `-t`, `--transformonly` | off | Run only the transformation part. |

## Pipeline parameter file

`Warp_sci.py -p PARAMETER_FILE` reads a text parameter file.  The parser matches
the labels below, so keeping the label text is important.  See
`TEST/WIDE/paramSample.txt` for a working example.

| Parameter label | Values | Default | Effect |
| --- | --- | --- | --- |
| `Apscatter` | `yes` / `no` | `yes` | Subtract scattered light. |
| `Manual Aperture` | `yes` / `no` | `no` | Read aperture limits from `ap=LOW:HIGH` in the input list. |
| `Background Subtraction` | `none`, `average`, `median`, `minimum`, `fit` | `none` | Subtract background spectra. Non-`none` enables background subtraction and uses `bg=REGION` from the input list. |
| `Cosmic Ray Correction` | `yes` / `no` | `yes` | Detect and interpolate cosmic rays. |
| `Cosmic ray threshold sigma` | float | `10.0` | Cosmic-ray detection threshold. |
| `Cosmic ray maximum sigma` | float | `20.0` | Maximum cosmic-ray threshold when variable thresholds are used. |
| `Cosmic ray Var/Ave ratio` | float | `2.0` | Variance/average threshold for cosmic-ray distribution. |
| `Cosmic ray ratio between slit positions` | float | `1.5` | Cosmic-ray count ratio threshold between slit positions. |
| `Cosmic ray fix sigma` | `yes` / `no` | `no` | Use a fixed cosmic-ray sigma threshold. |
| `Extract all orders` | `yes` / `no` | `yes` | Reduce all available echelle orders. |
| `Selected orders` | comma-separated order numbers, or `no` | `[]` | Orders to reduce when not reducing all orders. |
| `Set cut range` | comma-separated floats, or `no` | `1.05, 1.3` | Free spectral range factors used for final 1D spectra. |
| `Sky Emission` | `yes` / `no` | `no` | Extract spectra from sky frames. |
| `Measure Shift` | `yes` / `no` | `yes` | Measure wavelength shifts among object frames. |
| `Correct Shift` | `yes` / `no` | `yes` | Apply measured or manual wavelength shifts. |
| `Manual Shift` | `yes` / `no` | `no` | Use `ws=SHIFT` values from the input list. |
| `CUTRANSFORM flux` | `yes` / `no` | `no` | Passed to IRAF `transform`/`cutransform` as the flux-conservation option. |
| `Extract 2d spectrum` | `yes` / `no` | `no` | Also extract 2D wavelength/spatial spectra. |

These settings affect the numerical science products.  Reference-output
comparisons should record the parameter file, input list, calibration directory,
and command-line options used to generate the reference data.

## Versioning and releases

WARP does not have a separate production deployment environment.  The `main`
branch should therefore be treated as the latest released version that users may
pull and run.

Release policy:

- Do not push directly to `main`; use pull requests.
- Run the relevant tests before merging changes into `main`.
- Small refactoring or documentation pull requests do not need their own version
  bump.
- Bump the version when a reviewed set of changes is ready to become the new
  user-facing release.
- Mark release commits with annotated Git tags such as `v3.9.0`.
- Use tags, not long-lived version branches, as the normal way to refer to old
  released versions.

Example release tag:

```sh
git tag -a v3.9.0 -m "WARP 3.9.0"
git push origin v3.9.0
```

Users can return to a specific release with:

```sh
git checkout v3.9.0
```

## Tests

Run the lightweight tests with:

```sh
python3 -m pytest -q
```

These tests cover input-list and parameter-file parsing and do not require a
working IRAF/PyRAF installation.

Run the WARP science-pipeline smoke tests with:

```sh
./testWarpSci.sh
./testWarpSciFull.sh
```

Set `PYTHON` to test a specific Python executable:

```sh
PYTHON=python3.11 ./testWarpSci.sh
```
