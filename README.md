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
