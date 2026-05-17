from pathlib import Path
import sys

import pytest

from warp.config import config


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_config_import_does_not_import_pyraf():
    assert "pyraf" not in sys.modules


def test_wide_input_list_is_parsed_with_optional_fields():
    conf = config()

    conf.inputDataList(REPO_ROOT / "TEST/WIDE/4_Ari_list.txt")

    assert conf.objectlist == ["WINA00036701", "WINA00036703"]
    assert conf.skylist == ["WINA00036702", "WINA00036702"]
    assert conf.lowlim_input == [-7.0, -7.0]
    assert conf.upplim_input == [3.0, 3.0]
    assert conf.skysub_region == ["-22:-12,8:18", "-22:-12,8:18"]
    assert conf.waveshift_man == [0.0, 0.0]
    assert conf.objnum == 2
    assert conf.imnum == 3
    assert conf.imagelist.tolist() == ["WINA00036701", "WINA00036702", "WINA00036703"]


def test_long_wide_input_list_is_parsed_without_frame_limit_regression():
    conf = config()

    conf.inputDataList(REPO_ROOT / "TEST/WIDE/4_Ari_list_long.txt")

    assert conf.objnum == 100
    assert conf.imnum == 3
    assert len(conf.lowlim_input) == 100
    assert len(conf.upplim_input) == 100
    assert len(conf.skysub_region) == 100
    assert len(conf.waveshift_man) == 100


@pytest.mark.parametrize(
    "list_path, expected_objects, expected_skies, expected_imagelist",
    [
        (
            "TEST/HIRES-J/21_Peg_list.txt",
            ["WINA00027411", "WINA00027413"],
            ["WINA00027412", "WINA00027412"],
            ["WINA00027411", "WINA00027412", "WINA00027413"],
        ),
        (
            "TEST/HIRES-Y/HD_163336_list.txt",
            ["WINA00027326", "WINA00027327"],
            ["WINA00027327", "WINA00027326"],
            ["WINA00027326", "WINA00027327"],
        ),
    ],
)
def test_hires_input_lists_are_parsed(list_path, expected_objects, expected_skies, expected_imagelist):
    conf = config()

    conf.inputDataList(REPO_ROOT / list_path)

    assert conf.objectlist == expected_objects
    assert conf.skylist == expected_skies
    assert conf.objnum == len(expected_objects)
    assert conf.imagelist.tolist() == expected_imagelist


def test_input_list_requires_apertures_when_manual_aperture_is_enabled(tmp_path):
    input_list = tmp_path / "input.list"
    input_list.write_text("OBJ SKY\n")
    conf = config()
    conf.flag_manual_aperture = True

    with pytest.raises(SystemExit):
        conf.inputDataList(input_list)


def test_parameter_sample_updates_reduction_flags(capsys):
    conf = config()

    conf.readParamFile(REPO_ROOT / "TEST/WIDE/paramSample.txt")

    assert conf.flag_apscatter is True
    assert conf.flag_manual_aperture is True
    assert conf.flag_skysub is True
    assert conf.skysub_mode == "average"
    assert conf.flag_bpmask is True
    assert conf.CRthreshold == 10.0
    assert conf.CRmaxsigma == 20.0
    assert conf.CRvaratio == 2.0
    assert conf.CRslitposratio == 1.5
    assert conf.CRfixsigma is False
    assert conf.cutrange_list == [1.05, 1.3]
    assert conf.flag_skyemission is False
    assert conf.flag_wsmeasure is True
    assert conf.flag_wscorrect is True
    assert conf.flag_wsmanual is False
    assert conf.fluxinput == "no"

    capsys.readouterr()
