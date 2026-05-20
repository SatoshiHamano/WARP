import importlib
import sys


def import_without_pyraf(module_name):
    sys.modules.pop(module_name, None)
    sys.modules.pop("pyraf", None)

    importlib.import_module(module_name)

    assert "pyraf" not in sys.modules


def test_make_rawdata_list_import_does_not_import_pyraf():
    import_without_pyraf("make_rawdata_list")


def test_warp_package_import_does_not_import_pyraf():
    import_without_pyraf("warp")


def test_output_layout_import_does_not_import_pyraf():
    import_without_pyraf("warp.output_layout")


def test_report_import_does_not_import_pyraf():
    import_without_pyraf("warp.report")


def test_tex_source_maker_wrapper_import_does_not_import_pyraf():
    import_without_pyraf("tex_source_maker")


def test_tex_source_maker_wrapper_exports_report_api():
    tex_source_maker = importlib.import_module("tex_source_maker")
    report = importlib.import_module("warp.report")

    assert tex_source_maker.tex_source_make is report.tex_source_make
