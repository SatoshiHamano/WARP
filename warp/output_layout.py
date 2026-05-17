"""Directory names used by WARP reduction outputs."""

ONEDSPEC_DIRNAMES = ["AIR_flux", "AIR_norm", "AIR_cont", "VAC_flux", "VAC_norm", "VAC_cont"]
TWODSPEC_DIRNAMES = ["AIR", "VAC"]
SKYEMISSION_DIRNAMES = ["AIR", "VAC"]

INTERMEDIATE_OBJ_DIRNAMES = [
    "1-OBJ_sky_subs",
    "2-OBJ_scatter_subs",
    "3-OBJ_flat",
    "4-OBJ_mask",
    "5-OBJ_cut",
    "6-OBJ_transform",
    "7-OBJ_mask2",
]

INTERMEDIATE_OBJ_1DSPEC_DIRNAMES = [
    "1-OBJ-1DSPEC_extract",
    "2-OBJ-1DSPEC_truncate",
    "3-OBJ-1DSPEC_shift",
    "4-OBJ-1DSPEC_dispcor",
]

INTERMEDIATE_OBJ_2DSPEC_DIRNAMES = [
    "1-OBJ-2DSPEC_extract",
    "2-OBJ-2DSPEC_truncate",
    "3-OBJ-2DSPEC_shift",
]

INTERMEDIATE_SKY_DIRNAMES = [
    "1-SKY_flat",
    "2-SKY_mask",
    "3-SKY_cut",
    "4-SKY_transform",
    "5-SKY_extract",
    "6-SKY_truncate",
    "7-SKY_dispcor",
]
