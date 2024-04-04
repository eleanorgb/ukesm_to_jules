import iris
import ants
import numpy as np
import jules

PWDANCIL = "/hpc/projects/um1/ancil/atmos/"

JULES_ES_MAP = {
    0: 101,
    1: 102,
    2: 103,
    3: 201,
    4: 202,
    5: 3,
    6: 301,
    7: 302,
    8: 4,
    9: 401,
    10: 402,
    11: 501,
    12: 502,
    13: 6,
    14: 7,
    15: 8,
    16: 9,
}

JULES_ES_MAP_ICE = {
    **JULES_ES_MAP,
    17: 901,
    18: 902,
    19: 903,
    20: 904,
    21: 905,
    22: 906,
    23: 907,
    24: 908,
    25: 909,
    26: 910,
}

# #############################################################
def check_and_regrid(inCube, maskCube):
    """
    if JULES grid differs from landMask grid
            can the JULES latitudes be extended?
            can the JULES data be regrided?
    """
    outCube = None
    match = {"latitude": True, "longitude": True}
    inCube.coord("longitude").coord_system = maskCube.coord("longitude").coord_system
    inCube.coord("latitude").coord_system = maskCube.coord("latitude").coord_system
    inCube.coord("longitude").circular = True
    maskCube.coord("longitude").circular = True
    ants.utils.cube.guess_horizontal_bounds(inCube)
    ants.utils.cube.guess_horizontal_bounds(maskCube)

    # check longitudes
    for coord in ["latitude", "longitude"]:
        if len(inCube.coord(coord).points) != len(maskCube.coord(coord).points):
            match[coord] = False
        elif not (inCube.coord(coord).points == maskCube.coord(coord).points).all():
            match[coord] = False

    if not any(match.values()):
        print("REGRID: Neither latitude nor longitude match")
        outCube = inCube.regrid(maskCube, iris.analysis.Nearest())

    elif match["latitude"] is False and match["longitude"] is True:
        print("Longitude matches but latitude does not")
        # check if jules is on correct grid with reduced no of latitude
        if inCube.coord("latitude").points[0] in maskCube.coord("latitude").points:
            firstIdx = np.where(
                maskCube.coord("latitude").points == inCube.coord("latitude").points[0]
            )[0][0]
            if (
                inCube.coord("latitude").points[1]
                == maskCube.coord("latitude").points[firstIdx + 1]
            ):
                print("JULES on requested output grid, but not all latitudes: PAD")
                outCube = jules.pad_latitude(inCube, maskCube)
        else:
            print("REGRID: Neither latitude nor longitude match")
            outCube = inCube.regrid(maskCube, iris.analysis.Nearest())
    elif match["latitude"] and match["longitude"]:
        print("Both latitude and longitude match")
        outCube = inCube
    else:
        raise ValueError("Check grids of input and masks")

    return outCube


# #############################################################
def choose_mask(soilMask, iceMask, landMask, var_mask_choice=None):
    if var_mask_choice == None:
        mask = landMask
    elif var_mask_choice == "ice":
        mask = iceMask
    elif var_mask_choice == "soil":
        mask = soilMask
    else:
        raise ValueError("No option for mask =", var_mask_choice)
    return mask


# #############################################################
def define_masks(oceanGrid):
    """files used to define each output grid:"""
    if oceanGrid == "n96_ORCA025":
        # ORCA025 settings
        landMaskFile = PWDANCIL + "n96/orca025/land_sea_mask/etop01/v0/qrparm.mask"
        soilMaskFile = PWDANCIL + "n96/orca025/soil_parameters/hwsd_vg/v1/qrparm.soil"
    elif oceanGrid == "n96_ORCA1":
        # ORCA1 settings
        landMaskFile = PWDANCIL + "n96/orca1/land_sea_mask/etop01/v0/qrparm.mask"
        soilMaskFile = PWDANCIL + "n96/orca1/soil_parameters/hwsd_vg/v1/qrparm.soil"
    elif oceanGrid == "n96e_ORCA025":
        # ORCA025 settings
        landMaskFile = PWDANCIL + "n96e/orca025/land_sea_mask/etop01/v2/qrparm.mask"
        soilMaskFile = PWDANCIL + "n96e/orca025/soil_parameters/hwsd_vg//v3/qrparm.soil"
    elif oceanGrid == "n96e_ORCA1":
        # ORCA1 settings
        landMaskFile = PWDANCIL + "n96e/orca1/land_sea_mask/etop01/v0/qrparm.mask"
        soilMaskFile = PWDANCIL + "n96e/orca1/soil_parameters/hwsd_vg/v2/qrparm.soil"
    elif oceanGrid == "LGM":
        # last glacial maximum settings:
        landMaskFile = "/data/users/hadtq/lgm/qrparm.mask"
        soilMaskFile = "/data/users/hadtq/lgm/qrparm.soil"
    else:
        raise ValueError(
            "oceanGrid= " + oceanGrid + ": need landMaskFile and soilMaskFile"
        )
    return (landMaskFile, soilMaskFile)


# #############################################################
def get_masks(landMaskFile, soilMaskFile, soilMaskVariable="m01s00i043"):
    """using sm_sat as standard way of defining soil"""
    landMask = ants.load_cube(landMaskFile)
    soilData = ants.load_cube(
        soilMaskFile, iris.AttributeConstraint(STASH=soilMaskVariable)
    )

    # define icemask
    iceMask = landMask.copy()
    iceMask.data[:] = 0.0
    iceMask.data[soilData.data == 0.0] = 1.0
    iceMask.data[landMask.data == 0.0] = 0.0

    # define soilmask
    soilMask = landMask.copy()
    soilMask.data[:] = 0.0
    soilMask.data[landMask.data * soilData.data > 0.0] = 1.0

    # print sizes of arrays
    print(
        "number of sea, land, and total points",
        len(np.where(landMask.data == 0.0)[0]),
        np.sum(landMask.data),
        len(np.where(landMask.data == 0.0)[0]) + np.sum(landMask.data),
    )
    print(
        "number of soil, ice, and total points",
        np.sum(soilMask.data),
        np.sum(iceMask.data),
        np.sum(soilMask.data) + np.sum(iceMask.data),
    )
    # iris.save(landMask,"/scratch/hadea/landmask.nc")
    # iris.save(iceMask,"/scratch/hadea/icemask.nc")
    # iris.save(soilMask,"/scratch/hadea/soilmask.nc")
    return landMask, iceMask, soilMask
