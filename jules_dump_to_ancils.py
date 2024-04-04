"""
convert JULES dump into UM ancil
make generic enough for most JULES variables
"""

import jules
import iris
import ants
import ants.analysis
import numpy as np
import mule
import iris.fileformats
import iris.analysis
import ants.utils.cube
from define_jules_to_umdump_masks import (
    define_masks,
    get_masks,
    choose_mask,
    JULES_ES_MAP,
    JULES_ES_MAP_ICE,
    check_and_regrid,
)

import um_to_jules_stash_dict

stashall = [
    #     "m01s00i449",  # "frac_agr_prev"
    #     "m01s00i459",  # "frac_past_prev"
    #     "m01s00i287",  # "wood_prod_fast"
    #     "m01s00i288",  # "wood_prod_med",
    #     "m01s00i289",  # "wood_prod_slow"
    #     "m01s00i442",  # "ns_dpm"
    #     "m01s00i443",  # "ns_rpm"
    #     "m01s00i444",  # "ns_bio"
    #     "m01s00i445",  # "ns_hum"
    #     "m01s00i446",  # "n_inorg"
    #     "m01s00i466",  # "cs_dpm"
    #     "m01s00i467",  # "cs_rpm"
    #     "m01s00i468",  # "cs_bio"
    #     "m01s00i469",  # "cs_hum"
    #     "m01s00i281",  # "sthzw"
    #     "m01s00i278",  # "zw"
    #    # "m01s00i214",  # "sthu" have sum
    #    # "m01s00i215",  # "sthf" have sum
    #     "m01s00i020",  # "t_soil"
    #     "m01s00i240",  # "snow_tile"
    #     "m01s00i216",  # "frac"
    #     "m01s00i217",  # "lai"
    #     "m01s00i218",  # "canht"
    #     "m01s00i377",  # "rho_snow"
    #     "m01s00i376",  # "snow_depth"
    #     "m01s00i242",  # "snow_grnd"
    #     "m01s00i380",  # "nsnow"
    "m01s00i381",  # "snow_ds"
    "m01s00i382",  # "snow_ice"
    "m01s00i383",  # "snow_liq"
    "m01s00i384",  # "tsnow"
]


stashsel = ["m01s00i381"]  # , "m01s00i382", "m01s00i383", "m01s00i384", "m01s00i386"]
# ['m01s00i216']#'m01s00i384']#,'m01s00i384','m01s00i217','m01s00i218','m01s00i466','m01s00i467','m01s00i468','m01s00i469',
#'m01s00i458','m01s00i448',
#'m01s00i442','m01s00i443','m01s00i444','m01s00i445','m01s00i446']


# output grids
oceanGrid = "n96e_ORCA1"  # "LGM"

DICT_STASH = um_to_jules_stash_dict.get_jules_stash_dict()


# #############################################################
def process_frac(cube, maskCube, iceMask, landMask):
    """sort ice and bare soil points and ensure sum frac is 1"""

    bareSoilTileID = np.where(cube.coord("pseudo_level").points == 8)[0][0]
    iceTileID = [key for key, value in type_info.items() if str(value).startswith("9")]
    nonIceTileID = [
        key for key, value in type_info.items() if not str(value).startswith("9")
    ]
    # print(bareSoilTileID, iceTileID)

    cube = check_and_regrid(cube, maskCube)
    # deal with missing data points (not sure this is right)
    filler = ants.analysis.FillMissingPoints(cube)
    filler(cube)
    cube = cube.copy()  # odd requirement

    # where icemask is >0.5 set values in cube.data 0.0
    iceLoc = np.where(iceMask.data > 0.5)
    for itype in nonIceTileID:
        cube.data[itype, iceLoc[0], iceLoc[1]] = 0.0
    # check sum of ice tiles
    cubeTotalIce = cube[iceTileID].collapsed("pseudo_level", iris.analysis.SUM)
    print(
        "check sum ice frac", np.ma.max(cubeTotalIce.data), np.ma.min(cubeTotalIce.data)
    )
    if np.all(np.isclose(cubeTotalIce.data, 0)) or np.all(
        np.isclose(cubeTotalIce.data, 1)
    ):
        raise ValueError("look at regridding of frac - errors in ice")
    # cube.data[iceTileID, iceLoc[0], iceLoc[1]] = 1.0  # set ice to 1

    # where bare soil fraction is 1.0 set values in cube.data  to 0.0
    epsilon = 1e-4
    cubeBare = cube[bareSoilTileID].copy()
    bareSoilLoc = np.where(cubeBare.data > 1.0 - epsilon)
    cubeBare.data[bareSoilLoc[0], bareSoilLoc[1]] = 1.0
    cube.data[:, bareSoilLoc[0], bareSoilLoc[1]] = 0.0
    cube.data[bareSoilTileID, :, :] = cubeBare.data[:, :]

    # mask by land frac
    mask = landMask.data < 0.5
    broadcasted_mask = np.broadcast_to(mask, cube.data.shape)
    cube.data = np.ma.masked_array(cube.data, mask=broadcasted_mask)

    # set bare soil frac to 1 at points
    # what is frac summed over all tiles ?
    cubesum = cube.collapsed("type", iris.analysis.SUM)
    print("check sum frac", np.ma.max(cubesum.data), np.ma.min(cubesum.data))
    cubediff = cubesum.copy()
    cubediff.data = 1.0 - cubesum.data
    cube = cube.copy()  # odd requirement
    cube.data[bareSoilTileID, :, :] = (
        cube.data[bareSoilTileID, :, :] + cubediff.data[:, :]
    )
    cubesum = cube.collapsed("type", iris.analysis.SUM)
    print("recheck sum frac", np.ma.max(cubesum.data), np.ma.min(cubesum.data))
    return cube


# #############################################################
def add_tile_id_for_snow(cube):
    tileDimName = "pseudo_level"
    dim2_name = cube.dim_coords[1].long_name
    ndim2 = len(cube.coord(dim2_name).points)  # tile
    nsnow = len(cube.coord("snow").points)
    orderID = np.linspace(1, nsnow * ndim2, nsnow * ndim2).reshape(nsnow, ndim2)
    order_coord = iris.coords.AuxCoord(orderID, long_name="_pseudo_level_order")
    cube.add_aux_coord(order_coord, data_dims=[0, 1])

    pftSnowID = np.zeros([nsnow, ndim2])
    for i in range(0, nsnow):
        pftSnowID[i, :] = (
            cube.dim_coords[0].points[i] + cube.coord("pseudo_level").points * 1000.0
        )
    snowtile = iris.coords.AuxCoord(
        pftSnowID, long_name="snowtile", var_name="snowtile"
    )
    cube.add_aux_coord(snowtile, data_dims=[0, 1])

    outCubeList = iris.cube.CubeList([])
    for i in range(0, nsnow):
        # print("nsnow = ", i)
        typeCube = cube[i]
        iris.util.promote_aux_coord_to_dim_coord(typeCube, "_pseudo_level_order")
        typeCube.remove_coord("snow")
        typeCube.remove_coord(dim2_name)
        outCubeList.append(typeCube)
    outCube = outCubeList.concatenate_cube()
    print(outCube)

    # re-order the tiles
    dataIn = outCube.data
    tileDimID = outCube.coord(tileDimName).points
    for i in range(0, ndim2):
        j = np.where(
            cube.coord("pseudo_level").points == cube.coord("pseudo_level").points[i]
        )[0][0]
        outCube.data[i : i + nsnow, :, :] = dataIn[j : j + nsnow, :, :]
        outCube[i : i + nsnow].coord(tileDimName).points = tileDimID[j : j + nsnow]
    # print(outCube.coord(tileDimName).points)
    return outCube


# #############################################################
def add_tile_id(cube):
    dim1_name = cube.dim_coords[0].long_name
    ndim1 = len(cube.coord(dim1_name).points)
    # add an auxiliary coordinate containing tile IDs
    # re-order the tiles
    dataIn = cube.data
    for i in range(0, ndim1):
        j = np.where(
            cube.coord("pseudo_level").points == cube.coord("pseudo_level").points[i]
        )[0][0]
        cube.data[i, :, :] = dataIn[j, :, :]
    # rename the z-dimension
    cube.coord(dim1_name).long_name = "_pseudo_level_order"
    return cube


# #############################################################
def main(julesDump, pwdin, pwdout, stashsel="all"):
    """'run"""

    # if no stashcodes are specified make the ancils for stashall
    if stashsel == "all":
        stashsel = stashall
    print(stashsel)

    # get masks for land, soil and ice areas
    landMaskFile, soilMaskFile = define_masks(oceanGrid)
    landMask, iceMask, soilMask = get_masks(landMaskFile, soilMaskFile)

    # loop over each output variable
    for i, stash in enumerate(stashsel):
        print(stash, DICT_STASH[stash])

        # select mask:-- soil,ice or land (land default)
        maskCube = choose_mask(
            soilMask,
            iceMask,
            landMask,
            var_mask_choice=DICT_STASH[stash].get("land_soil_ice_mask"),
        )

        # read the input data from the jules dump
        name = DICT_STASH[stash]["name"]
        if name.startswith("ns_"):
            name = "ns"
        if name.startswith("cs_"):
            name = "cs"
        variable_cons = iris.Constraint(cube_func=lambda x: x.var_name == name)

        cube = jules.load_cube(
            pwdin + julesDump, variable_cons, missingdata=np.ma.masked
        )
        cube = iris.util.squeeze(cube)
        if DICT_STASH[stash].get("pool") is not None:
            pool_cons = iris.Constraint(
                scpool=lambda v: v == DICT_STASH[stash].get("pool") - 1
            )
            cube = cube.extract(pool_cons)

        all_coord_names = [coord.name() for coord in cube.coords()]
        if "type" in all_coord_names:
            mapped_points = [type_info[point] for point in cube.coord("type").points]
            aux_coord = iris.coords.AuxCoord(mapped_points, long_name="pseudo_level")
            cube.add_aux_coord(aux_coord, all_coord_names.index("type"))
        if "tile" in all_coord_names:
            mapped_points = [type_info[point] for point in cube.coord("tile").points]
            aux_coord = iris.coords.AuxCoord(mapped_points, long_name="pseudo_level")
            cube.add_aux_coord(aux_coord, all_coord_names.index("tile"))
        if "pft" in all_coord_names:
            mapped_points = [type_info[point] for point in cube.coord("pft").points]
            aux_coord = iris.coords.AuxCoord(mapped_points, long_name="pseudo_level")
            cube.add_aux_coord(aux_coord, all_coord_names.index("pft"))

        # fraction areas need special treatment
        if stash == "m01s00i216":
            cube = process_frac(cube, maskCube, iceMask, landMask)

        # for snow need to modify the tile ids
        if "snow" in all_coord_names:
            cube = add_tile_id_for_snow(cube)

        # if the variable is 3-dimensional add tile IDs and re-order the tile dimension
        if (cube.dim_coords[0].long_name in ["pft", "type", "tile"]) or (
            cube.dim_coords[1].long_name in ["pft", "type", "tile"]
        ):
            cube = add_tile_id(cube)

        # check if the input and output horizontal grids match
        # if not re-grid the data probably already done
        cube = check_and_regrid(cube, maskCube)

        # deal with missing data points  and mask out the sea
        filler = ants.analysis.FillMissingPoints(cube)
        filler(cube)
        cube.data = np.ma.masked_array(cube.data, mask=None)
        landMask.data = np.ma.masked_less(landMask.data, 0.5)
        cube.data.mask = landMask.data.mask.copy()

        # add key meta-data for ancils
        cube.attributes["STASH"] = iris.fileformats.pp.STASH.from_msi(stash)
        cube.attributes["grid_staggering"] = 6

        # save the data as an ancil and as a netCDF
        outFile = pwdout + julesDump[0:-3] + "." + oceanGrid + "." + stash
        ants.io.save.netcdf(cube, outFile + ".nc", fill_value=-9999)
        # ants.io.save.ancil(cube, outFile + ".anc")


# #############################################################
if __name__ == "__main__":
    """run"""
    pwdin = "/hpc/data/d00/hadea/jules_output/u-cc669_isimip3a_es/"
    # pwdin = "/home/h03/hadea/"
    pwdout = "/scratch/hadea/"
    julesDump = "isimip3a_ES_spinup_03.dump.18070101.0.nc"
    # julesDump = "dump27.nc"

    # get length of type array
    variable_cons = iris.Constraint(cube_func=lambda x: x.var_name == "frac")
    fraccube = jules.load_cube(
        pwdin + julesDump, variable_cons, missingdata=np.ma.masked
    )
    if len(fraccube.coord("type").points) == 17:
        type_info = JULES_ES_MAP
    elif len(fraccube.coord("type").points) == 27:
        type_info = JULES_ES_MAP_ICE
    else:
        raise ValueError("need to sort out landcover types")

    main(julesDump, pwdin, pwdout)  # , ["m01s00i443"])
