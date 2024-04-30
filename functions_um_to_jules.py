""" fucntions for extracting um data"""
import iris
import os
import numpy as np
import iris.coords as icoords

REGION_DICT = {
    "global": {"string": "global", "constraint": None},
    "noAntarctica": {
        "string": "noAntarctica",
        "constraint": iris.Constraint(latitude=lambda v: -58.0 <= v <= 85.0),
    },
    "nhlat": {
        "string": "nhlat",
        "constraint": iris.Constraint(latitude=lambda v: 55.0 <= v <= 85.0),
    },
}

# ############################################################


# ############################################################
def define_mapping_pseudo_type_to_type(ntypes):
    if ntypes == 13:   # pfts only, not other types
        mapping_pseudo_type_to_type = {
            3: 6,
            4: 9,
            101: 1,
            102: 2,
            103: 3,
            201: 4,
            202: 5,
            301: 7,
            302: 8,
            401: 10,
            402: 11,
            501: 12,
            502: 13,
        }
    elif ntypes == 27:   #elevation land ice tiles
        mapping_pseudo_type_to_type = {
            3: 6,
            4: 9,
            6: 14,
            7: 15,
            8: 16,
            9: 17,
            101: 1,
            102: 2,
            103: 3,
            201: 4,
            202: 5,
            301: 7,
            302: 8,
            401: 10,
            402: 11,
            501: 12,
            502: 13,
            901: 18,
            902: 19,
            903: 20,
            904: 21,
            905: 22,
            906: 23,
            907: 24,
            908: 25,
            909: 26,
            910: 27,
        }
    else:
        raise ValueError(
            "need to define mapping_pseudo_type_to_type for " + str(ntypes)
        )
    return mapping_pseudo_type_to_type


# ############################################################
def sort_by_month_year_key(filename):
    """sort driving data by year then month"""
    parts = filename.split(".")[-2]
    year = int(parts[-7:-3])
    month = parts[-3:]
    # Define a mapping of month names to their order
    month_order = {
        "jan": 1,
        "feb": 2,
        "mar": 3,
        "apr": 4,
        "may": 5,
        "jun": 6,
        "jul": 7,
        "aug": 8,
        "sep": 9,
        "oct": 10,
        "nov": 11,
        "dec": 12,
    }
    return (year, month_order.get(month.lower(), 0))


# ############################################################
def extract_constraint(stash, REGION_TO_EXTRACT):
    """define constraint for files"""
    stash_cons = iris.Constraint(cube_func=lambda x: x.attributes["STASH"] in stash)
    region_and_stash_cons = REGION_DICT[REGION_TO_EXTRACT]["constraint"] & stash_cons
    return region_and_stash_cons


# ############################################################
def stash_string_to_integer(stash_string=None):
    """e.g. m01s16i222 to 16222. Copied from myutils.py@8221"""

    if len(stash_string) != 10:
        print(len(stash_string), "~" + stash_string + "~")
        raise Exception("something wrong here")

    number_string = stash_string[1:3] + stash_string[4:6] + stash_string[7:]
    stash_integer = int(number_string) - 100000

    return stash_integer


# ############################################################
def make_output_file_name(input_filename, REGION_TO_EXTRACT, PWDOUT, UM_RUNID):
    output_filename = (
        PWDOUT
        + "/u-"
        + UM_RUNID
        + "/drive/"
        + os.path.splitext(os.path.basename(input_filename))[0]
    )
    if REGION_DICT[REGION_TO_EXTRACT]["string"] != None:
        output_filename += "_" + REGION_DICT[REGION_TO_EXTRACT]["string"] + ".nc"
    else:
        output_filename += ".nc"
    return output_filename


# ############################################################
def rename_cubes(DICT_STASH, cubelist):
    """this removes standard names may or may not be a good thing"""
    if isinstance(cubelist, iris.cube.CubeList):
        for cube in cubelist:
            if "STASH" in cube.attributes.keys():
                stash_int = stash_string_to_integer(str(cube.attributes["STASH"]))
            cube.long_name = (
                DICT_STASH[str(cube.attributes["STASH"])]["name"] + "_" + str(stash_int)
            )
            cube.units = DICT_STASH[str(cube.attributes["STASH"])]["units"]
            cube.standard_name = None
            cube.var_name = cube.long_name

        name_list = [cube.name() for cube in cubelist]
        if len(set(name_list)) != len(name_list):
            print(name_list)
            raise Exception(
                "check this case: may need to do some merging/concatenating"
            )
    else:
        if "STASH" in cubelist.attributes.keys():
            stash_int = stash_string_to_integer(str(cubelist.attributes["STASH"]))
        cubelist.long_name = (
            DICT_STASH[str(cubelist.attributes["STASH"])]["name"] + "_" + str(stash_int)
        )
        cubelist.units = DICT_STASH[str(cubelist.attributes["STASH"])]["units"]
        cubelist.standard_name = None
        cubelist.var_name = cubelist.long_name

    return cubelist


# ##############################################################################
def reorder_pseudo_type(cube):
    """"""
    cubelist = iris.cube.CubeList([])
    ntypes = len(cube.coord("pseudo_level").points)
    mapping_pseudo_type_to_type = define_mapping_pseudo_type_to_type(ntypes)
    for ijk in np.arange(0, ntypes):
        cube_tmp = cube[ijk].copy()
        cube_tmp.coord("pseudo_level").points = [
            mapping_pseudo_type_to_type.get(cube_tmp.coord("pseudo_level").points[0])
        ]
        cubelist.append(cube_tmp)
    cube = cubelist.merge_cube()
    return cube


# ##############################################################################
def sortout_initial_snow(snow_cube, nsnow, tile_coord):
    """"""
    cube = iris.cube.Cube(
        np.zeros(
            (
                nsnow,
                tile_coord.shape[0],
                snow_cube.coord("latitude").shape[0],
                snow_cube.coord("longitude").shape[0],
            )
        ),
        var_name=snow_cube.var_name,
        units=snow_cube.units,
        dim_coords_and_dims=[
            (iris.coords.DimCoord(np.arange(0, nsnow), long_name="snow"), 0),
            (tile_coord, 1),
            (snow_cube.coord("latitude"), 2),
            (snow_cube.coord("longitude"), 3),
        ],
    )
    for ijk in np.arange(0, tile_coord.shape[0]):
        cube.data[0:nsnow, ijk, :, :] = snow_cube.data[
            (0 + nsnow * ijk) : (nsnow + nsnow * ijk), :, :
        ]
    for attr_name, attr_value in snow_cube.attributes.items():
        cube.attributes[attr_name] = attr_value
    return cube


# ############################################################
def sortout_initial_cs_and_ns(pool_name, cubelist):
    """sort out cpools and npools only works for 1 layer currently"""
    cubelist_pools = iris.cube.CubeList([])
    for cube in cubelist:
        if cube.long_name.startswith(pool_name + "_"):
            if "dpm" in cube.long_name:
                ipool = 1
            elif "rpm" in cube.long_name:
                ipool = 2
            elif "bio" in cube.long_name:
                ipool = 3
            elif "hum" in cube.long_name:
                ipool = 4
            pool_coord = icoords.DimCoord(ipool, long_name="scpool", units=None)
            # add new dimenson to cube
            cube = iris.util.new_axis(cube)
            cube.var_name = pool_name
            cube.long_name = pool_name
            cube.add_dim_coord(pool_coord, 0)
            cubelist_pools.append(cube)
    iris.util.equalise_attributes(cubelist_pools)
    # print(cubelist_pools)
    cube = cubelist_pools.concatenate_cube()
    # add sclayer coordinate
    cube = iris.util.new_axis(cube)
    sclayer_coord = icoords.DimCoord(1, long_name="sclayer", units=None)
    cube.add_dim_coord(sclayer_coord, 0)
    cube.transpose([1, 0, 2, 3])

    return cube


# ############################################################
def sortout_initial_sth(cubelist):
    """sort out soil moisture"""
    cubelist_sthuf = cubelist.extract(
        iris.Constraint(
            cube_func=lambda x: x.attributes["STASH"] in ["m01s00i214", "m01s00i215"]
        )
    )
    if len(cubelist_sthuf) != 2:
        raise ValueError("need 2 cubes here but dont have that")
    cube = cubelist_sthuf[0] + cubelist_sthuf[1]
    cube.var_name = "sthuf"
    cube.long_name = "sthuf"
    return cube


# ############################################################
def rename_and_delete_dimensions(cubelist, l_remove_time=False):
    """sort out dimensions"""
    for cube in cubelist:
        all_coord_names = [coord.name() for coord in cube.coords()]
        if "realization" in all_coord_names:
            cube.remove_coord("realization")
        if l_remove_time:
            if "time" in all_coord_names:
                cube.remove_coord("time")
        if "pseudo_level" in all_coord_names:
            if len(cube.coord("pseudo_level").points) == 1:
                cube.remove_coord("pseudo_level")
            else:
                cube.coord("pseudo_level").points = np.arange(
                    0, cube.coord("pseudo_level").shape[0]
                )
        if "depth" in all_coord_names:
            cube.coord("depth").rename("soil")
        if "snow" in cube.var_name and "pseudo_level" in all_coord_names:
            cube.coord("pseudo_level").rename("tile")
        if "lai" in cube.var_name and "pseudo_level" in all_coord_names:
            cube.coord("pseudo_level").rename("pft")
        if "canht" in cube.var_name and "pseudo_level" in all_coord_names:
            cube.coord("pseudo_level").rename("pft")
        if "frac" in cube.var_name and "pseudo_level" in all_coord_names:
            cube.coord("pseudo_level").rename("type")
        # rest of cubes with pseudo in
    for cube in cubelist:
        all_coord_names = [coord.name() for coord in cube.coords()]
        if "pseudo_level" in all_coord_names:
            cube.coord("pseudo_level").rename("tile")
    return cubelist
