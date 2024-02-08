""" fucntions for extracting um data"""
import iris
import os
import iris.coords as icoords

REGION_DICT = {
    None: {"string": "global", "constraint": None},
    "noAntarctica": {
        "string": "noAntarctica",
        "constraint": iris.Constraint(latitude=lambda v: -58.0 <= v <= 85.0),
    },
    "nhlat": {
        "string": "nhlat",
        "constraint": iris.Constraint(latitude=lambda v: 55.0 <= v <= 85.0),
    },
}


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


# ############################################################
def sortout_initial_cs_and_ns(pool_name, cubelist):
    """sort out cpools and npools"""
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
            poolcoord = icoords.DimCoord(ipool, long_name="scpool", units=None)
            # add new dimenson to cube
            cube = iris.util.new_axis(cube)
            cube.var_name = pool_name
            cube.long_name = pool_name
            cube.add_dim_coord(poolcoord, 0)
            cubelist_pools.append(cube)
    iris.util.equalise_attributes(cubelist_pools)
    # print(cubelist_pools)
    cube = cubelist_pools.concatenate_cube()
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
        if "depth" in all_coord_names:
            cube.coord("depth").rename("soil")
        if "snow" in cube.var_name and "pseudo_level" in all_coord_names:
            cube.coord("pseudo_level").rename("tile")
        if "lai" in cube.var_name and "pseudo_level" in all_coord_names:
            cube.coord("pseudo_level").rename("pft")
        if "canht" in cube.var_name and "pseudo_level" in all_coord_names:
            cube.coord("pseudo_level").rename("pft")
        # rest of cubes with pseudo in
    for cube in cubelist:
        all_coord_names = [coord.name() for coord in cube.coords()]
        if "pseudo_level" in all_coord_names:
            cube.coord("pseudo_level").rename("tile")
    return cubelist
