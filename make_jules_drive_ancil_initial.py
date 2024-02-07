import glob
import iris
import numpy as np
import jules
from functions_um_to_jules import extract_constraint
from functions_um_to_jules import sort_by_month_year_key
from functions_um_to_jules import rename_cubes
from functions_um_to_jules import make_output_file_name
from functions_um_to_jules import sortout_initial_cs_and_ns
from functions_um_to_jules import sortout_initial_sth
from functions_um_to_jules import rename_and_delete_dimensions
import um_to_jules_stash_dict

L_USE_ROSE = False

if not L_USE_ROSE:
    STREAM = "a.pk"
    RUNID = "dc429"
    PWDIN = "/scratch/hadea/um_to_jules/data/"
    PWDOUT = "/scratch/hadea/um_to_jules/" + RUNID + "/"
    REGION_TO_EXTRACT = None  # "noAntarctica"  # "nhlat" "noAntarctica" None
    DUMPNAME = "cz700a.da20990101_00"


STASHDRIVE = [
    "m01s00i010",
    "m01s30i111",
    "m01s30i417",  # q,T,p
    "m01s02i207",
    "m01s01i235",  # lwf,swf
    "m01s04i203",
    "m01s04i204",
    "m01s05i205",
    "m01s05i206",
    "m01s05i216",
    "m01s05i215",  # pr not all available
    "m01s03i230",  # wind
    "m01s01i231",  # diffuse frac
]
STASHTOP = ["m01s00i274", "m01s00i275", "m01s00i276"]

STASHSOIL = [
    "m01s00i040",
    "m01s00i041",
    "m01s00i043",
    "m01s00i044",
    "m01s00i046",
    "m01s00i047",
    "m01s00i048",
    "m01s00i220",
    "m01s00i418",
    "m01s00i207",
]

STASHOROG = ["m01s00i004", "m01s00i033"]

STASHGRID = ["m01s00i505", "m01s00i157", "m01s00i030"]

STASHRIVERS = ["m01s00i151", "m01s00i152", "m01s00i153"]

STASHINIT = [
    "m01s00i442",
    "m01s00i443",
    "m01s00i444",
    "m01s00i445",
    "m01s00i446",
    "m01s00i466",
    "m01s00i467",
    "m01s00i468",
    "m01s00i469",
    "m01s00i287",
    "m01s00i288",
    "m01s00i289",
    "m01s00i449",
    "m01s00i459",  # ones before this are for triffid
    "m01s00i281",
    "m01s00i278",
    "m01s00i229",
    "m01s00i233",
    "m01s00i214",
    "m01s00i215",
    "m01s00i020",
    "m01s00i240",
    "m01s00i216",
    "m01s00i217",
    "m01s00i218",
    "m01s03i462",
    "m01s00i213",
    "m01s02i31",
    "m01s00i377",
    "m01s00i376",
    "m01s00i242",
    "m01s00i380",
    "m01s00i381",
    "m01s00i382",
    "m01s00i383",
    "m01s00i384",
    "m01s00i386",
]


DICT_STASH = um_to_jules_stash_dict.get_jules_stash_dict()

region_and_stash_cons = extract_constraint(STASHDRIVE, REGION_TO_EXTRACT)


# ##############################################################################
def make_driving():
    search_pattern = f"{PWDIN}{RUNID}/{STREAM}/{RUNID}{STREAM}*.pp"
    print("files searched: ", search_pattern)
    files = glob.glob(search_pattern)
    sorted_files = sorted(files, key=sort_by_month_year_key)

    for input_filename in sorted_files[0:121]:
        print("~~~")
        cubelist_met = iris.load(input_filename)
        cubelist_met = cubelist_met.extract(region_and_stash_cons)
        cubelist_met = rename_cubes(DICT_STASH, cubelist_met)
        # for cube in cubelist_met:
        #    print(cube.long_name)
        output_filename = make_output_file_name(
            input_filename, REGION_TO_EXTRACT, PWDOUT
        )
        iris.save(cubelist_met, output_filename)
        print(output_filename)


# ##############################################################################
def make_topmodel(cubelist_dump):
    """make topmodel ancillaries"""
    region_and_stash_cons = extract_constraint(STASHTOP, REGION_TO_EXTRACT)
    cubelist_top = cubelist_dump.extract(region_and_stash_cons)
    cubelist_top = rename_cubes(DICT_STASH, cubelist_top)
    iris.save(
        cubelist_top,
        PWDOUT + "ancils/" + RUNID + "_" + REGION_TO_EXTRACT + "_topmodel.nc",
    )
    print(cubelist_top)
    return


def make_soil(cubelist_dump):
    """make soil ancillaries"""
    region_and_stash_cons = extract_constraint(STASHSOIL, REGION_TO_EXTRACT)
    cubelist_soil = cubelist_dump.extract(region_and_stash_cons)
    cubelist_soil = rename_cubes(DICT_STASH, cubelist_soil)
    iris.save(
        cubelist_soil,
        PWDOUT + "ancils/" + RUNID + "_" + REGION_TO_EXTRACT + "_soilprops.nc",
    )
    print(cubelist_soil)
    return


def make_rivers(cubelist_dump):
    """make rivers ancillaries - not working yet"""
    region_and_stash_cons = extract_constraint(STASHRIVERS, REGION_TO_EXTRACT)
    cubelist_rivers = cubelist_dump.extract(region_and_stash_cons)
    cubelist_rivers = rename_cubes(DICT_STASH, cubelist_rivers)
    iris.save(
        cubelist_rivers,
        PWDOUT + "ancils/" + RUNID + "_" + REGION_TO_EXTRACT + "_rivers.nc",
    )
    print(cubelist_rivers)
    return


def make_model_height(cubelist_dump):
    """find  height of lowest model level - use stash 004
    https://www-avd/nbviewer/MetOffice/file_/home/h01/hadtq/lib/ipy/model_level_height.ipynb
    """
    region_and_stash_cons = extract_constraint(STASHOROG, REGION_TO_EXTRACT)
    region_and_stash_cons = region_and_stash_cons & iris.Constraint(
        cube_func=lambda c: c.cell_methods == ()
    )
    cubelist_orog = cubelist_dump.extract(region_and_stash_cons)
    cubelist_orog = rename_cubes(DICT_STASH, cubelist_orog)
    cube_theta = cubelist_orog.extract_cube(iris.Constraint(name="theta_4"))
    cube_theta1 = cube_theta.extract(
        iris.Constraint(model_level_number=1)
    )  ## extract level 1 height
    cube_height = cube_theta1.copy()
    cube_height.units = "m"
    cube_height.rename("model_level1_height")
    cube_height.data = (
        cube_theta1.coord("altitude").points
        - cube_theta1.coord("surface_altitude").points
    )
    iris.save(
        cube_height, PWDOUT + "ancils/" + RUNID + "_" + REGION_TO_EXTRACT + "_height.nc"
    )
    return


def make_model_grid(cubelist_dump):
    """get model grid
    something wrong with grid area"""

    region_and_stash_cons = extract_constraint(STASHGRID, REGION_TO_EXTRACT)
    region_and_stash_cons = region_and_stash_cons & iris.Constraint(
        cube_func=lambda c: c.cell_methods == ()
    )
    cubelist_grid = cubelist_dump.extract(region_and_stash_cons)
    cubelist_grid = rename_cubes(DICT_STASH, cubelist_grid)

    lonlat_cube = cubelist_grid.extract_cube(
        iris.Constraint(cube_func=lambda x: x.attributes["STASH"] == "m01s00i030")
    )
    lsmask = lonlat_cube.copy()
    lat2d = lonlat_cube.copy()
    lat2d.data = np.tile(
        lat2d.coord("latitude").points, (len(lat2d.coord("longitude").points), 1)
    ).transpose()
    lat2d.var_name = "lat2d"
    lat2d.long_name = "lat2d"
    lat2d.standard_name = None
    cubelist_grid.append(lat2d)

    lon2d = lonlat_cube.copy()
    lon2d.data = np.tile(
        lon2d.coord("longitude").points, (len(lon2d.coord("latitude").points), 1)
    )
    lon2d.var_name = "lon2d"
    lon2d.long_name = "lon2d"
    lon2d.standard_name = None
    cubelist_grid.append(lon2d)

    iris.save(
        cubelist_grid,
        PWDOUT + "ancils/" + RUNID + "_" + REGION_TO_EXTRACT + "_model_grid.nc",
    )
    return lsmask


def make_initial_conditions(cubelist_dump, lsmask):
    region_and_stash_cons = extract_constraint(STASHINIT, REGION_TO_EXTRACT)
    cubelist_init = cubelist_dump.extract(region_and_stash_cons)
    cubelist_init = rename_cubes(DICT_STASH, cubelist_init)
    # print(cubelist_init)
    # sort out soil moisture
    cube = sortout_initial_sth(cubelist_init)
    cubelist_init.append(cube)
    # sort out cpools and npools
    for pool_name in ["ns", "cs"]:
        cube = sortout_initial_cs_and_ns(pool_name, cubelist_init)
        cubelist_init.append(cube)
    l_remove_time = True
    cubelist_init = rename_and_delete_dimensions(cubelist_init, l_remove_time)

    jules.save(
        cubelist_init,
        PWDOUT + "dump/" + RUNID + "_" + REGION_TO_EXTRACT + "_dump.nc",
        lsmask=lsmask,
        landpointsonly=True,
    )
    return


if __name__ == "__main__":
    make_driving()

    cubelist_dump = iris.load(PWDIN + RUNID + "/dump/" + DUMPNAME)
    make_topmodel(cubelist_dump)
    make_rivers(cubelist_dump)
    make_soil(cubelist_dump)
    make_model_height(cubelist_dump)
    lsmask = make_model_grid(cubelist_dump)
    make_initial_conditions(cubelist_dump, lsmask)