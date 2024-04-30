import sys
import glob
import iris
import numpy as np
import jules
import os
import iris.coords as icoords
from functions_um_to_jules import extract_constraint
from functions_um_to_jules import sort_by_month_year_key
from functions_um_to_jules import rename_cubes
from functions_um_to_jules import make_output_file_name
from functions_um_to_jules import sortout_initial_cs_and_ns
from functions_um_to_jules import sortout_initial_sth
from functions_um_to_jules import sortout_initial_snow
from functions_um_to_jules import rename_and_delete_dimensions
from functions_um_to_jules import reorder_pseudo_type
from functions_um_to_jules import REGION_DICT
import um_to_jules_stash_dict

L_USE_ROSE = True
NSNOW = 10

if not L_USE_ROSE:
    UM_RUNID = "dc429"
    PWDUSE = "/scratch/hadea/um_to_jules/"
    REGION_TO_EXTRACT = "noAntarctica"  # "nhlat" "noAntarctica" "global"
    UM_DUMPFILENAME = "cz700a.da20990101_00"
    START_YEAR = 1
    END_YEAR = 1
    MAKE_INIT_ANCIL = True
    MAKE_INIT_DUMP = True

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
    "m01s05i214",  # pr not all available
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

STASHLIGHTNING = ["m01s50i082", "m01s21i097"]

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
    "m01s00i231",
    "m01s00i576",
    "m01s00i377",
    "m01s00i376",
    "m01s00i242",
    "m01s00i380",
    "m01s00i381",
    "m01s00i382",
    "m01s00i383",
    "m01s00i384",
    "m01s00i386",
    "m01s00i418",
]

STASHPRESOUT = ["m01s19i111", "m01s00i252"]

DICT_STASH = um_to_jules_stash_dict.get_jules_stash_dict()

FILE_WITH_ANCIL_PATHS = "ancilpaths.dat"

# check if region to extract is in dictionary here


# ##############################################################################
def read_parameters_from_rose():
    # Define parameter names
    parameter_names = [
        "UM_RUNID",
        "PWDUSE",
        "REGION_TO_EXTRACT",
        "UM_DUMPFILENAME",
        "START_YEAR",
        "END_YEAR",
        "MAKE_INIT_ANCIL",
        "MAKE_INIT_DUMP",
    ]

    # Check if correct number of arguments is provided
    if len(sys.argv) not in np.arange(
        len(parameter_names) - 2,
        len(parameter_names) + 2,
    ):  # +1 to account for script name
        print(
            f"Usage: python {sys.argv[0]} {' '.join(parameter_names)} [START_YEAR] [END_YEAR] [MAKE_INIT_ANCIL] [MAKE_INIT_DUMP]"
        )
        sys.exit(1)

    parameters_dict = {}
    for arg in sys.argv[1:]:
        key, value = arg.split("=")
        parameters_dict[key.strip()] = value.strip()

    UM_RUNID = parameters_dict["UM_RUNID"]
    PWDUSE = parameters_dict["PWDUSE"]
    REGION_TO_EXTRACT = parameters_dict["REGION_TO_EXTRACT"]
    UM_DUMPFILENAME = parameters_dict["UM_DUMPFILENAME"]
    START_YEAR = parameters_dict.get("START_YEAR", 1)
    END_YEAR = parameters_dict.get("END_YEAR", 1)
    MAKE_INIT_ANCIL_str = parameters_dict.get(
        "MAKE_INIT_ANCIL", "False"
    )  # Default to "False" if not provided
    MAKE_INIT_DUMP_str = parameters_dict.get(
        "MAKE_INIT_DUMP", "False"
    )  # Default to "False" if not provided

    MAKE_INIT_ANCIL = MAKE_INIT_ANCIL_str.lower() in ["true", "t", "yes", "y"]
    # Check if MAKE_INIT_ANCIL is a boolean
    if not isinstance(MAKE_INIT_ANCIL, bool):
        print("MAKE_INIT_ANCIL must be either True or False.")
        sys.exit(1)

    MAKE_INIT_DUMP = MAKE_INIT_DUMP_str.lower() in ["true", "t", "yes", "y"]
    # Check if MAKE_INIT_DUMP is a boolean
    if not isinstance(MAKE_INIT_DUMP, bool):
        print("MAKE_INIT_DUMP must be either True or False.")
        sys.exit(1)

    if (START_YEAR == 1 and END_YEAR != 1) or (START_YEAR != 1 and END_YEAR == 1):
        print(
            "START_YEAR and END_YEAR must either both be 1 or both be valid model years."
        )
        sys.exit(1)

    if END_YEAR < START_YEAR:
        print("END_YEAR must be greater than START_YEAR.")
        sys.exit(1)

    return (
        UM_RUNID,
        PWDUSE,
        REGION_TO_EXTRACT,
        UM_DUMPFILENAME,
        START_YEAR,
        END_YEAR,
        MAKE_INIT_ANCIL,
        MAKE_INIT_DUMP,
    )


# ##############################################################################
def make_prescribed_from_input(um_ancillary_filename):
    cubelist = iris.load(um_ancillary_filename)
    stashlist = []
    for cube in cubelist:
        stashlist.append(str(cube.attributes["STASH"]))
    region_and_stash_cons = extract_constraint(stashlist, REGION_TO_EXTRACT)
    cubelist = cubelist.extract(region_and_stash_cons)
    print(cubelist)
    cubelist = rename_cubes(DICT_STASH, cubelist)
    for cube in cubelist:
        print(cube.var_name)
        iris.save(
            cube,
            f"{PWDUSE}/u-{UM_RUNID}/ancils/{UM_RUNID}_{REGION_DICT[REGION_TO_EXTRACT]['string']}_{cube.var_name}.nc",
        )
    return


# ##############################################################################
def make_prescribed_from_output():
    stream = "a.pm"
    search_pattern = f"{PWDUSE}/u-{UM_RUNID}/{stream}/{UM_RUNID}{stream}*.pp"
    print("files searched: ", search_pattern)
    files = glob.glob(search_pattern)
    sorted_files = sorted(files, key=sort_by_month_year_key)
    print(sorted_files)
    region_and_stash_cons = extract_constraint(STASHPRESOUT, REGION_TO_EXTRACT)
    cubelist = iris.load(sorted_files, region_and_stash_cons)
    print(cubelist)
    cubelist = rename_cubes(DICT_STASH, cubelist)
    print(cubelist)
    for cube in cubelist:
        print(cube.var_name)
        # co2_mmr must have the same value everywhere
        if cube.var_name in ["co2_mmr_252"]:
            cube_global = cube.collapsed(["latitude", "longitude"], iris.analysis.MEAN)
            print(cube_global.data[:])
            cube.data[:, :, :] = cube_global.data[:, np.newaxis, np.newaxis]
        iris.save(
            cube,
            f"{PWDUSE}/u-{UM_RUNID}/ancils/{UM_RUNID}_{REGION_DICT[REGION_TO_EXTRACT]['string']}_{cube.var_name}.nc",
        )
    return


# ##############################################################################
def make_driving(year):
    """make the driving data from high temporal resolution"""
    stream = "a.pk"
    if year > 1:
        search_pattern = (
            f"{PWDUSE}/u-{UM_RUNID}/{stream}/{UM_RUNID}{stream}{str(year)}*.pp"
        )
    else:
        search_pattern = f"{PWDUSE}/u-{UM_RUNID}/{stream}/{UM_RUNID}{stream}*.pp"
    print("files searched: ", search_pattern)
    files = glob.glob(search_pattern)
    sorted_files = sorted(files, key=sort_by_month_year_key)

    region_and_stash_cons = extract_constraint(STASHDRIVE, REGION_TO_EXTRACT)

    for input_filename in sorted_files:
        print("~~~")
        cubelist_met = iris.load(input_filename, region_and_stash_cons)
        cubelist_met = rename_cubes(DICT_STASH, cubelist_met)
        # for cube in cubelist_met:
        #    print(cube.long_name)
        output_filename = make_output_file_name(
            input_filename, REGION_TO_EXTRACT, PWDUSE, UM_RUNID
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
        f"{PWDUSE}/u-{UM_RUNID}/ancils/{UM_RUNID}_{REGION_DICT[REGION_TO_EXTRACT]['string']}_topmodel.nc",
    )
    print(cubelist_top)
    return


# ##############################################################################
def make_soil(cubelist_dump):
    """make soil ancillaries"""
    region_and_stash_cons = extract_constraint(STASHSOIL, REGION_TO_EXTRACT)
    cubelist_soil = cubelist_dump.extract(region_and_stash_cons)
    cubelist_soil = rename_cubes(DICT_STASH, cubelist_soil)
    iris.save(
        cubelist_soil,
        f"{PWDUSE}/u-{UM_RUNID}/ancils/{UM_RUNID}_{REGION_DICT[REGION_TO_EXTRACT]['string']}_soilprops.nc",
    )
    print(cubelist_soil)
    return


# ##############################################################################
def make_rivers(cubelist_dump):
    """make rivers ancillaries - not working yet"""
    region_and_stash_cons = extract_constraint(STASHRIVERS, REGION_TO_EXTRACT)
    cubelist_rivers = cubelist_dump.extract(region_and_stash_cons)
    cubelist_rivers = rename_cubes(DICT_STASH, cubelist_rivers)
    iris.save(
        cubelist_rivers,
        f"{PWDUSE}/u-{UM_RUNID}/ancils/{UM_RUNID}_{REGION_DICT[REGION_TO_EXTRACT]['string']}_rivers.nc",
    )
    print(cubelist_rivers)
    return


# ##############################################################################
def make_c2g_lightning(lat2d):
    """
    Gets lightning for input to jules
    if 21097 is available use that WE THINK
    but if not use 50082, only need conversion for 50082
         Converts total flashes (#50082) to cloud-to-ground for INFERNO
    """

    stream = "a.pm"
    search_pattern = f"{PWDUSE}/u-{UM_RUNID}/{stream}/{UM_RUNID}{stream}*.pp"
    print("files searched: ", search_pattern)
    files = glob.glob(search_pattern)
    sorted_files = sorted(files, key=sort_by_month_year_key)
    # print(sorted_files)
    region_and_stash_cons = extract_constraint(STASHLIGHTNING, REGION_TO_EXTRACT)
    cubelist = iris.load(sorted_files, region_and_stash_cons)
    stashlist = []
    for cube in cubelist:
        stashlist.append(str(cube.attributes["STASH"]))
    if "m01s21i097" in stashlist:
        print("use this one")
        cube = cubelist.extract(
            iris.Constraint(cube_func=lambda x: x.attributes["STASH"] in ["m01s21i097"])
        )
    elif "m01s50i082" in stashlist:
        cube = cubelist.extract_cube(
            iris.Constraint(cube_func=lambda x: x.attributes["STASH"] in ["m01s50i082"])
        )

        # Conversion code taken from UM vn 12.0 src/atmosphere/UKCA/ukca_light.F90
        adh = (
            -6.64e-05 * (np.absolute(lat2d.data) ** 2)
            - 4.73e-03 * np.absolute(lat2d.data)
            + 7.34
        )
        az = (
            0.021 * (adh**4)
            - 0.648 * (adh**3)
            + 7.493 * (adh**2)
            - 36.54 * adh
            + 63.09
        )
        ap = 1.0 / (1.0 + az)
        cube.data = cube.data * ap
    else:
        print("no lightning data available in files given")
        return
    cube = rename_cubes(DICT_STASH, cube)
    cube.data = cube.data / 30.0  # convert month to day?

    iris.save(
        cube,
        f"{PWDUSE}/u-{UM_RUNID}/ancils/{UM_RUNID}_{REGION_DICT[REGION_TO_EXTRACT]['string']}_{cube.var_name}.nc",
    )
    return


# ##############################################################################
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
        cube_height,
        f"{PWDUSE}/u-{UM_RUNID}/ancils/{UM_RUNID}_{REGION_DICT[REGION_TO_EXTRACT]['string']}_height.nc",
    )
    return


# ##############################################################################
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

    cubelist_tmp = iris.cube.CubeList([])
    for cube in cubelist_grid:
        if np.ma.is_masked(cube.data):
            cube.data[cube.data.mask] = 0.0
        cubelist_tmp.append(cube)
    cubelist_grid = cubelist_tmp.copy()

    iris.save(
        cubelist_grid,
        f"{PWDUSE}/u-{UM_RUNID}/ancils/{UM_RUNID}_{REGION_DICT[REGION_TO_EXTRACT]['string']}_model_grid.nc",
    )
    return lsmask, lat2d


# ##############################################################################
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

    # sort out n_inorg only worked for non layered cs
    cubelist_tmp = iris.cube.CubeList([])
    for cube in cubelist_init:
        if cube.long_name.startswith("n_inorg_") or cube.long_name.startswith("clay_"):
            cube = iris.util.new_axis(cube)
            sclayer_coord = icoords.DimCoord(1, long_name="sclayer", units=None)
            cube.add_dim_coord(sclayer_coord, 0)
            cubelist_tmp.append(cube)
        else:
            cubelist_tmp.append(cube)
    cubelist_init = cubelist_tmp.copy()

    # sort out layered snow
    cube_frac = cubelist_dump.extract_cube(
        extract_constraint(["m01s00i216"], REGION_TO_EXTRACT)
    )
    ntiles = cube_frac.coord("pseudo_level").shape[0]
    cubelist_tmp = iris.cube.CubeList([])
    for cube in cubelist_init:
        all_coord_names = [coord.name() for coord in cube.coords()]
        if "pseudo_level" in all_coord_names:
            if cube.coord("pseudo_level").shape[0] > ntiles:
                cube = sortout_initial_snow(
                    cube, NSNOW, cube_frac.coord("pseudo_level")
                )
                cubelist_tmp.append(cube)
            elif cube.coord("pseudo_level").shape[0] == ntiles:
                cube = reorder_pseudo_type(cube)
                cubelist_tmp.append(cube)
            else:
                cubelist_tmp.append(cube)
    cubelist_init = cubelist_tmp.copy()

    # sort out tile order for pft
    cubelist_tmp = iris.cube.CubeList([])
    for cube in cubelist_init:
        if cube.var_name in ["canht_218", "lai_217"]:
            cube = reorder_pseudo_type(cube)
            cubelist_tmp.append(cube)
        else:
            cubelist_tmp.append(cube)
    cubelist_init = cubelist_tmp.copy()

    l_remove_time = True
    cubelist_init = rename_and_delete_dimensions(cubelist_init, l_remove_time)

    jules.save(
        cubelist_init,
        f"{PWDUSE}/u-{UM_RUNID}/dump/{UM_RUNID}_{REGION_DICT[REGION_TO_EXTRACT]['string']}_dump.nc",
        lsmask=lsmask,
        landpointsonly=True,
        lsmask_missingdata_str="zeros",
    )
    iris.save(
        cubelist_init,
        f"{PWDUSE}/u-{UM_RUNID}/ancils/{UM_RUNID}_{REGION_DICT[REGION_TO_EXTRACT]['string']}_latlon.nc",
    )
    return


# ##############################################################################
# ##############################################################################
if __name__ == "__main__":
    """run main routine"""

    if L_USE_ROSE:
        # python make_jules_drive_ancil_initial.py UM_RUNID=dc429 PWDUSE=/scratch/hadea/um_to_jules REGION_TO_EXTRACT=noAntarctica UM_DUMPFILENAMEcz700a.da20990101_00
        (
            UM_RUNID,
            PWDUSE,
            REGION_TO_EXTRACT,
            UM_DUMPFILENAME,
            START_YEAR,
            END_YEAR,
            MAKE_INIT_ANCIL,
            MAKE_INIT_DUMP,
        ) = read_parameters_from_rose()

    print("# #####################################################")
    print("UM_RUNID = " + UM_RUNID)
    print("PWDUSE = " + PWDUSE)
    print("REGION_TO_EXTRACT = " + REGION_TO_EXTRACT)
    print("UM_DUMPFILENAME = " + UM_DUMPFILENAME)
    print("START_YEAR = " + str(START_YEAR) + " (NB. 1 is all years)")
    print("END_YEAR = " + str(END_YEAR) + " (NB. 1 is all years)")
    print("MAKE_INIT_ANCIL = " + str(MAKE_INIT_ANCIL))
    print("MAKE_INIT_DUMP = " + str(MAKE_INIT_DUMP))
    print("# #####################################################")

    if not os.path.exists(f"{PWDUSE}/u-{UM_RUNID}"):
        os.makedirs(f"{PWDUSE}/u-{UM_RUNID}")
    for subdir in ["ancils", "dump", "drive"]:
        if not os.path.exists(f"{PWDUSE}/u-{UM_RUNID}/{subdir}"):
            os.makedirs(f"{PWDUSE}/u-{UM_RUNID}/{subdir}")

    if not MAKE_INIT_ANCIL and not MAKE_INIT_DUMP:
        for year in range(int(START_YEAR), int(END_YEAR) + 1):
            make_driving(year)

    if MAKE_INIT_ANCIL:
        # need to get these filenames from somewhere else and need for population
        if (
            os.path.isfile(FILE_WITH_ANCIL_PATHS)
            and os.path.getsize(FILE_WITH_ANCIL_PATHS) > 0
        ):
            with open(FILE_WITH_ANCIL_PATHS, "r") as file:
                for filename in file:
                    filename = (
                        filename.strip()
                    )  # Remove leading/trailing whitespace and newlines
                    make_prescribed_from_input(filename)
        else:
            print("Expecting file ancilpaths.dat in current directory")
            print("It either doesnt exist or has no filenames")
            print("There are no prescribed datasets derived from UM input data")

        # from monthly output from UM run
        make_prescribed_from_output()

    if MAKE_INIT_DUMP:
        cubelist_dump = iris.load(f"{PWDUSE}/u-{UM_RUNID}/a.da/{UM_DUMPFILENAME}")
        make_topmodel(cubelist_dump)
        make_rivers(cubelist_dump)
        make_soil(cubelist_dump)
        make_model_height(cubelist_dump)
        lsmask, lat2d = make_model_grid(cubelist_dump)
        # lat2d = iris.load_cube("/scratch/hadea/um_to_jules/dc429/ancils/dc429_noAntarctica_model_grid.nc","lat2d")
        make_c2g_lightning(lat2d)
        make_initial_conditions(cubelist_dump, lsmask)
