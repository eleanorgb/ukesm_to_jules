""" fucntions for extracting um data"""
import iris
import os

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
    if region != None:
        region_and_stash_cons = region.copy()
        region_and_stash_cons = (
            region_and_stash_cons[REGION_TO_EXTRACT]["constraint"] & stash_cons
        )
    else:
        region_and_stash_cons = stash_cons
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
def make_output_file_name(input_filename, REGION_TO_EXTRACT):
    output_filename = (
        PWDOUT + "drive/" + os.path.splitext(os.path.basename(input_filename))[0]
    )
    if region != None:
        output_filename += "_" + region[REGION_TO_EXTRACT]["string"] + ".nc"
    else:
        output_filename += ".nc"
    return output_filename


# ############################################################
def rename_cubes(DICT_STASH, cubelist):
    """this removes standard names may or may not be a good thing"""
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
        raise Exception("check this case: may need to do some merging/concatenating")
    return cubelist
