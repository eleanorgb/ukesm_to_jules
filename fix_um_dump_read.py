import itertools

from iris.cube import CubeList
from iris.fileformats.pp import load_pairs_from_fields
from iris.fileformats.um import um_to_pp
import iris
from functions_um_to_jules import stash_integer_to_string

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

# Keep only cubes whose STASH is in the desired sets.
_keep_stash = set(
    STASHPRESOUT
    + STASHTOP
    + STASHSOIL
    + STASHOROG
    + STASHGRID
    + STASHRIVERS
    + STASHLIGHTNING
    + STASHINIT
)

filepath = "/data/users/eleanor.burke/um_to_jules/u-cy838/a.da/cy838a.da19500101_00"
fields_iter = um_to_pp(filepath)

# Just for demonstration, take only first 500 fields as file is really big.
# n_fields_limit = 2000
# fields_iter = itertools.islice(fields_iter, n_fields_limit)

# Remove additional orography fields
selected_fields = []
found_orog = False
for field in fields_iter:
    stash_minor = field.lbuser[3]
    if stash_minor == 33:
        if found_orog:
            # Skip additional orography fields after the first one.
            continue
        found_orog = True
    if stash_integer_to_string(field.lbuser[3]) in _keep_stash:
        #print(field.lbuser[3], stash_integer_to_string(field.lbuser[3]))
        selected_fields.append(field)
print(f"Selected {len(selected_fields)} fields after filtering for orography and getting stash required.")

# Convert fields to 2d cubes, as from "iris.load_raw"
raw_cubes = CubeList([cube for cube, field in load_pairs_from_fields(selected_fields)])
print(f"Converted to {len(raw_cubes)} raw cubes.")



# Make the set of STASH names available for inspection/reuse.
stash_names = {str(cube.attributes.get("STASH")) for cube in raw_cubes}
print(f"Unique STASH names ({len(stash_names)}): {sorted(stash_names)}")

# Merge cubes STASH-by-STASH.
merged_by_stash = CubeList([])
for stash_name in sorted(stash_names):
    cubes_for_stash = CubeList(
        [cube for cube in raw_cubes if str(cube.attributes.get("STASH")) == stash_name]
    )
    if len(cubes_for_stash) == 1:
        merged_by_stash.append(cubes_for_stash[0])
        continue

    try:
        merged_stash_cubes = cubes_for_stash.merge()
        merged_by_stash.extend(merged_stash_cubes)
        print(
            f"Merged STASH {stash_name}: {len(cubes_for_stash)} -> {len(merged_stash_cubes)}"
        )
    except Exception as err:
        print(
            f"[WARNING] Could not merge STASH {stash_name} ({len(cubes_for_stash)} cubes): {err}"
        )
        # Fall back to keeping cubes unmerged for this STASH.
        # merged_by_stash.extend(cubes_for_stash)

raw_cubes = merged_by_stash
print(f"Prepared {len(raw_cubes)} cubes after per-STASH merging.")

# N.B. we now have "raw cubes" (including an orography cube)
#  - use a merge operation to replicate what a plain "iris.load" would do.
# print(raw_cubes)
# iris.save(raw_cubes, "/data/users/eleanor.burke/um_to_jules/u-cy838/a.da/cy838a.da19500101_00_raw.nc")

cubes = raw_cubes.merge()
cubes = CubeList([
    cube for cube in cubes
    if not any(cube.coords(c) for c in ("forecast_period_0", "forecast_period", "forecast_period_1"))
])
print(cubes)
iris.save(
    cubes,
    "/data/users/eleanor.burke/um_to_jules/u-cy838/a.da/cy838a.da19500101_00.sel.nc",
)

# Show results
print(f"Merged {len(raw_cubes)} -> {len(cubes)} cubes :")
print(cubes)
# Show one cube, which has model-levels + orography
print("\nExample result cube :")
print(cubes.extract_cube("air_potential_temperature"))
