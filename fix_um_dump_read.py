#import itertools

from iris.cube import CubeList
from iris.fileformats.pp import load_pairs_from_fields
from iris.fileformats.um import um_to_pp
import iris

filepath = '/data/users/eleanor.burke/um_to_jules/u-cy838/a.da/cy838a.da19500101_00'
fields_iter = um_to_pp(filepath)

# Just for demonstration, take only first 500 fields as file is really big.
#n_fields_limit = 500
#fields_iter = itertools.islice(fields_iter, n_fields_limit)

# Remove duplicate fields: keep first occurrence of each STASH code, skip duplicates
selected_fields = []
seen_stash = set()
for field in fields_iter:
    stash_code = field.lbuser[3]  # minor STASH code
    # Also consider major section for full uniqueness
    stash_section = field.lbuser[0]
    full_stash = (stash_section, stash_code)
    
    if full_stash in seen_stash:
        # Skip duplicate
        continue
    seen_stash.add(full_stash)
    selected_fields.append(field)

# Convert fields to 2d cubes, as from "iris.load_raw"
raw_cubes = CubeList([cube for cube, field in load_pairs_from_fields(selected_fields)])

# N.B. we now have "raw cubes" (including an orography cube)
#  - use a merge operation to replicate what a plain "iris.load" would do.
cubes = raw_cubes.merge()

iris.save(cubes, "/data/users/eleanor.burke/um_to_jules/u-cy838/a.da/cy838a.da19500101_00.nc")

# Show results
print(f"Merged {len(raw_cubes)} -> {len(cubes)} cubes :")
print(cubes)
# Show one cube, which has model-levels + orography
print("\nExample result cube :")
print(cubes.extract_cube("air_potential_temperature"))
