# Read Coulomb input files in the .inr format
# Important parameters: 
# 1. Poisson's Ratio
# 2. Coefficient of Friction
# 3. Depth
# 4. Fault params (x start, x finish, y start, y finish, Kode, rt. lat, reverse, dip angle, top, bottom, comment)
# 5. Grid parameters
# 6. Map info (min lon, max lon, zero lon, min lat, max lat, zero lat)

# In this case, we will have to convert slip/rake into strike slip / dip slip. 

from . import coulomb_collections

def read_inr(input_file):
    return [];
