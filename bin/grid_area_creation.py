#!/usr/bin/env python3
#%%
from os import walk, path, system
from pyprojroot.here import here

MSWEP_folder = here("Data/Data_MSWEP")
MSWEP_paths = [
    path.join(dirpath, f)
    for (dirpath, dirnames, filenames) in walk(MSWEP_folder)
    for f in filenames
]
# select one file
MSWEP_path = MSWEP_paths[0]
# %%
# create the grid cell area file
command = f"cdo gridarea {MSWEP_path} {here('Data')}/grid_area_MSWEP.nc"
system(command)

# %%
