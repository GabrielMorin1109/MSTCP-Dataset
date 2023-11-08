#!/usr/bin/env python3
# %%
from pyprojroot.here import here
import os
from datetime import datetime


def paths_to_rasters() -> dict:
    """
    Returns a list of paths representing the location where the projected rasters will be saved.
    """

    # Get the path from where the rasters is saved
    parent_path = here("Data/Data_MSWEP/")

    # List all files from subfolder of parent_path, and read out-of-memory all rasters
    nc_paths = [
        os.path.join(dirpath, f)
        for (dirpath, dirnames, filenames) in os.walk(parent_path)
        for f in filenames
    ]

    # Extract filenames with regex
    files_name = [
        os.path.splitext(os.path.basename(nc_path))[0] for nc_path in nc_paths
    ]

    # Variant:
    variant = [
        os.path.basename(os.path.normpath(os.path.dirname(nc_path)))
        for nc_path in nc_paths
    ]

    # Get corresponding date from filename
    date = [datetime.strptime(date_string, "%Y%j.%H") for date_string in files_name]

    # Create the dictionary containing all information
    OUT = {
        "origin": nc_paths,
        "files_name": files_name,
        "variant": variant,
        "time": date,
    }
    return OUT


# %%
