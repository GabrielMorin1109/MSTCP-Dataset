#!/usr/bin/env python3
# %%
# Packages ----------------------------------------------------------------
from pyprojroot.here import here

import geopandas as gpd
import pandas as pd
import numpy as np

from get_paths_to_rasters import paths_to_rasters


# %% =============================================================================================================================
def IBTrACS_read(
    path_IBTrACS: str = "Data/Data_IBTrACS/ibtracs.ALL.list.v04r00.csv",
) -> gpd.GeoDataFrame:
    """
    Import IBTrACS to memory as Geopandas.
    The time (ISO_TIME) is rounded to the nearest 3 hours to match MSWEP.
    For good measure, MSWEP is too rounded to the nearest 3 hours.

    Path to the corresponding raster is stored in the Geopandas DataFrame.
    """
    IBTrACS = pd.read_csv(here(path_IBTrACS), skiprows=[1])
    # Convert to geopandas format
    IBTrACS = gpd.GeoDataFrame(
        IBTrACS, geometry=gpd.points_from_xy(IBTrACS.LON, IBTrACS.LAT), crs=4326
    )
    # Convert ISO_TIME from string to datetime
    IBTrACS["ISO_TIME"] = pd.to_datetime(
        IBTrACS["ISO_TIME"], format="%Y-%m-%d %H:%M:%S"
    )
    # Round to the nearest 3 hours to match MSWEP
    IBTrACS["ISO_TIME"] = IBTrACS["ISO_TIME"].dt.round("3H")

    # Get path to the modified rasters
    dict_nc_paths = pd.DataFrame(paths_to_rasters())
    # Round to the nearest 3 hours to match IBTrACS
    dict_nc_paths["time"] = dict_nc_paths["time"].round("3H")

    # Add information in IBTrACS
    IBTrACS = pd.merge(
        left=IBTrACS.reset_index(),
        right=dict_nc_paths,
        left_on=["ISO_TIME"],
        right_on=["time"],
    )
    # Correct the coordinates outside of the bounding box of -180,180, -90, 90.
    IBTrACS["LON"] = IBTrACS["LON"].apply(lambda x: np.sign(x) * (np.abs(x) % 180))
    # set index name as 'row_id'
    IBTrACS.rename(columns={"index": "row_id"}, inplace=True)

    return IBTrACS
