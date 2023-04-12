#!/usr/bin/env python3
# %%
# Packages ----------------------------------------------------------------
from pyprojroot.here import here

# Geo stuff:
import geopandas as gpd
import pandas as pd
import xarray as xr
import rioxarray
import numpy as np

import shutil
import pickle
from glob import glob
import os
from gc import collect  # garbage collector
from tqdm import tqdm  # progress bar, tqdm.pandas() to use progress_apply with pandas

# My scripts
from geodesic_point_buffer import geodesic_point_buffer  # geosesic buffer
from raster_statistics import raster_statistics
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
    # set index name as 'row_id'
    IBTrACS.index.set_names("row_id", inplace=True)
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
    return IBTrACS


IBTrACS = IBTrACS_read(path_IBTrACS="Data/Data_IBTrACS/ibtracs.ALL.list.v04r00.csv")

# Correct the longitude coordinate outside of the bounding box of -180,180, -90, 90.
IBTrACS["LON"] = IBTrACS["LON"].apply(lambda x: np.sign(x) * (np.abs(x) % 180))

# For this exercise, select Harvey
# IBTrACS = IBTrACS[IBTrACS.SID == "2017228N14314"]

# %% =============================================================================================================================
# Temporary folder file to back up the for loop progression
os.makedirs(here("Data/Data_TCP_Temporary/"), exist_ok=True)

# list where the extractions will be stored
Spatial_average_MSWEP = []
iter = 0  # used for save result at each 1000 iterations

for index, row in tqdm(IBTrACS.iterrows(), total=IBTrACS.shape[0]):
    iter += 1
    path_file = row["origin"]
    lon = row["LON"]
    lat = row["LAT"]

    # Import the corresponding MSWEP file
    MSWEP = (
        xr.open_dataarray(path_file)
        .rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
        .rio.write_crs("epsg:4326", inplace=True)
        # .chunk(chunks={"lat": 100, "lon": 100}) # hard to make this work
    )

    # Compute the statistics
    Spatial_average_MSWEP.append(
        raster_statistics(
            DataArray=MSWEP, lat=lat, lon=lon, max_radius_km=500, ring_thickness_km=10
        )
    )
    # Close MSWEP
    MSWEP.close()

    # At each 1000 iterations, create a serialized back-up of Spatial_average_MSWEP object.
    # The back-up will be deleted at the end of the script.
    if iter % 1000 == 0:
        path = here("Data/Data_TCP_Temporary/")
        # if False:
        # remove old back-up file
        for fname in os.listdir(path):
            if fname.startswith("Spatial_average_MSWEP_back-up"):
                os.remove(os.path.join(path, fname))

        # Make a new back-up file
        fBackUp = os.path.join(
            path, f"Spatial_average_MSWEP_back-up_iteration:{iter}.pickle"
        )
        with open(fBackUp, "wb") as handle:
            pickle.dump(Spatial_average_MSWEP, handle, protocol=pickle.HIGHEST_PROTOCOL)
        # Get garbage collector to empty memory
        collect()

# %% =============================================================================================================================
# Set index as row_id
if IBTrACS.index.name != "row_id":
    IBTrACS.set_index("row_id", inplace=True)

# Create the Dataframe containing the non binned statistics -----------------------
pdf_Spatial_average_MSWEP = (
    pd.DataFrame(Spatial_average_MSWEP)
    .set_index(IBTrACS.index)
    .drop("bin_statistics", axis=1)
    .apply(pd.Series.explode)
)
# Keep the necessary columns in IBTrACS and add the non binned raster statistics
col_IBTrACS = ["SID", "NAME", "ISO_TIME", "LAT", "LON"]
col_select_non_binned = col_IBTrACS + pdf_Spatial_average_MSWEP.columns.to_list()
IBTrACS_non_binned = IBTrACS.join(pdf_Spatial_average_MSWEP, on="row_id")
IBTrACS_non_binned = IBTrACS_non_binned[col_select_non_binned]

# Create dataframe containing the binned statistics --------------------------------
bin_statistics = (
    pd.DataFrame(x["bin_statistics"] for x in Spatial_average_MSWEP)
    .set_index(IBTrACS.index)
    .apply(pd.Series.explode)
)
# Keep the necessary columns in IBTrACS and add the binned raster statistics
col_select_binned = col_IBTrACS + bin_statistics.columns.to_list()
IBTrACS_binned = IBTrACS.join(bin_statistics, on="row_id")
IBTrACS_binned = IBTrACS_binned[col_select_binned]


# Create directory if absent and save result to file -------------------------------
os.makedirs(here("Data/Data_TCP/"), exist_ok=True)
IBTrACS_non_binned.to_csv(
    here("Data/Data_TCP/non_binned_statistics.csv.gz"), compression="gzip"
)
IBTrACS_binned.to_csv(
    here("Data/Data_TCP/binned_statistics.csv.gz"), compression="gzip"
)

# Remove the back-up directory (Data_TCP_Temporary) when the script is finished ----
path = here("Data/Data_TCP_Temporary")
try:
    shutil.rmtree(path)
except OSError as e:
    print("Error: %s - %s." % (e.filename, e.strerror))
# %%
