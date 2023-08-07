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
import re

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
from IBTrACS_read import IBTrACS_read

# %% =============================================================================================================================
IBTrACS = IBTrACS_read(path_IBTrACS="Data/Data_IBTrACS/ibtracs.ALL.list.v04r00.csv")

# %% =============================================================================================================================
# Create a Backup folder to back up the for loop progression
Backup_path = here("Data/Data_TCP_Backup/")
os.makedirs(Backup_path, exist_ok=True)

# %%
# If Backup_path has file, then read it and start the loop at the last saved iteration:
Backup_file_name = os.listdir(Backup_path)
# test if there is a file in the Backup folder
if Backup_file_name:
    Spatial_average_MSWEP = []
    # full path to file in the Backup folder
    Backup_file_path = [os.path.join(Backup_path, f) for f in Backup_file_name][0]
    with open(Backup_file_path, "rb") as openfile:
        while True:
            try:
                Spatial_average_MSWEP.append(
                    pickle.load(openfile)
                )  # should have used extent...
            except EOFError:
                break
    # Correction of extending the list, since I used append...
    Spatial_average_MSWEP = [Spatial_average_MSWEP[0][0] + Spatial_average_MSWEP[0][1:]]
    iter_last_backup = len(Spatial_average_MSWEP[0])

    time_IBTrACS = [x["time"] for x in Spatial_average_MSWEP[0]]
    SID_IBTrACS = [x["SID"] for x in Spatial_average_MSWEP[0]]

    iterated_rows = pd.MultiIndex.from_arrays(
        [SID_IBTrACS, time_IBTrACS], names=("SID", "ISO_TIME")
    )
    all_rows_to_iterate = pd.MultiIndex.from_frame(IBTrACS.loc[:, ["SID", "ISO_TIME"]])
    rows_to_iterate = ~all_rows_to_iterate.isin(iterated_rows)

    # rows_to_iterate = ~(
    #     IBTrACS.ISO_TIME.isin(time_IBTrACS) & IBTrACS.SID.isin(SID_IBTrACS)
    # )

else:
    # list where the extractions will be stored
    Spatial_average_MSWEP = []
    iter_last_backup = 0
    rows_to_iterate = pd.Series(
        True for i in range(IBTrACS.shape[0])
    )  # select all rows


# %%
# Back-up the loop progression after an iteration into a Backup folder
def back_up(iter_value):
    path = here("Data/Data_TCP_Backup/")
    # remove old back-up file
    for fname in os.listdir(path):
        if fname.startswith("Spatial_average_MSWEP_back-up"):
            os.remove(os.path.join(path, fname))

    # Make a new back-up file
    fBackUp = os.path.join(
        path, f"Spatial_average_MSWEP_back-up_iteration:{iter_value}.pickle"
    )
    with open(fBackUp, "wb") as handle:
        pickle.dump(Spatial_average_MSWEP, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # Get garbage collector to empty memory
    collect()


# %%
if np.any(rows_to_iterate):
    iter = iter_last_backup  # used for save result at each 1000 iterations
    for index, row in tqdm(
        IBTrACS[rows_to_iterate].iterrows(),
        total=IBTrACS[rows_to_iterate].shape[0],
        position=0,
        leave=True,
        initial=iter_last_backup,
    ):
        iter += 1
        path_file = row["origin"]
        lon = row["LON"]
        lat = row["LAT"]
        # When there is more than one variable in the nc file,
        # then select only precipitation
        try:
            MSWEP = xr.open_dataarray(path_file)
        except ValueError:
            MSWEP = xr.open_dataset(path_file)["precipitation"]
        MSWEP.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
        MSWEP.rio.write_crs("epsg:4326", inplace=True)
        # MSWEP.chunk(chunks={"lat": 100, "lon": 100}) # hard to make this work

        # Compute the statistics
        Raster_analysis = raster_statistics(
            DataArray=MSWEP, lat=lat, lon=lon, max_radius_km=500, ring_thickness_km=10
        )
        Raster_analysis["time"] = row["time"]
        Raster_analysis["SID"] = row["SID"]
        Raster_analysis["row_id"] = row["row_id"]

        Spatial_average_MSWEP.append(Raster_analysis)
        # Close MSWEP
        MSWEP.close()

        # At each 1000 iterations, create a serialized back-up of Spatial_average_MSWEP object.
        # The back-up will be deleted at the end of the script.
        if iter % 1000 == 0:
            back_up(iter_value=iter)

    # When the loop finishes, back up one more time
    Spatial_average_MSWEP = (
        Spatial_average_MSWEP[0] + Spatial_average_MSWEP[1:]
    )  # again, should have use extent, correct the data format
    back_up(iter_value=iter)


# %% =============================================================================================================================
# Set index as row_id
if IBTrACS.index.name != "row_id":
    IBTrACS.set_index("row_id", inplace=True)
if len(Spatial_average_MSWEP) == 1:
    Spatial_average_MSWEP = Spatial_average_MSWEP[0]

# Create the Dataframe containing the non binned statistics -----------------------
import copy

pdf_Spatial_average_MSWEP = copy.deepcopy(
    pd.DataFrame.from_dict(Spatial_average_MSWEP, orient="columns")
    .set_index("row_id")
    .drop("bin_statistics", axis=1)
    .apply(pd.Series.explode)
    .rename(columns={"time": "ISO_TIME"})
)
# Keep the necessary columns in IBTrACS and add the non binned raster statistics
col_IBTrACS = ["SID", "NAME", "ISO_TIME", "LAT", "LON"]
col_select_non_binned = np.unique(
    col_IBTrACS + pdf_Spatial_average_MSWEP.columns.to_list()
)
IBTrACS_non_binned = IBTrACS.join(
    pdf_Spatial_average_MSWEP,
    on="row_id",
    how="outer",
    rsuffix="_right",
)
IBTrACS_non_binned = IBTrACS_non_binned[col_select_non_binned]

# Create dataframe containing the binned statistics --------------------------------
bin_statistics = copy.deepcopy(
    pd.DataFrame.from_dict(Spatial_average_MSWEP, orient="columns")
    .loc[:, ["row_id", "SID", "time", "bin_statistics"]]
    .rename(columns={"time": "ISO_TIME"})
)
bin_statistics = pd.concat(
    [
        bin_statistics.drop(["bin_statistics"], axis=1),
        pd.DataFrame(bin_statistics.bin_statistics.to_list()).apply(pd.Series.explode),
    ],
    axis=1,
)

# Keep the necessary columns in IBTrACS and add the binned raster statistics
col_select_binned = np.unique(col_IBTrACS + bin_statistics.columns.to_list())

IBTrACS_binned = pd.merge(
    IBTrACS, bin_statistics, on="row_id", how="outer", suffixes=("", "_right")
)
IBTrACS_binned = IBTrACS_binned[col_select_binned]

# Create directory if absent and save result to file -------------------------------
os.makedirs(here("Data/Data_TCP/"), exist_ok=True)
IBTrACS_non_binned.to_csv(
    here("Data/Data_TCP/non_binned_statistics.csv.gz"), compression="gzip"
)
IBTrACS_binned.to_csv(
    here("Data/Data_TCP/binned_statistics.csv.gz"), compression="gzip"
)

# Remove the back-up directory (Data_TCP_Backup) when the script is finished ----
if False:
    path = here("Data/Data_TCP_Backup")
    try:
        shutil.rmtree(path)
    except OSError as e:
        print("Error: %s - %s." % (e.filename, e.strerror))
# %%
