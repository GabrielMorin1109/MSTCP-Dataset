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
import json
import ast

import shutil
import pickle
from glob import glob
import os
from gc import collect  # garbage collector
from tqdm import tqdm


from geodesic_point_buffer import geodesic_point_buffer  # geosesic buffer
from raster_statistics import raster_statistics
from get_paths_to_rasters import paths_to_rasters
from IBTrACS_read import IBTrACS_read

# %% =============================================================================================================================
IBTrACS = IBTrACS_read(path_IBTrACS="Data/Data_IBTrACS/ibtracs.ALL.list.v04r00.csv")

# %% =============================================================================================================================
computed_extractions = [
    os.path.splitext(computed_extraction)[0]
    for computed_extraction in os.listdir(here("Data/Data_TCP/json/"))
]

IBTrACS = IBTrACS[~IBTrACS.row_id.isin([computed_extractions])]

# %% =============================================================================================================================
for index, row in tqdm(
    IBTrACS.iterrows(),
    total=IBTrACS.shape[0],
    position=0,
    leave=True,
):
    path_file = row["origin"]
    variant = row["variant"]
    lon = row["LON"]
    lat = row["LAT"]
    row_id = row["row_id"]
    # When there is more than one variable in the nc file,
    # then select only precipitation
    try:
        MSWEP = xr.open_dataarray(path_file)
    except ValueError:
        MSWEP = xr.open_dataset(path_file)["precipitation"]
    MSWEP.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
    MSWEP.rio.write_crs("epsg:4326", inplace=True)

    # Compute the statistics
    Raster_analysis = raster_statistics(
        DataArray=MSWEP, lat=lat, lon=lon, max_radius_km=500, ring_thickness_km=10
    )
    Raster_analysis["time"] = row["time"]
    Raster_analysis["SID"] = row["SID"]
    Raster_analysis["variant"] = row["variant"]
    Raster_analysis["row_id"] = row_id
    Raster_analysis_path = here(f"Data/Data_TCP/{row_id}.csv")
    pd.DataFrame.from_dict(Raster_analysis, orient="index").transpose().to_csv(
        Raster_analysis_path
    )
    # Close MSWEP
    MSWEP.close()

# %% =============================================================================================================================
# combine statistics
mypath = str(here("Data/Data_TCP/"))
all_files = [
    mypath + f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))
]
combined_csv = pd.concat([pd.read_csv(f) for f in all_files])
combined_csv.to_csv(str(here("Data/Data_TCP.csv")))


# %%
# return to the old files format
df_TCP = pd.read_csv(str(here("Data/Data_TCP.csv")))
df_TCP.rename(columns={"time": "ISO_TIME"}, inplace=True)
df_TCP = df_TCP.iloc[:, 2:]
df_TCP = pd.merge(
    df_TCP,
    IBTrACS.loc[:, ["row_id", "LAT", "LON", "NAME"]],
    on="row_id",
    how="left",
    suffixes=("", "_right"),
)
cols_str_to_list = [
    "radius_of_maximum_rain",
    "lon_max_precipitation",
    "lat_max_precipitation",
    "max_precipitation",
    "area_averaged_TCP",
    "bin_sequence",
    "binned_area_averaged_TCP",
]

for column in cols_str_to_list:
    df_TCP[column] = df_TCP[column].apply(lambda x: ast.literal_eval(x))

# binned
new_binned = df_TCP.loc[
    :,
    [
        "SID",
        "ISO_TIME",
        "NAME",
        "LAT",
        "LON",
        "variant",
        "bin_sequence",
        "binned_area_averaged_TCP",
        "row_id",
    ],
].explode(["bin_sequence", "binned_area_averaged_TCP"], ignore_index=True)
new_binned.to_csv(str(here("Data/binned_Data_TCP.csv")), index=False)
# non_binned
non_binned_cols = [
    "row_id",
    "SID",
    "ISO_TIME",
    "NAME",
    "LAT",
    "LON",
    "variant",
    "RA_over_threshold",
    "area_averaged_TCP",
    "lat_max_precipitation",
    "lon_max_precipitation",
    "max_precipitation",
    "radius_of_maximum_rain",
]
cols_list_to_first = [
    "radius_of_maximum_rain",
    "lon_max_precipitation",
    "lat_max_precipitation",
    "max_precipitation",
    "area_averaged_TCP",
]
for column in cols_list_to_first:
    df_TCP[column] = df_TCP[column].apply(lambda x: x[0])

df_TCP.loc[:, non_binned_cols].to_csv(
    str(here("Data/non_binned_Data_TCP.csv")), index=False
)

# %%
