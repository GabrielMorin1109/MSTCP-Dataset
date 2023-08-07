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

import matplotlib.pyplot as plt

import seaborn as sns
import cartopy.crs as ccrs

# My scripts
from geodesic_point_buffer import geodesic_point_buffer  # geosesic buffer
from raster_statistics import raster_statistics
from get_paths_to_rasters import paths_to_rasters
from IBTrACS_read import IBTrACS_read

# Set radius of the geodesic buffer
max_radius_km = 500

# %%
IBTrACS = IBTrACS_read(path_IBTrACS="Data/Data_IBTrACS/ibtracs.ALL.list.v04r00.csv")

# Select Harvey:
IBTrACS = IBTrACS[IBTrACS.SID == "2017228N14314"]
IBTrACS.USA_WIND = IBTrACS.USA_WIND.astype(int)

# Some track wa recorded by multiple station. Keep the first occurrence
IBTrACS = IBTrACS[~IBTrACS.ISO_TIME.duplicated(keep="first")]

IBTrACS["np_ISO_TIME"] = IBTrACS.ISO_TIME.to_numpy()
IBTrACS.reset_index(inplace=True)


# %%
class spatial_manipulations:
    def __init__(self, gdf):
        self.gdf = gdf

    def buffers_geometry(self):
        buffers = self.gdf.apply(
            lambda row: geodesic_point_buffer(
                lat=row["LAT"], lon=row["LON"], radius_km=max_radius_km
            ),
            axis=1,
        )
        return buffers

    def xr_maps_from_gpd(self):
        # Trajectory buffer
        geometries_gds = self.buffers_geometry()

        # coordinate of the bounding box surrounding Harvey track
        minx, miny, maxx, maxy = geometries_gds.unary_union.bounds

        maps_xr = (
            xr.open_mfdataset(self.gdf.origin)
            .rio.set_spatial_dims(x_dim="lon", y_dim="lat")
            .rio.write_crs("epsg:4326")
            # Clip map to the bounding box
            .rio.clip_box(
                minx=minx,
                miny=miny,
                maxx=maxx,
                maxy=maxy,
                crs="EPSG:4326",
            )
        )
        maps_xr = (maps_xr.precipitation) / 3
        return maps_xr

    def map_clip_by_time(self):
        maps_xr = self.xr_maps_from_gpd()
        geometries_gds = self.buffers_geometry()
        OUT = []
        for idx in self.gdf.index:
            time_i = self.gdf.np_ISO_TIME.loc[idx]
            clipped_maps_xr = maps_xr.sel(time=time_i).rio.clip(
                geometries_gds.loc[idx : (idx + 1)].geometry.values
            )
            OUT.append(clipped_maps_xr)
        clipped_maps_xr = (
            xr.concat(OUT, "time")
            .rio.set_spatial_dims(x_dim="lon", y_dim="lat")
            .rio.write_crs("epsg:4326")
        )
        return clipped_maps_xr

    def xr_maps_func(self, np_func, dim: tuple = ()):
        func = getattr(np, np_func)
        maps_xr = self.map_clip_by_time().reduce(func, dim=dim)
        return maps_xr

    def apply_fun(self, np_func, dim: tuple = ()):
        self_stats = self.xr_maps_func(dim=dim, np_func=np_func).rio.clip(
            geometries_gds.geometry.values, geometries_gds.crs
        )
        return self_stats

    def plot(self, np_func, legend_title, path_plot_to_file, levels_input, dim):
        from shapely.geometry import LineString
        from matplotlib import cm
        from matplotlib.colors import ListedColormap, LinearSegmentedColormap
        import cartopy.feature as cf

        # get the buffers as one geometry object
        geometries_gds = self.buffers_geometry()

        projection = ccrs.Mercator()
        crs = ccrs.PlateCarree()

        # color:
        top = cm.get_cmap("PuBu")
        bottom = cm.get_cmap("gist_heat_r", 128)  # combine it all
        newcolors = np.vstack(
            (top(np.linspace(0, 1, 128)), bottom(np.linspace(0, 1, 128)[1:]))
        )
        cmap = ListedColormap(newcolors, name="OrangeBlue")

        plt.figure(dpi=150)
        ax = plt.axes(projection=projection, frameon=True)

        # Add land borders
        ax.add_feature(cf.COASTLINE.with_scale("50m"), lw=0.5, color="gray")
        ax.add_feature(cf.BORDERS.with_scale("50m"), lw=0.3, color="gray")
        cbar_kwargs = {
            "orientation": "horizontal",
            "shrink": 0.6,
            "pad": 0.05,
            "aspect": 40,
            "label": legend_title,
            "ticks": levels_input,
        }

        # Add precipitation data
        contourf_plot = self.apply_fun(dim=dim, np_func=np_func).plot.contourf(
            ax=ax,
            transform=ccrs.PlateCarree(),
            cbar_kwargs=cbar_kwargs,
            levels=levels_input,
            cmap=cmap,
            extend="max",
        )

        # reduce text size of cbar
        cbar = contourf_plot.colorbar
        cbar.ax.tick_params(labelsize=8)

        # Add IBTrACS TC track
        gpd.GeoSeries(LineString(self.gdf.geometry.values)).plot(
            ax=ax,
            transform=ccrs.PlateCarree(),
            color="black",
            alpha=1,
        )

        # Other customization
        gl = ax.gridlines(
            crs=crs,
            draw_labels=True,
            linewidth=0.6,
            color="gray",
            alpha=0.3,
            linestyle="-.",
        )
        gl.xlabel_style = {"size": 7}
        gl.ylabel_style = {"size": 7}
        ax.set_title("")

        # coordinate of the bounding box surrounding Harvey track
        minx, miny, maxx, maxy = geometries_gds.unary_union.bounds
        ax.set_extent([minx, maxx, miny, maxy], crs=crs)

        # Save plot
        if path_plot_to_file is not None:
            plt.savefig(path_plot_to_file)
        plt.show()


# %%
# tropical storm status
TCS = IBTrACS.loc[IBTrACS.USA_WIND.gt(34).idxmax() :]

# First graph :
spatial_manipulations(TCS).plot(
    np_func="nansum",
    dim=("time"),
    legend_title="Total precipitation",
    path_plot_to_file=here("plots/USA_WIND_gt_34_TC_precipitation.png"),
    levels_input=[0, 1, 5, 10, 20, 30, 50, 80, 130, 210, 300],
)

# %%
world_gdf = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres")).to_crs(
    "EPSG:4326"
)
USA = world_gdf[world_gdf.name == "United States of America"].geometry.unary_union
# Overland
geometries_gds = spatial_manipulations(IBTrACS).buffers_geometry()
OL = IBTrACS[geometries_gds.intersects(USA)]

# the plot:
spatial_manipulations(OL).plot(
    np_func="nanmax",
    dim=("time"),
    legend_title="Maximum precipitation",
    path_plot_to_file=here("plots/max_TC_precipitation_rates.png"),
    levels_input=[0, 1, 5, 10, 15, 20, 25, 30, 35],
)
# %%
spatial_manipulations(OL).plot(
    np_func="nansum",
    dim=("time"),
    legend_title="Total precipitation",
    path_plot_to_file=here("plots/USA_TC_precipitation.png"),
    levels_input=[0, 1, 5, 10, 20, 30, 50, 80, 130, 210, 300],
)

# %%
