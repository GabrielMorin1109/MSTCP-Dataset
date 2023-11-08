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
ring_thickness_km = 100

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

    def buffers_geometry(self, max_radius_km):
        buffers = self.gdf.apply(
            lambda row: geodesic_point_buffer(
                lat=row["LAT"], lon=row["LON"], radius_km=max_radius_km
            ),
            axis=1,
        )
        return buffers

    def xr_maps_from_gpd(self):
        # Trajectory buffer
        geometries_gds = self.buffers_geometry(max_radius_km=max_radius_km)

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
        geometries_gds = self.buffers_geometry(max_radius_km=max_radius_km)
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
        # get the buffers as one geometry object
        geometries_gds = self.buffers_geometry(max_radius_km=max_radius_km)
        self_stats = self.xr_maps_func(dim=dim, np_func=np_func).rio.clip(
            geometries_gds.geometry.values, geometries_gds.crs
        )
        return self_stats

    def plot(
        self,
        np_func,
        legend_title,
        path_plot_to_file,
        levels_input,
        dim,
        stats_plot_b,
        ring_buffer_plot_b,
        rainfall_area_plot_b,
    ):
        from shapely.geometry import LineString, Point, Polygon
        from matplotlib import cm
        from matplotlib.colors import ListedColormap, LinearSegmentedColormap
        import cartopy.feature as cf

        # get the buffers as one geometry object
        geometries_gds = self.buffers_geometry(max_radius_km=max_radius_km)

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

        # Add precipitation data ~
        # Stats plot is different than other plots, a Boolean is used to switch the procedure.
        if stats_plot_b is True and self.gdf.geometry.shape[0] == 1:
            # cmap.set_bad("gray")
            precip_xar = self.xr_maps_from_gpd()
            if rainfall_area_plot_b is True:
                # ax.set_facecolor("gray")
                precip_xar = precip_xar.where(precip_xar >= 0.5)

            xarray_plot = precip_xar[0].plot.contourf(
                ax=ax,
                transform=ccrs.PlateCarree(),
                cbar_kwargs=cbar_kwargs,
                levels=levels_input,
                cmap=cmap,
                extend="max",
            )
            # xarray_plot = precip_xar.plot(
            #     ax=ax,
            #     transform=ccrs.PlateCarree(),
            #     cbar_kwargs=cbar_kwargs,
            #     levels=levels_input,
            #     cmap=cmap,
            #     extend="max",
            # )
            if ring_buffer_plot_b is True:
                # Add ring buffers
                tmp = gpd.GeoSeries()
                # tmp = self.buffers_geometry(max_radius_km=20)[0].difference(
                #     self.buffers_geometry(max_radius_km=10)[0]
                # )
                for bin in range(0, max_radius_km, ring_thickness_km):
                    tmp.loc[len(tmp)] = self.buffers_geometry(
                        max_radius_km=bin + ring_thickness_km
                    )[0]
                tmp.exterior.plot(
                    ax=ax,
                    transform=ccrs.PlateCarree(),
                    color="black",
                    alpha=0.6,
                    linewidth=0.4,
                )
                # tmp_to_remove = tmp.envelope.difference(tmp)
                # gpd.GeoDataFrame(geometry=[tmp_to_remove], crs="EPSG:4326").plot(
                #     ax=ax,
                #     transform=ccrs.PlateCarree(),
                #     color="gray",
                #     alpha=0.9,
                #     linewidth=0.4,
                # )

                # for bin in range(20, 40, ring_thickness_km):
                #     tmp.loc[len(tmp)] = self.buffers_geometry(
                #         max_radius_km=bin + ring_thickness_km
                #     )[0]
                # tmp.exterior.plot(
                #     ax=ax,
                #     transform=ccrs.PlateCarree(),
                #     color="black",
                #     alpha=0.6,
                #     linewidth=0.4,
                # )
            elif ring_buffer_plot_b is False:
                # Add background
                self.buffers_geometry(max_radius_km=max_radius_km).exterior.plot(
                    ax=ax,
                    transform=ccrs.PlateCarree(),
                    color="black",
                    alpha=1,
                )
                lon, lat = (
                    precip_xar[0]
                    .where(precip_xar[0] == precip_xar[0].max(), drop=True)
                    .indexes.values()
                )

                if rainfall_area_plot_b is False:
                    # TODO
                    # Add distance between the center and the maximum Precipitation value
                    lon_center, lat_center = self.gdf.geometry[0].xy
                    coords_center = [lon_center[0], lat_center[0]]
                    coords_max_precip = [lon[0], lat[0]]
                    dist_max = LineString([coords_center, coords_max_precip])
                    ax.plot(
                        dist_max.xy[0],
                        dist_max.xy[1],
                        color="black",
                        transform=ccrs.PlateCarree(),
                    )
                    # Add maximum precipitation value to the plot
                    ax.scatter(
                        lon[0],
                        lat[0],
                        marker="x",
                        color="black",
                        s=100,
                        transform=ccrs.PlateCarree(),
                    )
        elif stats_plot_b is False or stats_plot_b is None:
            if np_func == "nansum":
                xarray_plot = self.apply_fun(dim=dim, np_func=np_func) * 3
            else:
                xarray_plot = self.apply_fun(dim=dim, np_func=np_func)
            xarray_plot = xarray_plot.plot.contourf(
                ax=ax,
                transform=ccrs.PlateCarree(),
                cbar_kwargs=cbar_kwargs,
                levels=levels_input,
                cmap=cmap,
                extend="max",
            )

        # reduce text size of cbar
        cbar = xarray_plot.colorbar
        cbar.ax.tick_params(labelsize=8)

        # Add IBTrACS TC track
        if self.gdf.geometry.shape[0] > 1:
            geo_object = gpd.GeoSeries(LineString(self.gdf.geometry.values))
        elif self.gdf.geometry.shape[0] == 1:
            geo_object = self.gdf.geometry
        if rainfall_area_plot_b is False:
            geo_object.plot(
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
        # if ring_buffer_plot_b is True:
        #     minx, miny, maxx, maxy = tmp.bounds
        ax.set_extent([minx, maxx, miny, maxy], crs=crs)
        if stats_plot_b is True:
            cbar.remove()

        # Save plot
        if path_plot_to_file is not None:
            plt.savefig(path_plot_to_file)
        plt.show()

    ######################################################################
    def stats_plot(
        self,
        np_func,
        legend_title,
        path_plot_to_file,
        levels_input,
        dim,
    ):
        from shapely.geometry import LineString, Point, Polygon
        from matplotlib import cm
        from matplotlib.colors import ListedColormap, LinearSegmentedColormap
        import cartopy.feature as cf

        # get the buffers as one geometry object
        geometries_gds = self.buffers_geometry(max_radius_km=max_radius_km)

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
        precip_xar = self.xr_maps_from_gpd()
        xarray_plot = (
            precip_xar.where(precip_xar >= 1.5, other=0)
            .pipe(lambda x: x.where(x <= 1.5, other=1.5))
            .pipe(lambda x: x.where(x == 0, other=1e8))[0]
            .plot.contour(
                ax=ax,
                transform=ccrs.PlateCarree(),
                cmap=cmap,
                levels=[0, 1.5],
                # colors="black",
                extend="max",
                linewidths=1,
            )
        )
        #     ax=ax,
        #     transform=ccrs.PlateCarree(),
        #     # cbar_kwargs={
        #     #     "orientation": "horizontal",
        #     #     "shrink": 0.6,
        #     #     "pad": 0.05,
        #     #     "aspect": 40,
        #     # },
        #     # cmap=cmap,
        #     levels=[0, 1.5],
        #     colors="black",
        #     extend="max",
        #     linewidths=1,
        #     # levels=[1]
        # )
        ################################
        # Add ring buffers
        tmp = gpd.GeoSeries()
        for bin in range(0, max_radius_km, ring_thickness_km):
            tmp.loc[len(tmp)] = self.buffers_geometry(
                max_radius_km=bin + ring_thickness_km
            )[0]
        tmp.exterior.plot(
            ax=ax,
            transform=ccrs.PlateCarree(),
            color="black",
            alpha=0.6,
            linewidth=0.4,
        )

        # Add background
        # self.buffers_geometry(max_radius_km=max_radius_km).exterior.plot(
        #     ax=ax,
        #     transform=ccrs.PlateCarree(),
        #     color="black",
        #     alpha=1,
        # )
        lon, lat = (
            precip_xar[0]
            .where(precip_xar[0] == precip_xar[0].max(), drop=True)
            .indexes.values()
        )

        xarray_plot = self.apply_fun(dim=dim, np_func=np_func)
        xarray_plot = xarray_plot.plot.contourf(
            ax=ax,
            transform=ccrs.PlateCarree(),
            cbar_kwargs=cbar_kwargs,
            levels=levels_input,
            cmap=cmap,
            extend="max",
        )

        # reduce text size of cbar
        cbar = xarray_plot.colorbar
        cbar.ax.tick_params(labelsize=8)
        # Add IBTrACS TC track
        if self.gdf.geometry.shape[0] > 1:
            geo_object = gpd.GeoSeries(LineString(self.gdf.geometry.values))
        elif self.gdf.geometry.shape[0] == 1:
            geo_object = self.gdf.geometry

        geo_object.plot(
            ax=ax,
            transform=ccrs.PlateCarree(),
            color="black",
            alpha=1,
        )
        ################################
        # Add distance between the center and the maximum Precipitation value
        lon_center, lat_center = self.gdf.geometry[0].xy
        coords_center = [lon_center[0], lat_center[0]]
        coords_max_precip = [lon[0], lat[0]]
        dist_max = LineString([coords_center, coords_max_precip])
        ax.plot(
            dist_max.xy[0],
            dist_max.xy[1],
            color="black",
            transform=ccrs.PlateCarree(),
        )
        # Add maximum precipitation value to the plot
        ax.scatter(
            lon[0],
            lat[0],
            marker="x",
            color="black",
            s=100,
            transform=ccrs.PlateCarree(),
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
        plt.savefig(path_plot_to_file)
        plt.show()
        # return self.xr_maps_from_gpd().where(precip_xar >= 0.5)[0]


# %%
spatial_manipulations(IBTrACS.iloc[[0]]).stats_plot(
    np_func="nansum",
    dim=("time"),
    legend_title="Total precipitation",
    path_plot_to_file=here("plots/all_stats.png"),
    levels_input=[0, 0.5, 1, 1.5, 2, 2.5],
)
# %%
# Illustration of the statistics using one observation
spatial_manipulations(IBTrACS.iloc[[0]]).plot(
    np_func="nansum",
    dim=("time"),
    legend_title="Total precipitation",
    path_plot_to_file=here("plots/statistics_representation_non_binned.png"),
    levels_input=[0, 0.5, 1, 1.5, 2, 2.5],
    stats_plot_b=True,
    ring_buffer_plot_b=False,
    rainfall_area_plot_b=False,
)
# binned statistics
spatial_manipulations(IBTrACS.iloc[[0]]).plot(
    np_func="nansum",
    dim=("time"),
    legend_title="Total precipitation",
    path_plot_to_file=here("plots/statistics_representation_binned.png"),
    levels_input=[0, 0.5, 1, 1.5, 2, 2.5],
    stats_plot_b=True,
    ring_buffer_plot_b=True,
    rainfall_area_plot_b=False,
)
# rainfall area:
spatial_manipulations(IBTrACS.iloc[[0]]).plot(
    np_func="nansum",
    dim=("time"),
    legend_title="Total precipitation",
    path_plot_to_file=here("plots/statistics_representation_binned_rainfall_area.png"),
    levels_input=[0, 1.5],  # [0, 0.5, 1, 1.5, 2, 2.5],
    stats_plot_b=True,
    ring_buffer_plot_b=False,
    rainfall_area_plot_b=True,
)

# ----------------------------------------------------------------
# %%
# tropical storm status
TCS = IBTrACS.loc[IBTrACS.USA_WIND.gt(34).idxmax() :]

# First graph :
spatial_manipulations(TCS).plot(
    np_func="nansum",
    dim=("time"),
    legend_title="Total precipitation",
    path_plot_to_file=here("plots/USA_WIND_gt_34_TC_precipitation.png"),
    levels_input=[x * 3 for x in [0, 1, 5, 10, 20, 30, 50, 80, 130, 210, 300]],
    stats_plot_b=False,
    ring_buffer_plot_b=False,
    rainfall_area_plot_b=False,
)

# %%
world_gdf = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres")).to_crs(
    "EPSG:4326"
)
USA = world_gdf[world_gdf.name == "United States of America"].geometry.unary_union
# Overland
geometries_gds = spatial_manipulations(IBTrACS).buffers_geometry(
    max_radius_km=max_radius_km
)
OL = IBTrACS[geometries_gds.intersects(USA)]

# the plot:
spatial_manipulations(OL).plot(
    np_func="nanmax",
    dim=("time"),
    legend_title="Maximum precipitation",
    path_plot_to_file=here("plots/max_TC_precipitation_rates.png"),
    levels_input=[0, 1, 5, 10, 15, 20, 25, 30, 35],
    stats_plot_b=False,
    ring_buffer_plot_b=False,
    rainfall_area_plot_b=False,
)
# %%
spatial_manipulations(OL).plot(
    np_func="nansum",
    dim=("time"),
    legend_title="Total precipitation",
    path_plot_to_file=here("plots/USA_TC_precipitation.png"),
    levels_input=[x * 3 for x in [0, 1, 5, 10, 20, 30, 50, 80, 130, 210, 300]],
    stats_plot_b=False,
    ring_buffer_plot_b=False,
    rainfall_area_plot_b=False,
)
# %%
