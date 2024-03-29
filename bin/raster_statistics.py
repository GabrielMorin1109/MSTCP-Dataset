#!/usr/bin/env python3
from pyprojroot.here import here
from pyproj import CRS
import geopandas as gpd
import pandas as pd
import xarray as xr
import rioxarray
import numpy as np
from geodesic_point_buffer import geodesic_point_buffer  # geosesic buffer
from haversine import haversine

# Area weights:
grid_area_MSWEP = xr.open_dataarray(here("Data/grid_area_MSWEP.nc"))


# Rasterstats:
def raster_statistics(
    DataArray,
    lon,
    lat,
    max_radius_km: int = 500,
    ring_thickness_km: int = 10,
    grid_area: xr.DataArray = grid_area_MSWEP,
    crs: int = 4326,
) -> dict:
    # Non binned raster statistics --------------------------------
    ## Clip raster with geodesic buffer at coordinates
    geometry = geodesic_point_buffer(lat=lat, lon=lon, radius_km=max_radius_km)
    DataArray_clipped = DataArray.rio.clip([geometry], CRS.from_epsg(crs), drop=False)

    area_averaged_TCP = (
        DataArray_clipped.weighted(grid_area).mean(("lon", "lat")).values / 3
    )

    # Select element over threshold and sum the area over threshold of 0.5 mm/h (thus 0.5 mm/h * 3/3 = 1.5 mm/3h)
    over_threshold = ~np.isnan(DataArray_clipped.where(DataArray_clipped >= 0.5 * 3))
    RA_over_threshold = grid_area.where(over_threshold).sum().values / (
        1000 * 1000
    )  # conversion from m2 to km2

    # Maximum and location of maximum precipitation within a 500km radius buffer
    xr_max_precipitation = DataArray_clipped.isel(
        DataArray_clipped.argmax(dim=["lon", "lat"])
    )
    lon_max_precipitation = xr_max_precipitation["lon"].values
    lat_max_precipitation = xr_max_precipitation["lat"].values
    max_precipitation = (
        xr_max_precipitation.values / 3
    )  # Divided by 3 to convert from mm/3h to mm/h
    radius_of_maximum_rain = haversine(
        lon1=lon, lat1=lat, lon2=lon_max_precipitation, lat2=lat_max_precipitation
    )

    # Binned raster statistics ------------------------------------
    bin_sequence = [*range(0, max_radius_km, ring_thickness_km)]
    binned_area_averaged_TCP = []

    for radius_bin in bin_sequence:
        # Clip raster with the polygon
        inner_bin_geometry = geodesic_point_buffer(
            lat=lat, lon=lon, radius_km=radius_bin
        )
        outer_bin_geometry = geodesic_point_buffer(
            lat=lat, lon=lon, radius_km=radius_bin + ring_thickness_km
        )
        bin_geometry = outer_bin_geometry.difference(inner_bin_geometry)

        bin_DataArray_clipped = DataArray_clipped.rio.clip(
            [bin_geometry], CRS.from_epsg(crs), drop=False
        )
        # Add grid cell area as weight and calculate the weighted average.
        bin_DataArray_clipped_weighted = bin_DataArray_clipped.weighted(grid_area)

        binned_area_averaged_TCP.extend(
            bin_DataArray_clipped_weighted.mean(("lon", "lat")).values / 3
        )

    OUT = {
        "RA_over_threshold": RA_over_threshold,
        "radius_of_maximum_rain": radius_of_maximum_rain,
        "lon_max_precipitation": lon_max_precipitation,
        "lat_max_precipitation": lat_max_precipitation,
        "max_precipitation": max_precipitation,
        "area_averaged_TCP": area_averaged_TCP,
        "bin_sequence": bin_sequence,
        "binned_area_averaged_TCP": binned_area_averaged_TCP,
    }
    return OUT
