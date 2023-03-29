#!/usr/bin/env python3
#%% Packages ----------------------------------------------------------------
from pyprojroot.here import here

# Geo stuff:
import geopandas as gpd
import pandas as pd
import dask_geopandas as dask_gpd  # parallel geopandas


#%%
def gdf_ring(gdf, npartitions=10, inner_radius=0, ring_thickness=10 * 1e3):
    """
    A ring buffer can be thought of as the shadows of a doughnut on a 2D surface.

    Goal:
    The function takes a Geopandas Dataframe (the gdf object) that contains the coordinates
    of the location where a circular ring buffer should be created.

    Parameters:
    The inner_radius parameter represents the radius of the empty hole.
    The ring_thickness parameter represents the thickness of the ring.

    # Output:
    The function returned the Geopandas Dataframe, but the ring buffer replace the existing coordinates.
    """
    # print message about projection units
    unit_CRS = gdf.crs.axis_info[1].unit_name
    print(f"Check if ring distance units are in {unit_CRS}")

    # Error checking
    if inner_radius < 0 or ring_thickness <= 0:
        raise ValueError(
            "Check if inner_radius is negative or if ring thickness is less or equal to zero."
        )

    # Convert geopandas to dask geopandas
    dask_gdf = dask_gpd.from_geopandas(gdf, npartitions=npartitions)

    # Create ring buffer in parallel:
    outer_radius = inner_radius + ring_thickness
    ring_geometry = (
        dask_gdf.buffer(outer_radius)
        .difference(dask_gdf.buffer(inner_radius), align=True)
        .compute()
    )

    # Create the data structure
    gdf_rings = gpd.GeoDataFrame(
        {"distance_from_eye_to_InnerBin": inner_radius, "geometry": ring_geometry}
    )

    # index names & values, and set to the output dataframe
    idx = pd.MultiIndex.from_product([[inner_radius], dask_gdf.index.values.compute()])

    # MultiIndex to store the distance values and row indices
    idx.set_names(["distance_from_eye_to_InnerBin", "row_id"], inplace=True)
    gdf_rings = gdf_rings.set_index(idx)

    # Merge key information to the ring buffer:
    # gdf_rings = gdf_rings.join(gdf.loc[:, ("SID", "ISO_TIME")])
    gdf_rings = gdf_rings.join(gdf.drop(columns="geometry"))

    # return back to 4326
    gdf_rings = gdf_rings.to_crs(4326)

    # Print information
    return gdf_rings
