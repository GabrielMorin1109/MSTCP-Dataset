#!/usr/bin/env python3
# %%
from pyproj import CRS, Transformer
from shapely.geometry import Point, Polygon, box as Box
from shapely.ops import transform
import numpy as np
import warnings
from shapely.errors import ShapelyDeprecationWarning


def geodesic_point_buffer(lat, lon, radius_km: int) -> Polygon:
    """
    Used a dynamic azimuthal equidistant projection to do a geodesic buffer.
    Script adapted from:
    - https://gis.stackexchange.com/a/289923
    - https://stackoverflow.com/a/38707971/13205929
    """
    # Azimuthal equidistant projection
    aeqd_proj = CRS.from_proj4(f"+proj=aeqd +lat_0={lat} +lon_0={lon} +x_0=0 +y_0=0")
    tfmr = Transformer.from_proj(aeqd_proj, aeqd_proj.geodetic_crs)
    buf = Point(0, 0).buffer(radius_km * 1000).buffer(0)  # distance in metres
    poly = transform(tfmr.transform, buf).exterior.coords[:]
    disk = Polygon(poly)
    # If the disk is without coordinate then do not go further
    if radius_km == 0:
        return disk

    # Fix up segments that cross the coordinate singularity at longitude ±180.
    # We do this unconditionally because it may or may not create a non-simple
    # polygon, depending on where the initial point was.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
        boundary = np.array(disk.boundary)
    i = 0
    while i < boundary.shape[0] - 1:
        if abs(boundary[i + 1, 0] - boundary[i, 0]) > 180:
            assert (boundary[i, 1] > 0) == (boundary[i, 1] > 0)
            vsign = -1 if boundary[i, 1] < 0 else 1
            hsign = -1 if boundary[i, 0] < 0 else 1
            boundary = np.insert(
                boundary,
                i + 1,
                [
                    [hsign * 180, boundary[i, 1]],
                    [hsign * 180, vsign * 90],
                    [-hsign * 180, vsign * 90],
                    [-hsign * 180, boundary[i + 1, 1]],
                ],
                axis=0,
            )
            i += 5
        else:
            i += 1
    disk = Polygon(boundary).buffer(0)

    # If the fixed-up polygon doesn't contain the origin point, invert it.
    if not disk.contains(Point(lon, lat)):
        disk = Box(-180, -90, 180, 90).difference(disk)

    assert disk.is_valid
    assert disk.boundary.is_simple
    assert disk.contains(Point(lon, lat))
    return disk


# %%
# import pyproj
# from shapely.geometry import Point, Polygon, box as Box
# from shapely.ops import transform as sh_transform
# from functools import partial
# import numpy as np


# def geodesic_point_buffer(lat, lon, radius_km):
#     """
#     Used a dynamic azimuthal equidistant projection to do a geodesic buffer.

#     Generate a shapely.Polygon object representing a disk on the
#     surface of the Earth, containing all points within RADIUS meters
#     of latitude/longitude LAT/LON.

#     Script adapted from: https://stackoverflow.com/a/38707971/13205929
#     """

#     wgs84_globe = pyproj.Proj(proj="latlong", ellps="WGS84")

#     aeqd = pyproj.Proj(proj="aeqd", ellps="WGS84", datum="WGS84", lat_0=lat, lon_0=lon)
#     disk = sh_transform(
#         partial(pyproj.transform, aeqd, wgs84_globe),
#         Point(0, 0).buffer(radius_km * 1000),
#     )

#     # Fix up segments that cross the coordinate singularity at longitude ±180.
#     # We do this unconditionally because it may or may not create a non-simple
#     # polygon, depending on where the initial point was.
#     boundary = np.array(disk.boundary)
#     i = 0
#     while i < boundary.shape[0] - 1:
#         if abs(boundary[i + 1, 0] - boundary[i, 0]) > 180:
#             assert (boundary[i, 1] > 0) == (boundary[i, 1] > 0)
#             vsign = -1 if boundary[i, 1] < 0 else 1
#             hsign = -1 if boundary[i, 0] < 0 else 1
#             boundary = np.insert(
#                 boundary,
#                 i + 1,
#                 [
#                     [hsign * 180, boundary[i, 1]],
#                     [hsign * 180, vsign * 90],
#                     [-hsign * 180, vsign * 90],
#                     [-hsign * 180, boundary[i + 1, 1]],
#                 ],
#                 axis=0,
#             )
#             i += 5
#         else:
#             i += 1
#     disk = Polygon(boundary).buffer(0)

#     # If the fixed-up polygon doesn't contain the origin point, invert it.
#     if not disk.contains(Point(lon, lat)):
#         disk = Box(-180, -90, 180, 90).difference(disk)

#     assert disk.is_valid
#     assert disk.boundary.is_simple
#     assert disk.contains(Point(lon, lat))
#     return disk
