# # Extent embedding

"""
    embed_extent(obj)

Recursively wrap the object with a GeoInterface.jl geometry,
calculating and adding an `Extents.Extent` to all objects.

This can improve performance when extents need to be checked multiple times,
such when needing to check if many points are in geometries, and using their extents
as a quick filter for obviously exterior points.

# Keywords

$THREADED_KEYWORD
$CRS_KEYWORD
"""
embed_extent(x; threaded=false, crs=nothing) = 
    apply(identity, GI.PointTrait(), x; calc_extent=true, threaded, crs)
