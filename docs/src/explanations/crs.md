# Coordinate Reference Systems

[Coordinate Reference System](https://en.wikipedia.com/Spatial_reference_system)s are simply descriptions of what some set of coordinates really mean in reference to some standard.

In a mathematical sense, coordinate reference systems can be thought of defining a _space_, with associated transformations from and to latitude-longitude space (plate-carree, long-lat, WGS84) which is the default CRS we assume.

## Geographic CRS

If a CRS is _geographic_, that means that it refers to coordinates on a sphere.  Such coordinates should ideally be handled using a spherical geometry library like Google's s2.  GeometryOps does not currently handle spherical geometry computations except in special cases ([`perimeter`](@ref), [`GeodesicSegments`](@ref) in `segmentize`, [`GeodesicDistance`](@ref)).

A non-geographic CRS is assumed to be in Cartesian space.

## Projected CRS

