# Coordinate Reference Systems

[Coordinate Reference System](https://en.wikipedia.com/Spatial_reference_system)s are simply descriptions of what some set of coordinates really mean in reference to some standard.

In a mathematical sense, coordinate reference systems can be thought of defining a _space_, with associated transformations from and to latitude-longitude space (plate-carree, long-lat, WGS84) which is the default CRS we assume.

## Geographic CRS

If a CRS is _geographic_, that means that it refers to coordinates on a sphere.  Such coordinates should ideally be handled using a spherical geometry library like Google's s2.  GeometryOps does not currently handle spherical geometry computations except in special cases ([`perimeter`](@ref), [`GeodesicSegments`](@ref) in `segmentize`, [`GeodesicDistance`](@ref)).

A non-geographic CRS is assumed to be in Cartesian space.

## Projected CRS

Projected CRS are generally treated as Cartesian.

## Ways to describe CRS

The geographic community seems to be standardizing on [Well Known Text]() as the "best" CRS identifier.  This is quite verbose, but is unambiguous and easy enough to read once you get the hang of it.

To indicate the type of CRS definition you're using, you can wrap a string in its corresponding `GeoFormatTypes` type.

## CRS format table
<!-- TODO: convert this to a Markdown table-->
- Proj-strings: a brief but powerful way to describe a set of known CRS + some transformations to them.  Really useful when plotting and interactively adjusting CRS.  See the Proj docs.
- EPSG codes: a short way to refer to a known coordinate system in the database of the European Petroleum Survey Group.  Example: `EPSG:4236`.
- ESRI codes: similar to EPSG codes, but referring to CRS known to ESRI instead.  Example: `ESRI:12345`
- ProjJSON
- KML
- Mapinfo CoordSys