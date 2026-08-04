# Coordinate Reference Systems

[Coordinate Reference System](https://en.wikipedia.com/wiki/Spatial_reference_system)s are simply descriptions of what some set of coordinates really mean in reference to some standard.

In a mathematical sense, coordinate reference systems can be thought of defining a _space_, with associated transformations from and to latitude-longitude space (plate-carree, long-lat, WGS84) which is the default CRS we assume.

## Geographic CRS

If a CRS is _geographic_, its coordinates describe longitude and latitude rather than a flat plane.  `GO.area` defaults to `AutoManifold()`: when Proj is available and recognizes the CRS, it uses an ellipsoid-aware `Geodesic` calculation and returns square metres.  Without Proj, geographic coordinates use a degree longitude/latitude `Spherical` calculation.

Use [`Planar`](@ref) explicitly when area in the geometry's native coordinate units is required.

## Projected CRS

Projected geometries, geometries with an unknown CRS trait, and geometries without a CRS use `Planar` calculations in native square units.  A projected map area represents physical area only for an equal-area projection.

## Ways to describe CRS

Completely separate from the _meaning_ of the CRS is the way you describe or define it.  There are a [dizzying array of ways](@ref crs-format-table) to do this, but two easy ones are Proj strings and Well Known Text.

The geographic community seems to be standardizing on [Well Known Text]() as the "best" CRS identifier.  This is quite verbose, but is unambiguous and easy enough to read once you get the hang of it.

To indicate the type of CRS definition you're using, you can wrap a string in its corresponding `GeoFormatTypes` type.

## [CRS format table](@id crs-format-table)
<!-- TODO: convert this to a Markdown table-->
- Proj-strings: a brief but powerful way to describe a set of known CRS + some transformations to them.  Really useful when plotting and interactively adjusting CRS.  See the Proj docs.
- EPSG codes: a short way to refer to a known coordinate system in the database of the European Petroleum Survey Group.  Example: `EPSG:4236`.
- ESRI codes: similar to EPSG codes, but referring to CRS known to ESRI instead.  Example: `ESRI:12345`
- ProjJSON: a more structured way to express Proj-strings using JSON.  
- KML: key-markup language, an XML extension, used in web feature services
- Mapinfo CoordSys: 
