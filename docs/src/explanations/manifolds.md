# Manifolds

A manifold is, mathematically, a description of some space that is locally Euclidean (i.e., locally flat).  
All geographic projections, and the surface of the sphere and ellipsoid, fall under this category of space - 
and these are all the spaces that are relevant to geographic geometry.

## What manifolds are available?

GeometryOps has four [`Manifold`](@ref) types: [`Planar`](@ref), [`Spherical`](@ref), [`Geodesic`](@ref), and [`AutoManifold`](@ref).

- `Planar()` is, as the name suggests, a perfectly Cartesian, usually 2-dimensional, space.  The shortest path from one point to another is a straight line.
- `Spherical(; radius)` describes points on the surface of a sphere of a given radius.  
  The most convenient sphere for geometry processing is the unit sphere, but one can also use 
  the sphere of the Earth for e.g. projections.
- `Geodesic(; semimajor_axis, inv_flattening)` describes points on the surface of a flattened ellipsoid, 
  similar to the Earth.  The parameters describe the curvature and shape of the ellipsoid, and are equivalent 
  to the flags `+a` and `+f` in Proj's ellipsoid specification.  The default values are the values of the WGS84
  ellipsoid.

  For `Geodesic`, we need an `AbstractGeodesic` that can wrap representations from Proj.jl and SphericalGeodesics.jl.
- `AutoManifold()` selects a manifold from a geometry's CRS when an operation supports automatic selection.

The idea here is that the manifold describes how the geometry needs to be treated.  

## Why this is needed

The classical problem this is intended to solve is that in GIS, latitude and longitude coordinates 
are often treated as planar coordinates, when they in fact live on the sphere/ellipsoid, and must be 
treated as such.  For example, computing the area of the USA on the lat/long plane yields a result of `1116`,
which is plainly nonsensical.  

## How this is done

In order to avoid this, GeometryOps combines CRS traits with manifolds.

1. GeoInterface's `crstrait`, which describes the CRS type of a geometry.
2. GeometryOps's `Manifold` type, which defines the surface on which to perform operations.
3. Proj, when loaded, which recognizes CRS definitions, supplies ellipsoid parameters, and converts projected linear units.


`GO.area` uses `AutoManifold()` by default.  With Proj, a recognized geographic CRS selects an ellipsoid-aware `Geodesic` calculation in square metres.  Without Proj, a geographic CRS selects a degree longitude/latitude `Spherical` calculation.  Projected, unknown, and CRS-less geometries use `Planar` native square units.

Passing `Planar()` explicitly always computes in native coordinate units.  A map-plane area is a physical area only in an equal-area projection.

## Algorithms and manifolds

Algorithms define what operation is performed on the geometry; however, the choice of algorithm can also depend on the manifold.  L'Huilier's algorithm for the area of a polygon is not applicable to the plane, but is applicable to either the sphere or ellipsoid, for example.