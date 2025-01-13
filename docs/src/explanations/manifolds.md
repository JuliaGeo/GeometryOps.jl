# Manifolds

A manifold is, mathematically, a description of some space that is locally Euclidean (i.e., locally flat).  
All geographic projections, and the surface of the sphere and ellipsoid, fall under this category of space - 
and these are all the spaces that are relevant to geographic geometry.

## What manifolds are available?

GeometryOps has three [`Manifold`](@ref) types: [`Planar`](@ref), [`Spherical`](@ref), and [`Geodesic`](@ref).

- `Planar()` is, as the name suggests, a perfectly Cartesian, usually 2-dimensional, space.  The shortest path from one point to another is a straight line.
- `Spherical(; radius)` describes points on the surface of a sphere of a given radius.  
  The most convenient sphere for geometry processing is the unit sphere, but one can also use 
  the sphere of the Earth for e.g. projections.
- `Geodesic(; semimajor_axis, inv_flattening)` describes points on the surface of a flattened ellipsoid, 
  similar to the Earth.  The parameters describe the curvature and shape of the ellipsoid, and are equivalent 
  to the flags `+a` and `+f` in Proj's ellipsoid specification.  The default values are the values of the WGS84
  ellipsoid.

  For `Geodesic`, we need an `AbstractGeodesic` that can wrap representations from Proj.jl and SphericalGeodesics.jl.

The idea here is that the manifold describes how the geometry needs to be treated.  

## Why this is needed

The classical problem this is intended to solve is that in GIS, latitude and longitude coordinates 
are often treated as planar coordinates, when they in fact live on the sphere/ellipsoid, and must be 
treated as such.  For example, computing the area of the USA on the lat/long plane yields a result of `1116`,
which is plainly nonsensical.  

## How this is done

In order to avoid this, we've introduced three complementary CRS-related systems to the JuliaGeo ecosystem.  

1. GeoInterface's `crstrait`.  This is a method that returns the ideal CRS _type_ of a geometry, either Cartesian or Geographic.
2. Proj's `PreparedCRS` type, which extracts ellipsoid parameters and the nature of the projection from a coordinate reference system, and
   caches the results in a struct.  This allows GeometryOps to quickly determine the correct manifold to use for a given geometry.
3. GeometryOps's `Manifold` type, which defines the surface on which to perform operations.  This is what allows GeometryOps to perform
   calculations correctly depending on the nature of the geometry.


The way this flow works, is that when you load a geometry using GeoDataFrames, its CRS is extracted and parsed into a `PreparedCRS` type.
This is then used to determine the manifold to use for the geometry, and the geometry is converted to the manifold's coordinate system.

There is a table of known geographic coordinate systems in GeoFormatTypes.jl, and anything else is assumed to be 
a Cartesian or planar coordinate system.  CRStrait is used as the cheap determinant, but PreparedCRS is more general and better to use if possible.

When GeometryOps sees a geometry, it first checks its CRS to see if it is a geographic coordinate system.  If it is, it uses the `PreparedCRS`, or falls back to `crstrait` and geographic defaults to determine the manifold.

## Algorithms and manifolds

Algorithms define what operation is performed on the geometry; however, the choice of algorithm can also depend on the manifold.  L'Huilier's algorithm for the area of a polygon is not applicable to the plane, but is applicable to either the sphere or ellipsoid, for example.