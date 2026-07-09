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

## Extents and spatial indexing

The manifold also changes what a bounding box is.  On `Planar()`, the extent of a geometry
is the coordinate-wise minimum and maximum of its vertices — what `GI.extent` returns — and
it contains the whole geometry, because straight edges stay inside the box around their
endpoints.  On the sphere, the lon/lat box of the vertices under-covers in two ways:

- Edges are great-circle arcs, which bulge away from the chord between their endpoints.
  An arc between two points at latitude 60° runs poleward of the 60° parallel, leaving
  the box around its vertices.
- A ring can enclose a pole without touching it.  A cell whose vertices all sit at
  latitude 80° contains the pole and every longitude, and no union of per-edge boxes
  will discover that.

(Longitude also wraps at the antimeridian, so for a geometry crossing it the lon/lat
box isn't even well defined without splitting or special-casing.)

`Extents.extent(m::Manifold, geom)` computes extents *on the manifold*.  On `Planar()`
it is `GI.extent(geom)`.  On `Spherical()` it returns a 3D Cartesian
`Extent{(:X, :Y, :Z)}`: the box in ℝ³ around the geometry as a region on the unit
sphere.  Cartesian boxes have no antimeridian or pole singularities, and compose with
the ordinary `Extents.union` and `Extents.intersects`.  The box covers arc bulge and
enclosed poles; rings and polygons follow S2's convention of counterclockwise winding
with the interior on the left, so a clockwise ring is read as enclosing the (huge)
region outside it.

```@example manifold-extents
import GeometryOps as GO, GeoInterface as GI
import Extents

cap = GI.Polygon([[(lon, 80.0) for lon in 0.0:30.0:360.0]])  # a cell around the north pole

GI.extent(cap)   # the lon/lat box of the vertices — can't tell the pole is inside
```

```@example manifold-extents
GO.extent(GO.Spherical(), cap)   # the 3D box on the unit sphere — Z reaches 1
```

Spatial index construction accepts a manifold the same way: `RTree(m::Manifold, algorithm, geoms)`
and `NaturalIndex(m::Manifold, geoms)` build the tree over each geometry's extent on `m`.
On `Spherical()` the leaf boxes are the 3D boxes above, and queries work in that space —
with a 3D extent, or with any predicate over extents, like a `SphericalCap` through
`Extents.intersects`:

```@example manifold-extents
import GeometryOps.FlexibleRTrees: RTree, HPR
import GeometryOps.SpatialTreeInterface as STI
using GeometryOps.UnitSpherical: SphericalCap, UnitSphericalPoint

# 12 cells in a band from latitude 60° to 80°, plus the polar cap above them
band = [GI.Polygon([[(lon, 60.0), (lon + 30.0, 60.0), (lon + 30.0, 80.0), (lon, 80.0), (lon, 60.0)]])
        for lon in 0.0:30.0:330.0]
tree = RTree(GO.Spherical(), HPR(), vcat(band, [cap]))

polecap = SphericalCap(UnitSphericalPoint(0.0, 0.0, 1.0), 0.05)  # 0.05 rad around the pole
STI.query(tree, Base.Fix1(Extents.intersects, polecap))
```

Only geometry 13 — the polar cap — is returned: the band cells' boxes stop short of the
pole even though their vertices reach latitude 80°, so the tree prunes them.

The one thing to keep straight is that a tree and its queries must live on the same
manifold.  `GI.extent(geom)` of a geographic geometry is a lon/lat box; handed to a
spherical tree, `Extents.intersects` will happily compare its `X` (longitude, up to 180)
against the tree's `X` (a Cartesian coordinate, up to 1) and return nonsense.  Convert
the query the same way the tree was built: `Extents.extent(Spherical(), geom)`.