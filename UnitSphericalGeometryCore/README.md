# UnitSphericalGeometryCore.jl

The essential types for geometry on the unit sphere, with as few dependencies as possible
(GeoInterface, Extents, CoordinateTransformations, StaticArrays and GeoFormatTypes).

- `UnitSphericalPoint`: a point on the unit 2-sphere as Cartesian `xyz`, usable as a static
  3-vector and as a GeoInterface point.
- `UnitSphereFromGeographic` / `GeographicFromUnitSphere`: CoordinateTransformations between
  (longitude, latitude) in degrees and unit-sphere `xyz`.
- `SphericalCap`: a circular region of the sphere, with the Extents.jl predicates
  (`intersects`, `contains`, `within`, `union`, `grow`) and conversion to a Cartesian
  `Extent{(:X, :Y, :Z)}` for spatial indexing.
- `spherical_bounding_cap(geom)` and `to_unit_spherical_points(geom)`: from any GeoInterface
  geometry to a bounding cap, or to its vertices on the unit sphere.
- `spherical_distance`, `slerp`, `robust_cross_product`: the numerically careful primitives
  (ported from Google's S2) everything else is built on.

```julia
using UnitSphericalGeometryCore
import GeoInterface as GI

p = UnitSphericalPoint(GI.Point(-74.0, 40.7))       # lon, lat in degrees -> unit xyz
GeographicFromUnitSphere()(p)                        # back to (lon, lat)

cap = SphericalCap(p, 0.01)                          # 0.01 rad around p
ring = GI.LineString([(-75.0, 40.0), (-73.0, 40.0), (-73.0, 41.0)])
Extents.intersects(cap, spherical_bounding_cap(ring))
```

This package is developed in the [GeometryOps.jl](https://github.com/JuliaGeo/GeometryOps.jl)
monorepo, where `GeometryOps.UnitSpherical` builds the orientation predicates and arc
intersection on top of it.
