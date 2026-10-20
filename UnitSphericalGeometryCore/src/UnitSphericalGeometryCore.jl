"""
    UnitSphericalGeometryCore

The essential types for geometry on the unit sphere, with as few dependencies as possible:

- [`UnitSphericalPoint`](@ref): a point on the unit 2-sphere, stored as Cartesian `xyz`
  and usable as a static 3-vector and as a GeoInterface point.
- [`UnitSphereFromGeographic`](@ref) and [`GeographicFromUnitSphere`](@ref): the
  transformations between (longitude, latitude) in degrees and unit-sphere `xyz`.
- [`SphericalCap`](@ref): a circular region of the sphere, with the Extents.jl predicates
  (`intersects`, `contains`, `within`, `union`, `grow`) and conversion to a Cartesian
  `Extent{(:X, :Y, :Z)}` for use in spatial indices.
- [`spherical_bounding_cap`](@ref) and [`to_unit_spherical_points`](@ref): go from any
  GeoInterface geometry to a bounding cap, or to its vertices as unit-spherical points.
- [`spherical_distance`](@ref), [`slerp`](@ref) and [`robust_cross_product`](@ref): the
  numerically careful primitives (ported from S2) that the above are built on.

Predicates (orientation, arc intersection) and algorithms live in packages that build on
this one, such as `GeometryOps.UnitSpherical`.
"""
module UnitSphericalGeometryCore

using CoordinateTransformations
using StaticArrays, LinearAlgebra
import GeoInterface as GI, GeoFormatTypes as GFT
import Extents
import Random

export UnitSphericalPoint, UnitSphereFromGeographic, GeographicFromUnitSphere,
       slerp, SphericalCap, spherical_distance,
       spherical_bounding_cap, to_unit_spherical_points,
       robust_cross_product

@static if VERSION >= v"1.11"
    eval(Meta.parse("public circumcenter_on_unit_sphere, angle_between, RobustCrossProduct"))
end

include("point.jl")

include("robustcrossproduct/RobustCrossProduct.jl")
using .RobustCrossProduct: robust_cross_product, min_stable_norm

include("coordinate_transforms.jl")
include("slerp.jl")
include("cap.jl")
include("bounding_cap.jl")

end
