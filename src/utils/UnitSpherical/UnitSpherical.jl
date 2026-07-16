module UnitSpherical

using CoordinateTransformations
using StaticArrays, LinearAlgebra
import GeoInterface as GI, GeoFormatTypes as GFT
import Extents

import Random

# using TestItems # this is a thin package that allows TestItems.@testitem to be parsed.

include("point.jl")

include("robustcrossproduct/RobustCrossProduct.jl")
# Re-export from RobustCrossProduct
using .RobustCrossProduct: robust_cross_product
export robust_cross_product

include("coordinate_transforms.jl")
include("slerp.jl")
include("cap.jl")
include("predicates.jl")
include("arc_intersection.jl")
include("arc_extent.jl")

export UnitSphericalPoint, UnitSphereFromGeographic, GeographicFromUnitSphere,
       slerp, SphericalCap, spherical_distance, spherical_orient, point_on_spherical_arc,
       spherical_ring_contains, spherical_ring_encloses, spherical_exterior_anchor,
       spherical_arc_intersection, ArcIntersectionResult,
       arc_cross, arc_hinge, arc_overlap, arc_disjoint,
       spherical_arc_extent,
       to_unit_spherical_points

"""
    to_unit_spherical_points(ring) -> Vector{<:UnitSphericalPoint}

Convert a ring (linear ring or any GeoInterface point iterator) to a vector of
UnitSphericalPoints, treating geographic input as (longitude, latitude).
`UnitSphericalPoint`s pass through unchanged.
"""
function to_unit_spherical_points(ring)
    return [UnitSphericalPoint(p) for p in GI.getpoint(ring)]
end

end