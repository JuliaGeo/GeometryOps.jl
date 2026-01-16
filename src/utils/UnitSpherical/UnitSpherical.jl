module UnitSpherical

using CoordinateTransformations
using StaticArrays, LinearAlgebra
import GeoInterface as GI, GeoFormatTypes as GFT

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

export UnitSphericalPoint, UnitSphereFromGeographic, GeographicFromUnitSphere,
       slerp, SphericalCap, spherical_distance, spherical_orient, point_on_spherical_arc,
       spherical_arc_intersection, ArcIntersectionResult,
       arc_cross, arc_hinge, arc_overlap, arc_disjoint,
       to_unit_spherical_points

"""
    to_unit_spherical_points(ring) -> Vector{UnitSphericalPoint{Float64}}

Convert a ring (linear ring or any GeoInterface point iterator) to a vector of UnitSphericalPoints.
Uses UnitSphereFromGeographic which is a no-op for already-converted points.
"""
function to_unit_spherical_points(ring)
    transform = UnitSphereFromGeographic()
    return [transform((GI.x(p), GI.y(p))) for p in GI.getpoint(ring)]
end

end