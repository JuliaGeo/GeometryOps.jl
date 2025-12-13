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
       arc_cross, arc_hinge, arc_overlap, arc_disjoint

end