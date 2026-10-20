"""
    UnitSpherical

Geometry on the unit sphere.  The types and primitives (`UnitSphericalPoint`,
`SphericalCap`, the geographic transformations, `slerp`, `spherical_distance`,
`robust_cross_product`, bounding caps) come from the UnitSphericalGeometryCore sublibrary
and are re-exported here; this module adds the orientation predicates, great-circle arc
intersection and arc extents on top of them.
"""
module UnitSpherical

using UnitSphericalGeometryCore
# Internals that GeometryOps and its tests reach through `UnitSpherical.`
using UnitSphericalGeometryCore: circumcenter_on_unit_sphere, angle_between,
    _intersects, _disjoint, _contains, _unit_spherical_point
using UnitSphericalGeometryCore.RobustCrossProduct: min_stable_norm
const RobustCrossProduct = UnitSphericalGeometryCore.RobustCrossProduct

using StaticArrays, LinearAlgebra
import GeoInterface as GI
import Extents
import ExactPredicates

include("predicates.jl")
include("arc_intersection.jl")
include("arc_extent.jl")

# Re-export the sublibrary's API, so `using GeometryOps.UnitSpherical` is unchanged
export UnitSphericalPoint, UnitSphereFromGeographic, GeographicFromUnitSphere,
       slerp, SphericalCap, spherical_distance, spherical_bounding_cap,
       to_unit_spherical_points, robust_cross_product
# This module's own API
export spherical_orient, exact_spherical_orient,
       point_on_spherical_arc,
       spherical_ring_contains, spherical_ring_encloses, spherical_exterior_anchor,
       spherical_arc_intersection, ArcIntersectionResult,
       arc_cross, arc_hinge, arc_overlap, arc_disjoint,
       spherical_arc_extent

end
