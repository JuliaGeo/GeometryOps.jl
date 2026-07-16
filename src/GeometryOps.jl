# # GeometryOps.jl

module GeometryOps

import GeometryOpsCore
import GeometryOpsCore:
                TraitTarget,
                Manifold, Planar, Spherical, Geodesic, AutoManifold, WrongManifoldException,
                manifold, best_manifold,
                Algorithm, AutoAlgorithm, ManifoldIndependentAlgorithm, SingleManifoldAlgorithm, NoAlgorithm,
                BoolsAsTypes, True, False, booltype, istrue,
                TaskFunctors,
                WithTrait,
                WithXY, WithXYZ, WithXYM, WithXYZM,
                apply, applyreduce,
                flatten, reconstruct, rebuild, unwrap, _linearring,
                get_geometries,
                APPLY_KEYWORDS, THREADED_KEYWORD, CRS_KEYWORD, CALC_EXTENT_KEYWORD

export TraitTarget, Manifold, Planar, Spherical, Geodesic, apply, applyreduce, flatten, reconstruct, rebuild, unwrap, get_geometries 

using GeoInterface
using LinearAlgebra, Statistics, Random

using StaticArrays

import Tables, DataAPI
import StaticArrays
import DelaunayTriangulation # for convex hull and triangulation
import ExactPredicates
import Base.@kwdef
import GeoInterface.Extents: Extents
import SortTileRecursiveTree
import SortTileRecursiveTree: STRtree

const GI = GeoInterface
const DelTri = DelaunayTriangulation

const TuplePoint{T} = Tuple{T, T} where T <: AbstractFloat
const Edge{T} = Tuple{TuplePoint{T},TuplePoint{T}} where T

include("types.jl")
include("primitives.jl")
include("not_implemented_yet.jl")

include("utils/LoopStateMachine/LoopStateMachine.jl")
include("utils/SpatialTreeInterface/SpatialTreeInterface.jl")
include("utils/UnitSpherical/UnitSpherical.jl")

# Load utility modules in
using .LoopStateMachine, .SpatialTreeInterface, .UnitSpherical

include("utils/utils.jl")

include("utils/NaturalIndexing.jl")
using .NaturalIndexing

include("utils/FlexibleRTrees/FlexibleRTrees.jl")
using .FlexibleRTrees

# Load utility modules in
using .NaturalIndexing, .SpatialTreeInterface, .LoopStateMachine

include("methods/angles.jl")
include("methods/area.jl")
include("methods/barycentric.jl")
include("methods/buffer.jl")
include("methods/centroid.jl")
include("methods/convex_hull.jl")
include("methods/distance.jl")
include("methods/equals.jl")
include("methods/extent.jl")
include("methods/perimeter.jl")
include("methods/clipping/predicates.jl")
include("methods/clipping/clipping_processor.jl")
include("methods/clipping/coverage.jl")
include("methods/clipping/cut.jl")
include("methods/clipping/intersection.jl")
include("methods/clipping/difference.jl")
include("methods/clipping/union.jl")
include("methods/clipping/sutherland_hodgman.jl")
include("methods/geom_relations/contains.jl")
include("methods/geom_relations/coveredby.jl")
include("methods/geom_relations/covers.jl")
include("methods/geom_relations/crosses.jl")
include("methods/geom_relations/disjoint.jl")
include("methods/geom_relations/geom_geom_processors.jl")
include("methods/geom_relations/intersects.jl")
include("methods/geom_relations/overlaps.jl")
include("methods/geom_relations/touches.jl")
include("methods/geom_relations/within.jl")
include("methods/geom_relations/common.jl")
include("methods/geom_relations/relateng/de9im.jl")
include("methods/geom_relations/relateng/topology_predicate.jl")
include("methods/geom_relations/relateng/relate_predicates.jl")
# Kernel files: after de9im.jl (they use the `LOC_` constants) and before the
# topology-layer files (Task 6+) that will call the kernel functions.
include("methods/geom_relations/relateng/kernel.jl")
include("methods/geom_relations/relateng/kernel_planar.jl")
include("methods/geom_relations/relateng/kernel_spherical.jl")
# Node sections: after the kernel (uses `NodeKey`), before the point locator
# (AdjacentEdgeLocator builds NodeSections).
include("methods/geom_relations/relateng/node_sections.jl")
# Polygon node converter: after node sections (rewrites a polygon's
# NodeSection group into maximal-ring structure for `create_node`).
include("methods/geom_relations/relateng/polygon_node_converter.jl")
# Node-edge topology: after node sections and the converter (`create_node`
# assembles a RelateNode from converted sections).
include("methods/geom_relations/relateng/relate_node.jl")
# Indexed point-in-area location: after the kernel (uses `rk_orient`,
# `_node_points` and the `LOC_` codes), before the point locator (a prepared
# RelatePointLocator caches these locators per polygonal element).
include("methods/geom_relations/relateng/indexed_point_in_area.jl")
# Point location: after the kernel (uses `_node_point` and de9im constants).
include("methods/geom_relations/relateng/point_locator.jl")
# Input facade: after the point locator (RelateGeometry wraps a lazy
# RelatePointLocator) and node sections (RelateSegmentString creates them).
include("methods/geom_relations/relateng/relate_geometry.jl")
# Topology computer: after the input facade and node topology (it drives
# RelateGeometry locates, NodeSections grouping and RelateNode evaluation).
include("methods/geom_relations/relateng/topology_computer.jl")
# Edge intersector: feeds segment-pair intersections into the computer.
include("methods/geom_relations/relateng/edge_intersector.jl")
# The RelateNG engine: drives all of the above through the phased evaluation.
include("methods/geom_relations/relateng/relate_ng.jl")

# OverlayNG noding substrate (phase 1): geometries → exactly-noded arrangement.
include("methods/clipping/overlayng/noding/noded_arrangement.jl")
include("methods/clipping/overlayng/noding/collect.jl")
include("methods/clipping/overlayng/noding/node_identity.jl")
include("methods/clipping/overlayng/noding/split.jl")
include("methods/clipping/overlayng/noding/emit.jl")

# OverlayNG engine core (phase 2a): the half-edge graph over the arrangement.
include("methods/clipping/overlayng/overlay_label.jl")
include("methods/clipping/overlayng/edge_source.jl")
include("methods/clipping/overlayng/half_edge.jl")
include("methods/clipping/overlayng/overlay_graph.jl")

# OverlayNG engine core (phase 2b): labeller, result builders, and the driver.
include("methods/clipping/overlayng/overlay_labeller.jl")
include("methods/clipping/overlayng/maximal_edge_ring.jl")
include("methods/clipping/overlayng/polygon_builder.jl")
include("methods/clipping/overlayng/line_builder.jl")
include("methods/clipping/overlayng/intersection_point_builder.jl")
include("methods/clipping/overlayng/overlay_ng.jl")

include("methods/orientation.jl")
include("methods/polygonize.jl")
include("methods/minimum_bounding_circle.jl")
include("methods/voronoi.jl")

include("transformations/extent.jl")
include("transformations/flip.jl")
include("transformations/reproject.jl")
include("transformations/segmentize.jl")
include("transformations/simplify.jl")
include("transformations/smooth.jl")
include("transformations/tuples.jl")
include("transformations/transform.jl")
include("transformations/forcedims.jl")
include("transformations/correction/geometry_correction.jl")
include("transformations/correction/closed_ring.jl")
include("transformations/correction/intersecting_polygons.jl")
include("transformations/correction/antipodal_edge_split.jl")
include("transformations/correction/crossing_edge_split.jl")

# Import all names from GeoInterface and Extents, so users can do `GO.extent` or `GO.trait`.
for name in names(GeoInterface)
    @eval using GeoInterface: $name
end
for name in names(Extents)
    @eval using GeoInterface.Extents: $name
end

function __init__()
    # Handle all available errors!
    Base.Experimental.register_error_hint(_reproject_error_hinter, MethodError)
    Base.Experimental.register_error_hint(_geodesic_segments_error_hinter, MethodError)
    Base.Experimental.register_error_hint(_buffer_error_hinter, MethodError)
end

include("precompile.jl")

end
