# # GeometryOps.jl

module GeometryOps

import GeometryOpsCore
import GeometryOpsCore: 
                TraitTarget,
                Manifold, Planar, Spherical, Geodesic, AutoManifold, WrongManifoldException,
                Algorithm, AutoAlgorithm, ManifoldIndependentAlgorithm, SingleManifoldAlgorithm, NoAlgorithm,
                BoolsAsTypes, True, False, booltype,
                apply, applyreduce, 
                flatten, reconstruct, rebuild, unwrap, _linearring,
                APPLY_KEYWORDS, THREADED_KEYWORD, CRS_KEYWORD, CALC_EXTENT_KEYWORD

export TraitTarget, Manifold, Planar, Spherical, Geodesic, apply, applyreduce, flatten, reconstruct, rebuild, unwrap 

using GeoInterface
using GeometryBasics
using LinearAlgebra, Statistics

using GeometryBasics.StaticArrays

import Tables, DataAPI
import DelaunayTriangulation # for convex hull and triangulation
import ExactPredicates
import Base.@kwdef
import GeoInterface.Extents: Extents
import SortTileRecursiveTree
import SortTileRecursiveTree: STRtree

const GI = GeoInterface
const GB = GeometryBasics

const TuplePoint{T} = Tuple{T, T} where T <: AbstractFloat
const Edge{T} = Tuple{TuplePoint{T},TuplePoint{T}} where T

include("types.jl") # backend / algorithm types like GEOS, PROJ, etc.
include("primitives.jl") # moved to GeometryOpsCore
include("not_implemented_yet.jl") # functions that are not implemented yet, but we want stubs for, or the implementations might be in extensions

# Include utility modules first!
include("utils/LoopStateMachine/LoopStateMachine.jl") # Utils for functions that can tell the loop they run in to do something via the return value
include("utils/SpatialTreeInterface/SpatialTreeInterface.jl") # Utils for spatial trees
include("utils/NaturalIndexing.jl") # Utils for natural indexing
include("utils/UnitSpherical/UnitSpherical.jl") # Utils for unit spherical
include("utils/utils.jl") # More general utility functions

# Load utility modules in
using .NaturalIndexing, .SpatialTreeInterface, .LoopStateMachine, .UnitSpherical

include("methods/angles.jl")
include("methods/area.jl")
include("methods/barycentric.jl")
include("methods/buffer.jl")
include("methods/centroid.jl")
include("methods/convex_hull.jl")
include("methods/distance.jl")
include("methods/equals.jl")
include("methods/clipping/predicates.jl")
include("methods/clipping/clipping_processor.jl")
include("methods/clipping/coverage.jl")
include("methods/clipping/cut.jl")
include("methods/clipping/intersection.jl")
include("methods/clipping/difference.jl")
include("methods/clipping/union.jl")
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
include("methods/orientation.jl")
include("methods/polygonize.jl")

include("transformations/extent.jl")
include("transformations/flip.jl")
include("transformations/reproject.jl")
include("transformations/segmentize.jl")
include("transformations/simplify.jl")
include("transformations/tuples.jl")
include("transformations/transform.jl")
include("transformations/forcedims.jl")
include("transformations/correction/geometry_correction.jl")
include("transformations/correction/closed_ring.jl")
include("transformations/correction/intersecting_polygons.jl")

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

end
