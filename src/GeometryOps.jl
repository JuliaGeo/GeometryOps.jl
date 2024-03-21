# # GeometryOps.jl

module GeometryOps

using GeoInterface
using GeometryBasics
using GeometryBasics.StaticArrays
import Proj
using LinearAlgebra
import Proj.CoordinateTransformations.StaticArrays
import Base.@kwdef

using GeoInterface.Extents: Extents

const GI = GeoInterface
const GB = GeometryBasics

const TuplePoint{T} = Tuple{T, T} where T <: AbstractFloat
const Edge{T} = Tuple{TuplePoint{T},TuplePoint{T}} where T

include("primitives.jl")
include("utils.jl")

include("methods/angles.jl")
include("methods/area.jl")
include("methods/barycentric.jl")
include("methods/centroid.jl")
include("methods/distance.jl")
include("methods/equals.jl")
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
include("transformations/simplify.jl")
include("transformations/tuples.jl")
include("transformations/transform.jl")

end
