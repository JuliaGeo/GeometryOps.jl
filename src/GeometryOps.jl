# # GeometryOps.jl

module GeometryOps

using GeoInterface
using GeometryBasics
import Proj
using LinearAlgebra
import ExactPredicates
import Proj.CoordinateTransformations.StaticArrays

using GeoInterface.Extents: Extents

const GI = GeoInterface
const GB = GeometryBasics

const TuplePoint = Tuple{Float64,Float64}
const Edge = Tuple{TuplePoint,TuplePoint}

include("primitives.jl")
include("utils.jl")

include("methods/area.jl")
include("methods/barycentric.jl")
include("methods/bools.jl")
include("methods/centroid.jl")
include("methods/distance.jl")
include("methods/equals.jl")
include("methods/geom_relations/contains.jl")
include("methods/geom_relations/coveredby.jl")
include("methods/geom_relations/crosses.jl")
include("methods/geom_relations/disjoint.jl")
include("methods/geom_relations/geom_geom_processors.jl")
include("methods/geom_relations/intersects.jl")
include("methods/geom_relations/overlaps.jl")
include("methods/geom_relations/within.jl")
include("methods/polygonize.jl")

include("transformations/extent.jl")
include("transformations/flip.jl")
include("transformations/reproject.jl")
include("transformations/simplify.jl")
include("transformations/tuples.jl")
include("transformations/transform.jl")

end
