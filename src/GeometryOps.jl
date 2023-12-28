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

include("methods/bools.jl")
include("methods/distance.jl")
include("methods/area.jl")
include("methods/centroid.jl")
include("methods/intersects.jl")
include("methods/contains.jl")
include("methods/crosses.jl")
include("methods/disjoint.jl")
include("methods/overlaps.jl")
include("methods/within.jl")
include("methods/polygonize.jl")
include("methods/barycentric.jl")
include("methods/equals.jl")

include("transformations/extent.jl")
include("transformations/flip.jl")
include("transformations/simplify.jl")
include("transformations/reproject.jl")
include("transformations/tuples.jl")
include("transformations/transform.jl")

end
