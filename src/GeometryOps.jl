# # GeometryOps.jl

module GeometryOps

using GeoInterface
using GeometryBasics
import Proj
using LinearAlgebra

const GI = GeoInterface
const GB = GeometryBasics

include("primitives.jl")
include("utils.jl")

include("methods/signed_distance.jl")
include("methods/signed_area.jl")
include("methods/centroid.jl")
include("methods/contains.jl")
include("methods/polygonize.jl")
include("methods/barycentric.jl")

include("transformations/flip.jl")
include("transformations/simplify.jl")
include("transformations/reproject.jl")

end
