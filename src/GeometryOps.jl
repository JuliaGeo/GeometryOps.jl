# # GeometryOps.jl

module GeometryOps

using GeoInterface
using GeometryBasics
import Proj
using LinearAlgebra

using GeoInterface.Extents: Extents

const GI = GeoInterface
const GB = GeometryBasics

include("primitives.jl")
include("utils.jl")

include("methods/bools.jl")
include("methods/signed_distance.jl")
include("methods/signed_area.jl")
include("methods/centroid.jl")
include("methods/intersects.jl")
include("methods/contains.jl")
include("methods/crosses.jl")
include("methods/disjoint.jl")
include("methods/overlaps.jl")
include("methods/within.jl")
include("methods/polygonize.jl")
include("methods/barycentric.jl")

include("transformations/flip.jl")
include("transformations/simplify.jl")
include("transformations/reproject.jl")

end
