module GeometryOps

using GeoInterface
using GeometryBasics
import Proj

const GI = GeoInterface

include("primitives.jl")
include("methods/signed_distance.jl")
include("methods/signed_area.jl")
include("methods/centroid.jl")
include("methods/contains.jl")
include("transformations/reproject.jl")
include("transformations/flip.jl")

end
