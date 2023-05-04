module GeometryOps

using GeoInterface
using GeometryBasics

const GI = GeoInterface

include("methods/signed_distance.jl")
include("methods/signed_area.jl")
include("methods/centroid.jl")
include("methods/contains.jl")

end
