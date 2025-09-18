# # Proj extension

# This is the main module for the GeometryOps extension on Proj.jl.
module GeometryOpsProjExt

using GeometryOps, Proj

include("reproject.jl")
include("segmentize.jl")
include("perimeter_area.jl")

end