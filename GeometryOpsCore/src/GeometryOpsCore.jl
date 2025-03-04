module GeometryOpsCore

using Base.Threads: nthreads, @threads, @spawn

import GeoInterface
import GeoInterface as GI
import GeoInterface: Extents

# Import all names from GeoInterface and Extents, so users can do `GO.extent` or `GO.trait`.
for name in names(GeoInterface)
    @eval using GeoInterface: $name
end
for name in names(Extents)
    @eval using GeoInterface.Extents: $name
end

using Tables
using DataAPI
import StableTasks

include("keyword_docs.jl")
include("types.jl")

include("apply.jl")
include("applyreduce.jl")
include("other_primitives.jl")
include("geometry_utils.jl")

end