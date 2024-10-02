module GeometryOpsCore

using Base.Threads: nthreads, @threads, @spawn

import GeoInterface as GI

using Tables
using DataAPI

include("keyword_docs.jl")
include("types.jl")

include("apply.jl")
include("applyreduce.jl")
include("other_primitives.jl")

end