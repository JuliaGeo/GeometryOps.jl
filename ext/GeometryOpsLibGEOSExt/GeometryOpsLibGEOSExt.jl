module GeometryOpsLibGEOSExt

import GeometryOps as GO, LibGEOS as LG
import GeometryOps: GI

using GeometryOps
for name in names(GeometryOps; all = true)
    @eval using GeometryOps: $name
end

"""
    enforce(alg::GO.GEOS, kw::Symbol, f)

Enforce the presence of a keyword argument in a `GEOS` algorithm, and return `alg.params[kw]`.

Throws an error if the key is not present, and mentions `f` in the error message (since there isn't 
a good way to get the name of the function that called this method).
"""
function enforce(alg::GO.GEOS, kw::Symbol, f)
    if haskey(alg.params, kw)
        return alg.params[kw]
    else
        error("$(f) requires a `$(kw)` keyword argument to the `GEOS` algorithm, which was not provided.")
    end
end

include("buffer.jl")
include("segmentize.jl")
include("simple_overrides.jl")

end