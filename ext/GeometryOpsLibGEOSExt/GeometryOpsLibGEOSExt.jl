module GeometryOpsLibGEOSExt

import GeometryOps as GO, LibGEOS as LG
import GeometryOps: GI

import GeometryOps: GEOS, enforce

using GeometryOps
for name in filter(!in((:var"#eval", :eval, :var"#include", :include)), names(GeometryOps; all = true))
    @eval using GeometryOps: $name
end

include("buffer.jl")
include("segmentize.jl")
include("simple_overrides.jl")

end