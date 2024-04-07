module GeometryOpsLibGEOSExt

import GeometryOps as GO, LibGEOS as LG
import GeometryOps: GI

using GeometryOps
for name in names(GeometryOps; all = true)
    @eval using GeometryOps: $name
end

include("buffer.jl")
include("segmentize.jl")
include("simple_overrides.jl")

end