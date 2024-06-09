module GeometryOpsLibGEOSExt

import GeometryOps as GO, LibGEOS as LG
import GeometryOps: GI

import GeometryOps: GEOS, enforce

using GeometryOps
# The filter statement is required because in Julia, each module has its own versions of these
# functions, which serve to evaluate or include code inside the scope of the module.
# However, if you import those from another module (which you would with `all=true`),
# that creates an ambiguity which causes a warning during precompile/load time.
# In order to avoid this, we filter out these special functions.
for name in filter(!in((:var"#eval", :eval, :var"#include", :include)), names(GeometryOps; all = true))
    @eval using GeometryOps: $name
end

include("buffer.jl")
include("segmentize.jl")
include("simplify.jl")
include("simple_overrides.jl")

end