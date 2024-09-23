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

"""
    _wrap(geom; crs, calc_extent)

Wraps `geom` in a GI wrapper geometry of its geometry trait.  This allows us
to attach CRS and extent info to geometry types which otherwise could not hold
those, like LibGEOS and WKB geometries.

Returns a GI wrapper geometry, for which `parent(result) == geom`.
"""
function _wrap(geom; crs=GI.crs(geom), calc_extent = true)
    return GI.geointerface_geomtype(GI.geomtrait(geom))(geom; crs, extent = GI.extent(geom, calc_extent))
end

include("buffer.jl")
include("segmentize.jl")
include("simplify.jl")

include("simple_overrides.jl")

end
