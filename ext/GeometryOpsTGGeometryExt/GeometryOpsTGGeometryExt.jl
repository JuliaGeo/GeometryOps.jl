module GeometryOpsTGGeometryExt

using GeometryOps: TG
import GeometryOps as GO

using TGGeometry

for jl_fname in TGGeometry.TG_PREDICATES
    @eval GO.$jl_fname(::TG, geom1, geom2) = TGGeometry.$jl_fname(geom1, geom2)
end

end