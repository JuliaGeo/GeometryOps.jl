import GeometryOps: GI, GeoInterface, reproject, apply, transform, _is3d, True, False
import Proj

function reproject(geom;
    source_crs=nothing, target_crs=nothing, transform=nothing, kw...
)
    if isnothing(transform)
        if isnothing(source_crs) 
            source_crs = if GI.trait(geom) isa Nothing && geom isa AbstractArray
                GeoInterface.crs(first(geom))
            else
                GeoInterface.crs(geom)
            end
        end

        # If its still nothing, error
        isnothing(source_crs) && throw(ArgumentError("geom has no crs attached. Pass a `source_crs` keyword"))

        # Otherwise reproject
        reproject(geom, source_crs, target_crs; kw...)
    else
        reproject(geom, transform; kw...)
    end
end
function reproject(geom, source_crs, target_crs;
    time=Inf,
    always_xy=true,
    transform=nothing,
    kw...
)
    transform = if isnothing(transform) 
        s = source_crs isa Proj.CRS ? source_crs : convert(String, source_crs)
        t = target_crs isa Proj.CRS ? target_crs : convert(String, target_crs)
        Proj.Transformation(s, t; always_xy)
    else
        transform
    end
    reproject(geom, transform; time, target_crs, kw...)
end
function reproject(geom, transform::Proj.Transformation; time=Inf, target_crs=nothing, kw...)
    if _is3d(geom)
        return apply(GI.PointTrait(), geom; crs=target_crs, kw...) do p
            transform(GI.x(p), GI.y(p), GI.z(p))
        end
    else
        return apply(GI.PointTrait(), geom; crs=target_crs, kw...) do p
            transform(GI.x(p), GI.y(p))
        end
    end
end
