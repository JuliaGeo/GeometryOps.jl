

"""
    reproject(geometry; source_crs, target_crs, transform, always_xy, time)
    reproject(geometry, source_crs, target_crs; always_xy, time)
    reproject(geometry, transform; always_xy, time)

Reproject any GeoInterface.jl compatible `geometry` from `source_crs` to `target_crs`.

The returned object will be constructed from `GeoInterface.WrapperGeometry`
geometries, wrapping views of a `Vector{Proj.Point{D}}`, where `D` is the dimension.

# Arguments

- `geometry`: Any GeoInterface.jl compatible geometries.
- `source_crs`: the source coordinate referece system, as a GeoFormatTypes.jl object or a string.
- `target_crs`: the target coordinate referece system, as a GeoFormatTypes.jl object or a string.

If these a passed as keywords, `transform` will take priority.
Without it `target_crs` is always needed, and `source_crs` is
needed if it is not retreivable from the geometry with `GeoInterface.crs(geometry)`.

# Keywords

-`always_xy`: force x, y coordinate order, `true` by default.
    `false` will expect and return points in the crs coordinate order.
-`time`: the time for the coordinates. `Inf` by default.
"""
function reproject(geom;
    source_crs=nothing, target_crs=nothing, transform=nothing, kw...
)
    if isnothing(transform)
        source_crs = isnothing(source_crs) ? GeoInterface.crs(geom) : source_crs
        isnothing(source_crs) && throw(ArgumentError("geom has no crs attatched. Pass a `source_crs` keyword"))
        reproject(geom, source_crs, target_crs; kw...)
    else
        reproject(geom, transform; kw...)
    end
end
function reproject(geom, source_crs, target_crs;
    time=Inf,
    always_xy=true,
    transform=Proj.Transformation(Proj.CRS(source_crs), Proj.CRS(target_crs); always_xy),
)
    reproject(geom, transform; time, target_crs)
end
function reproject(geom, transform::Proj.Transformation; time=Inf, target_crs=nothing)
    if _is3d(geom)
        return apply(PointTrait, geom; crs=target_crs) do p
            transform(GI.x(p), GI.y(p), GI.z(p))
        end
    else
        return apply(PointTrait, geom; crs=target_crs) do p
            transform(GI.x(p), GI.y(p))
        end
    end
end
