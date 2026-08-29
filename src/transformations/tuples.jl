# # Tuple conversion

"""
    tuples(obj)

Convert all points in `obj` to `Tuple`s, wherever the are nested.

Returns a similar object or collection of objects using GeoInterface.jl
geometries wrapping `Tuple` points.

# Keywords

- `threaded`: `true` or `false`. Whether to use multithreading. Defaults to `false`.
- `crs`: The CRS to attach to geometries. Defaults to `nothing`.
- `calc_extent`: `true` or `false`. Whether to calculate the extent. Defaults to `true`.
"""
function tuples(geom, ::Type{T} = Float64; calc_extent = true, kw...) where T
    if _is3d(geom)
        return apply(PointTrait(), geom; calc_extent, kw...) do p
            (T(GI.x(p)), T(GI.y(p)), T(GI.z(p)))
        end
    else
        return apply(PointTrait(), geom; calc_extent, kw...) do p
            (T(GI.x(p)), T(GI.y(p)))
        end
    end
end
