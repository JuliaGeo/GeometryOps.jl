# # Tuple conversion

"""
    tuples(obj)

Convert all points in `obj` to `Tuple`s, wherever the are nested.

Returns a similar object or collection of objects using GeoInterface.jl
geometries wrapping `Tuple` points.

# Keywords

$APPLY_KEYWORDS
"""
function tuples(geom, ::Type{T} = Float64; kw...) where T
    if _is3d(geom)
        return apply(PointTrait(), geom; kw...) do p
            (T(GI.x(p)), T(GI.y(p)), T(GI.z(p)))
        end
    else
        return apply(PointTrait(), geom; kw...) do p
            (T(GI.x(p)), T(GI.y(p)))
        end
    end
end
