# # Tuple conversion

"""
    tuples(obj)

Convert all points in `obj` to `Tuple`s, wherever the are nested.

Returns a similar object or collection of objects using GeoInterface.jl
geometries wrapping `Tuple` points.

# Keywords

$APPLY_KEYWORDS
"""
function tuples(geom; kw...) 
    if _is3d(geom)
        return apply(PointTrait(), geom; kw...) do p
            (Float64(GI.x(p)), Float64(GI.y(p)), Float64(GI.z(p)))
        end
    else
        return apply(PointTrait(), geom; kw...) do p
            (Float64(GI.x(p)), Float64(GI.y(p)))
        end
    end
end
