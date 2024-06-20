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
    if _ismeasured(geom)
        return apply(PointTrait(), geom; kw...) do p
            TuplePoint_4D(p, T)
        end
    elseif _is3d(geom)
        return apply(PointTrait(), geom; kw...) do p
            TuplePoint_3D(p, T)
        end
    else
        return apply(PointTrait(), geom; kw...) do p
            TuplePoint_2D(p, T)
        end
    end
end
