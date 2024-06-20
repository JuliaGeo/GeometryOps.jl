# # SVector Points conversion

"""
    svpoints(obj)

Convert all points in `obj` to SVPoints, which are a subtype of a StaticVector, wherever
they are nested.

Returns a similar object or collection of objects using GeoInterface.jl geometries wrapping
SVPoints.

# Keywords

$APPLY_KEYWORDS
"""
function svpoints(geom, ::Type{T} = Float64; kw...) where T
    if _ismeasured(geom)
        return apply(PointTrait(), geom; kw...) do p
            SVPoint_4D(p)
        end
    elseif _is3d(geom)
        return apply(PointTrait(), geom; kw...) do p
            SVPoint_3D(p)
        end
    else
        return apply(PointTrait(), geom; kw...) do p
            SVPoint_2D(p)
        end
    end
end
