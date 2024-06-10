# # SVector Points conversion

"""
    svpoints(obj)

Convert all points in `obj` to GI.Points wrapping StaticVectors, wherever the are nested.

Returns a similar object or collection of objects using GeoInterface.jl geometries wrapping
GI.Points containing a StaticVector.

# Keywords

$APPLY_KEYWORDS
"""
function svpoints(geom, ::Type{T} = Float64; kw...) where T
    if _ismeasured(geom)
        return apply(PointTrait(), geom; kw...) do p
            SVPoint{4, T, _True, _True}((GI.x(p), GI.y(p), GI.z(p), GI.m(p)))
        end
    elseif _is3d(geom)
        return apply(PointTrait(), geom; kw...) do p
            SVPoint{3, T, _True, _False}((GI.x(p), GI.y(p), GI.z(p)))
        end
    else
        return apply(PointTrait(), geom; kw...) do p
            SVPoint{2, T, _False, _False}((GI.x(p), GI.y(p)))
        end
    end
end
