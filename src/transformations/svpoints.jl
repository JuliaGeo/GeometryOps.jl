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
            GI.Point(SA.SVector{4}(T(GI.x(p)), T(GI.y(p)), T(GI.z(p)), T(GI.m(p))))
        end
    elseif _is3d(geom)
        return apply(PointTrait(), geom; kw...) do p
            GI.Point(SA.SVector{3}(T(GI.x(p)), T(GI.y(p)), T(GI.z(p))))
        end
    else
        return apply(PointTrait(), geom; kw...) do p
            GI.Point(SA.SVector{2}(T(GI.x(p)), T(GI.y(p))))
        end
    end
end
