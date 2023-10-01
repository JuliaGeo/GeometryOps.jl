# # Tuple conversion

"""
    tuples(obj)

Convert all points on obj to `Tuple`s.
"""
function tuples(geom; kw...) 
    if _is3d(geom)
        return apply(PointTrait, geom; kw...) do p
            (Float64(GI.x(p)), Float64(GI.y(p)), Float64(GI.z(p)))
        end
    else
        return apply(PointTrait, geom; kw...) do p
            (Float64(GI.x(p)), Float64(GI.y(p)))
        end
    end
end
