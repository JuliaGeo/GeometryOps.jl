# # Coordinate flipping

# This is a simple example of how to use the `apply` functionality in a function,
# by flipping the x and y coordinates of a geometry.

"""
    flip(obj)

Swap all of the x and y coordinates in obj, otherwise
keeping the original structure (but not necessarily the
original type).

## Keywords 

$APPLY_KEYWORDS
"""
function flip(geom; kw...)
    if _ismeasured(geom)
        return apply(PointTrait(), geom; kw...) do p
            SVPoint_4D((GI.y(p), GI.x(p), GI.z(p), GI.m(p)))
        end
    elseif _is3d(geom)
        return apply(PointTrait(), geom; kw...) do p
            SVPoint_3D((GI.y(p), GI.x(p), GI.z(p)))
        end
    else
        return apply(PointTrait(), geom; kw...) do p
            SVPoint_2D((GI.y(p), GI.x(p)))
        end
    end
end
