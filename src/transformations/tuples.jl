"""
    tuples(obj)

Convert all points on obj to `Tuple`s.
"""
function tuples(geom) 
    if _is3d(geom)
        return apply(PointTrait, geom) do p
            (GI.x(p), GI.y(p), GI.z(p))
        end
    else
        return apply(PointTrait, geom) do p
            (GI.x(p), GI.y(p))
        end
    end
end
