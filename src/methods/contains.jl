# More GeometryBasics code

_cross(p1, p2, p3) = (p1[1] - p3[1]) * (p2[2] - p3[2]) - (p2[1] - p3[1]) * (p1[2] - p3[2])

# Implementation of a point-in-polygon algorithm
# from Luxor.jl.  This is the Hormann-Agathos (2001) algorithm.

# For the source, see https://github.com/JuliaGraphics/Luxor.jl/blob/66d60fb51f6b1bb38690fe8dcc6c0084eeb80710/src/polygons.jl#L190-L229.
function contains(ls::GeometryBasics.LineString{2, T1}, point::Point{2, T2}) where {T1, T2}
    pointlist = decompose(Point{2, promote_type(T1, T2)}, ls)
    c = false
    @inbounds for counter in eachindex(pointlist)
        q1 = pointlist[counter]
        # if reached last point, set "next point" to first point
        if counter == length(pointlist)
            q2 = pointlist[1]
        else
            q2 = pointlist[counter + 1]
        end
        if q1 == point
            # allowonedge || error("isinside(): VertexException a")
            continue
        end
        if q2[2] == point[2]
            if q2[1] == point[1]
                # allowonedge || error("isinside(): VertexException b")
                continue
            elseif (q1[2] == point[2]) && ((q2[1] > point[1]) == (q1[1] < point[1]))
                # allowonedge || error("isinside(): EdgeException")
                continue
            end
        end
        if (q1[2] < point[2]) != (q2[2] < point[2]) # crossing
            if q1[1] >= point[1]
                if q2[1] > point[1]
                    c = !c
                elseif ((_cross(q1, q2, point) > 0) == (q2[2] > q1[2]))
                    c = !c
                end
            elseif q2[1] > point[1]
                if ((_cross(q1, q2, point) > 0) == (q2[2] > q1[2]))
                    c = !c
                end
            end
        end
    end
    return c

end

function contains(poly::Polygon{2, T1}, point::Point{2, T2}) where {T1, T2}
    c = contains(poly.exterior, point)
    for interior in poly.interiors
        if contains(interior, point)
            return false
        end
    end
    return c
end

# TODOs: implement contains for mesh, simplex, and 3d objects (eg rect, triangle, etc.)

contains(mp::MultiPolygon{2, T1}, point::Point{2, T2}) where {T1, T2} = any((contains(poly, point) for poly in mp.polygons))
