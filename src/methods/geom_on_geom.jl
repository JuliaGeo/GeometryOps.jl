export point_on_line, line_on_line

"""
    point_on_line(point::Point, line::LineString; ignore_end_vertices::Bool=false)::Bool

Return true if a point is on a line. Accept a optional parameter to ignore the
start and end vertices of the linestring.

## Examples

```jldoctest
import GeoInterface as GI, GeometryOps as GO

point = (1, 1)
line = GI.LineString([(0, 0), (3, 3), (4, 4)])
GO.point_on_line(point, line)

# output
true
```
"""
function point_on_line(point, line; ignore_end_vertices::Bool=false)::Bool
    line_points = tuple_points(line)
    n = length(line_points)

    exclude_boundary = :none
    for i in 1:n - 1
        if ignore_end_vertices
            if i === 1 
                exclude_boundary = :start
            elseif i === n - 2 
                exclude_boundary = :end
            elseif (i === 1 && i + 1 === n - 1) 
                exclude_boundary = :both
            end
        end
        if point_on_segment(point, line_points[i], line_points[i + 1]) 
            return true
        end
    end
    return false
end

function point_on_segment(point, start, stop)
    # Parse out points
    x, y = GI.x(point), GI.y(point)
    x1, y1 = GI.x(start), GI.y(start)
    x2, y2 = GI.x(stop), GI.y(stop)
    Δx = x2 - x1
    Δy = y2 - y1
    #=
    Determine if the point is on the segment -> see if cross product of line and
    vector from line start to point is zero -> vectors are parallel. Then, check
    point is between segment endpoints. 
    =#
    on_line = _isparallel(Δx, Δy, (x - x1), (y - y1))
    between_endpoints = (x2 > x1 ? x1 <= x <= x2 : x2 <= x <= x1) &&
        (y2 > y1 ? y1 <= y <= y2 : y2 <= y <= y1)
    return on_line && between_endpoints
end


line_on_line(line1, line2) = line_on_line(trait(line1), line1, trait(line2), line2)
function line_on_line(t1::GI.AbstractCurveTrait, line1, t2::AbstractCurveTrait, line2)
    for p in GI.getpoint(line1)
        # FIXME: all points being on the line doesn't
        # actually mean the whole line is on the line...
        point_on_line(p, line2) || return false
    end
    return true
end
